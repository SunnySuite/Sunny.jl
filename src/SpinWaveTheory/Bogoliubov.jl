import LinearAlgebra: BlasInt

# Storage for `bogoliubov!`, holding the transformation `T` and energies `ε`
# that it produces, together with the LAPACK scratch space that `eigen!` would
# otherwise allocate on every call. The precision is `FT`.
struct BogoliubovWorkspace{FT}
    T     :: Matrix{Complex{FT}}   # Transformation, as documented by `excitations!`
    ε     :: Vector{FT}            # Signed energies, in the same grouping
    Tscr  :: Matrix{Complex{FT}}   # Scratch for the reordering below
    εscr  :: Vector{FT}
    work  :: Vector{Complex{FT}}   # LAPACK scratch for zhegvd or chegvd
    rwork :: Vector{FT}
    iwork :: Vector{BlasInt}
    info  :: Base.RefValue{BlasInt}
end

function BogoliubovWorkspace(L::Int, ::Type{FT}=Float64) where FT
    n = 2L
    CT = Complex{FT}
    # Minimum scratch sizes documented for ZHEGVD with jobz = 'V' and n > 1
    # (always true here, since n = 2L).
    return BogoliubovWorkspace(zeros(CT, n, n), zeros(FT, n), zeros(CT, n, n), zeros(FT, n),
                               zeros(CT, 2n + n^2), zeros(FT, 1 + 5n + 2n^2), zeros(BlasInt, 3 + 5n),
                               Ref(BlasInt(0)))
end

# The generalized eigenproblem A t = λ B t, overwriting A with the eigenvectors
# t and B with its Cholesky factor. This is `eigen!(Hermitian(A),
# Hermitian(B))` with every allocation hoisted into `ws`.
for (hegvd, FT) in ((:zhegvd_, :Float64), (:chegvd_, :Float32))
    @eval function hegvd!(ws::BogoliubovWorkspace{$FT}, A, B)
        (; ε, work, rwork, iwork, info) = ws
        n = size(A, 1)
        ccall((BLAS.@blasfunc($hegvd), LinearAlgebra.libblastrampoline), Cvoid,
              (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
               Ptr{Complex{$FT}}, Ref{BlasInt}, Ptr{Complex{$FT}}, Ref{BlasInt},
               Ptr{$FT}, Ptr{Complex{$FT}}, Ref{BlasInt}, Ptr{$FT},
               Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
               Clong, Clong),
              1, 'V', 'U', n, A, n, B, n, ε, work, length(work), rwork,
              length(rwork), iwork, length(iwork), info, 1, 1)
        info[] < 0 && throw(ArgumentError("Invalid argument #$(-info[]) to LAPACK $($(string(hegvd)))"))
        info[] > 0 && throw(PosDefException(info[]))
    end
end

# Bogoliubov transformation that diagonalizes a quadratic bosonic Hamiltonian,
# allowing for anomalous terms. The general procedure derives from Colpa,
# Physica A, 93A, 327-353 (1978). Overwrites data in H, if it has the precision
# of `ws`, and otherwise rounds a copy. The returned energies alias `ws.ε`, and
# so are valid only until the next call on the same workspace.
function bogoliubov!(ws::BogoliubovWorkspace{FT}, H::Matrix{<:Complex}) where FT
    (; T, ε, Tscr, εscr) = ws
    n = size(T, 1)
    L = div(n, 2)
    @assert size(H) == (n, n)
    # `eigen!` rejects a non-finite matrix before reaching LAPACK, which is not
    # guaranteed to terminate on one. Retain that check.
    all(isfinite, H) || throw(ArgumentError("Hamiltonian contains non-finite entries"))

    # Initialize T to the para-unitary identity Ĩ = diagm([ones(L), -ones(L)])
    T .= 0
    for i in 1:L
        T[i, i] = 1
        T[i+L, i+L] = -1
    end

    # Solve generalized eigenvalue problem, Ĩ t = λ H t, for columns t of T. A
    # copy of H in the precision of `ws` can use `Tscr`, which is needed only
    # after the solve.
    hegvd!(ws, T, H isa Matrix{Complex{FT}} ? H : copyto!(Tscr, H))

    # By Sylvester's theorem, "inertia" (sign signature) is invariant under a
    # congruence transform Ĩ → √H Ĩ √H, so exactly L of the λ are negative. LAPACK
    # returns them in ascending order, whereas they are wanted with the positive values
    # first and otherwise ascending in absolute value. That reordering is therefore the
    # fixed permutation [L+1:2L; L:-1:1], and no sort is required. It is not a product
    # of disjoint transpositions, so both T and ε are permuted out of scratch copies.
    #
    # Degenerate λ are the one case where this differs from sorting: `eigen!` breaks ties
    # with an unstable QuickSort, so the order within a degenerate block used to be
    # arbitrary, whereas it is now inherited from LAPACK. Either choice is a valid
    # eigenbasis, differing by a rotation within the block, and the permutation below is
    # at least reproducible.
    @assert ε[L] < 0 < ε[L+1]
    copyto!(Tscr, T)
    copyto!(εscr, ε)
    for j in 1:L
        for i in 1:n
            T[i, j]   = Tscr[i, L+j]
            T[i, L+j] = Tscr[i, L+1-j]
        end
        ε[j]   = εscr[L+j]
        ε[L+j] = εscr[L+1-j]
    end

    # Normalize columns of T so that para-unitarity holds, T† Ĩ T = Ĩ.
    for j in 1:n
        c = 1 / sqrt(abs(ε[j]))
        for i in 1:n
            T[i, j] *= c
        end
    end

    # Inverse of λ are eigenvalues of Ĩ H, or equivalently, of √H Ĩ √H. The first L
    # elements are positive, while the next L elements are negative. Their absolute
    # values are excitation energies for the wavevectors q and -q, respectively.
    @. ε = 1 / ε

    return ε
end

# Like `bogoliubov!`, but a Hamiltonian that is not positive definite throws an
# `InstabilityError` with the message `msg()`, built only on failure
function bogoliubov_or_throw!(msg, ws::BogoliubovWorkspace, H::Matrix{<:Complex})
    try
        return bogoliubov!(ws, H)
    catch err
        err isa PosDefException || rethrow()
        rethrow(InstabilityError(msg()))
    end
end

# Variant for callers outside a hot loop, which allocates a workspace per call. Unlike
# the method above, the returned energies are freshly allocated, and so remain valid.
function bogoliubov!(T::Matrix{ComplexF64}, H::Matrix{ComplexF64})
    L = div(size(H, 1), 2)
    @assert size(T) == size(H) == (2L, 2L)
    ws = BogoliubovWorkspace(L)
    bogoliubov!(ws, H)
    copyto!(T, ws.T)
    return ws.ε
end

