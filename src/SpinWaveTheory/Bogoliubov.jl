import LinearAlgebra: BlasInt

# Storage for `bogoliubov!`, holding the transformation `T` and energies `ε`
# that it produces, together with the LAPACK scratch space that `eigen!` would
# otherwise allocate on every call.
struct BogoliubovWorkspace
    T     :: Matrix{ComplexF64}   # Transformation, as documented by `excitations!`
    ε     :: Vector{Float64}      # Signed energies, in the same grouping
    Tscr  :: Matrix{ComplexF64}   # Scratch for the reordering below
    εscr  :: Vector{Float64}
    work  :: Vector{ComplexF64}   # LAPACK scratch for zhegvd
    rwork :: Vector{Float64}
    iwork :: Vector{BlasInt}
    info  :: Base.RefValue{BlasInt}
end

function BogoliubovWorkspace(L::Int)
    n = 2L
    # Minimum scratch sizes documented for ZHEGVD with jobz = 'V' and n > 1
    # (always true here, since n = 2L).
    return BogoliubovWorkspace(zeros(ComplexF64, n, n), zeros(n), zeros(ComplexF64, n, n), zeros(n),
                               zeros(ComplexF64, 2n + n^2), zeros(1 + 5n + 2n^2), zeros(BlasInt, 3 + 5n),
                               Ref(BlasInt(0)))
end

# Bogoliubov transformation that diagonalizes a quadratic bosonic Hamiltonian,
# allowing for anomalous terms. The general procedure derives from Colpa,
# Physica A, 93A, 327-353 (1978). Overwrites data in H. The returned energies
# alias `ws.ε`, and so are valid only until the next call on the same workspace.
function bogoliubov!(ws::BogoliubovWorkspace, H::Matrix{ComplexF64})
    (; T, ε, Tscr, εscr, work, rwork, iwork, info) = ws
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

    # Solve generalized eigenvalue problem, Ĩ t = λ H t, for columns t of T. This is
    # `eigen!(Hermitian(T), Hermitian(H))` with every allocation hoisted into `ws`.
    ccall((BLAS.@blasfunc(zhegvd_), LinearAlgebra.libblastrampoline), Cvoid,
          (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
           Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
           Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
           Clong, Clong),
          1, 'V', 'U', n, T, n, H, n, ε, work, length(work), rwork,
          length(rwork), iwork, length(iwork), info, 1, 1)
    info[] < 0 && throw(ArgumentError("Invalid argument #$(-info[]) to LAPACK zhegvd"))
    info[] > 0 && throw(PosDefException(info[]))

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
function bogoliubov_or_throw!(msg, ws::BogoliubovWorkspace, H::Matrix{ComplexF64})
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

