# Reference implementation of the generalized Lanczos method, without
# consideration of speed. The input A can be any Hermitian matrix. The optional
# argument S is an inner-product that must be both Hermitian and positive
# definite. Intuitively, Lanczos finds a low-rank approximant A S ≈ Q T Q† S,
# where T is tridiagonal and Q is orthogonal in the sense that Q† S Q = I. The
# first column of Q is the initial vector v, which must be normalized. The
# Lanczos method is useful in two ways. First, eigenvalues of T are
# representative of extremal eigenvalues of A (in an appropriate S-measure
# sense). Second, one obtains a very good approximation to the matrix-vector
# product, f(A S) v ≈ Q f(T) e₁, valid for any function f [1].
#
# The generalization of Lanczos to an arbitrary measure S is implemented as
# follows: Each matrix-vector product A v is replaced by A S v, and each vector
# inner product w† v is replaced by w† S v. Similar generalizations of Lanczos
# have been considered in [2] and [3].
#
# 1. N. Amsel, T. Chen, A. Greenbaum, C. Musco, C. Musco, Near-optimal
#    approximation of matrix functions by the Lanczos method, (2013)
#    https://arxiv.org/abs/2303.03358v2.
# 2. H. A. van der Vorst, Math. Comp. 39, 559-561 (1982),
#    https://doi.org/10.1090/s0025-5718-1982-0669648-0
# 3. M. Grüning, A. Marini, X. Gonze, Comput. Mater. Sci. 50, 2148-2156 (2011),
#    https://doi.org/10.1016/j.commatsci.2011.02.021.
function lanczos_ref(A, S, v; niters)
    αs = Float64[]
    βs = Float64[]
    vs = typeof(v)[]

    c = sqrt(real(dot(v, S, v)))
    @assert c ≈ 1 "Initial vector not normalized"
    v /= c

    w = A * S * v
    α = real(dot(w, S, v))
    w = w - α * v
    push!(αs, α)
    push!(vs, v)

    for _ in 2:niters
        β = sqrt(real(dot(w, S, w)))
        iszero(β) && break
        vp = w / β
        w = A * S * vp
        α = real(dot(w, S, vp))
        w = w - α * vp - β * v
        v = vp
        push!(αs, α)
        push!(βs, β)
        push!(vs, v)
    end

    Q = reduce(hcat, vs)
    T = SymTridiagonal(αs, βs)
    return (; Q, T)
end


# An optimized implementation of the generalized Lanczos method. See
# `lanczos_ref` for explanation, and a reference implementation. No explicit
# orthogonalization is performed. This optimized implementation follows "the
# most stable" of the four variants considered by Paige [1], as listed on
# Wikipedia. Here, however, we allow for an arbitary measure S. In practice,
# this means that each matrix-vector product A v is replaced by A S v, and each
# inner product w† v is replaced by w† S v.
#
# In this implementation, each Lanczos iteration requires only one matrix-vector
# multiplication involving A and S, respectively. 
#
# Details:
#
# * The return value is a pair, (T, lhs† Q). The symbol `lhs` denotes "left-hand
#   side" column vectors, packed horizontally into a matrix. If the `lhs`
#   argument is ommitted, its default value will be a matrix of width 0.
# * The initial vector `v` must be normalized such that `√v S v ≈ 1`. (An
#   alternative implementation might absorb this scale into the return value
#   `lhs† Q`.)
# * After `min_iters` iterations, this procedure will estimate the spectral
#   bandwidth Δϵ between extreme eigenvalues of T. Then `Δϵ / resolution` will
#   be an upper bound on the total number of iterations (which includes the
#   initial `min_iters` iterations as well as subsequent ones).
# * The number of iterations will never exceed length(v). If this limit is hit
#   then, mathematically, the Krylov subspace should be complete, and the matrix
#   T should be exactly similar to the matrix A S. In practice, however,
#   numerical error at finite precision may be severe.
# * Accumulation of numerical roundoff error can lead to significant artifacts
#   in some cases. See https://github.com/SunnySuite/Sunny.jl/pull/339 for
#   examples. An explicit orthogonalization step may be considered for future
#   work. But note that there can be mixing between near-degenerate eigenvalues
#   even when the number of Lanczos iterations is small.
#
# 1. C. C. Paige, IMA J. Appl. Math., 373-381 (1972),
#    https://doi.org/10.1093%2Fimamat%2F10.3.373.

function lanczos(mulA!, mulS!, v; min_iters, resolution=Inf, lhs=zeros(length(v), 0), verbose=false)
    αs = Float64[]
    βs = Float64[]
    lhs_adj_Q = Vector{ComplexF64}[]

    v = copy(v)
    vp  = zero(v)
    Sv  = zero(v)
    Svp = zero(v)
    w   = zero(v)
    Sw  = zero(v)

    mulS!(Sv, v)
    c = sqrt(real(dot(v, Sv)))
    @assert c ≈ 1 "Initial vector not normalized"
    @. v /= c
    @. Sv /= c

    mulA!(w, Sv)
    α = real(dot(w, Sv))
    @. w = w - α * v
    mulS!(Sw, w)
    push!(αs, α)
    push!(lhs_adj_Q, lhs' * v)

    # The maximum number of iterations is length(v)-1, because that is the
    # dimension of the full vector space. Break out early if niters has been
    # reached, which will be set according to the spectral bounds Δϵ after
    # iteration i == min_iters.
    niters = typemax(Int)
    for i in 1:length(v)-1
        if i == min_iters
            T = SymTridiagonal(αs, βs)
            Δϵ = eigmax(T) - eigmin(T)
            niters = max(min_iters, fld(Δϵ, resolution))
            niters += mod(niters, 2) # Round up to nearest even integer
            if verbose
                println("Δϵ=$Δϵ, niters=$niters")
            end
        end

        i >= niters && break

        # Let β = w† S w. If β < 0 then S is not a valid positive definite
        # measure. If β = 0, this would formally indicate that v is a linear
        # combination of only a subset of the eigenvectors of A S. In this case,
        # it is valid to break early for purposes of approximating the
        # matrix-vector product f(A S) v. 
        β² = real(dot(w, Sw))
        iszero(β²) && break
        β² < 0 && error("S is not a positive definite measure")

        β = sqrt(β²)
        @. vp = w / β
        @. Svp = Sw / β
        mulA!(w, Svp)
        α = real(dot(w, Svp))
        @. w = w - α * vp - β * v
        mulS!(Sw, w)
        @. v = vp
        push!(αs, α)
        push!(βs, β)
        push!(lhs_adj_Q, lhs' * v)
    end

    return SymTridiagonal(αs, βs), reduce(hcat, lhs_adj_Q)
end


"""
    eigbounds(swt, niters)

Returns estimates of the extremal eigenvalues of the matrix (Ĩ D)^2 using the
Lanczos algorithm. Here Ĩ is the Bogoliubov identity and D is the dynamical
matrix from LSWT for the wavevector q_reshaped. `niters` should be given a value
smaller than the dimension of `A`.
"""
function eigbounds(swt, q_reshaped, niters)
    q_reshaped = Vec3(q_reshaped)
    L = nbands(swt)

    # w = Ĩ v
    function mulA!(w, v)
        @views w[1:L]    = +v[1:L]
        @views w[L+1:2L] = -v[L+1:2L]
    end

    # w = D v
    function mulS!(w, v)
        mul_dynamical_matrix!(swt, reshape(w, 1, :), reshape(v, 1, :), [q_reshaped])
    end

    v = randn(ComplexF64, 2L)
    Sv = zero(v)
    mulS!(Sv, v)
    v /= sqrt(v' * Sv)
    T, _ = lanczos(mulA!, mulS!, v; min_iters=niters)
    return eigmin(T), eigmax(T)
end
