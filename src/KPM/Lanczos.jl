# Reference implementation, without consideration of speed.
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
    return Q, T
end


# Apply Lanczos iterations to a Hermitian matrix A, with the option to specify
# an inner-product measure S that must be both Hermitian and positive definite.
# Loosely speaking, Lanczos finds a low-rank approximant A S ≈ Q T Q† S, where T
# is tridiagonal and Q is orthogonal in the measure S, i.e. Q† S Q = I. The
# first column of Q is the initial vector v. The approximant is useful purposes
# of multiplying by v. For example, f(A S) v ≈ Q f(T) e₁ is a near-optimal
# approximation within the Krylov space, for any function f [1]
#
# Returns T, as well as Q† lhs, if left-hand side vectors are specified.
#
# The implementation below follows the notation of Wikipedia, and is the most
# stable of the four variants considered by Paige [2].
#
# The generalization to non-identity matrices S is implemented as follows: Each
# matrix-vector product A v becomes A S v, and each vector inner product w† v
# becomes w† S v. Each Lanczos iteration requires only one matrix-vector
# multiplication for A and S, respectively. Similar generalizations of Lanczos
# have been considered in [3] and [4].
#
# 1. N. Amsel, T. Chen, A. Greenbaum, C. Musco, C. Musco, Near-optimal
#    approximation of matrix functions by the Lanczos method, (2013)
#    https://arxiv.org/abs/2303.03358v2.
# 2. C. C. Paige, IMA J. Appl. Math., 373-381 (1972),
#    https://doi.org/10.1093%2Fimamat%2F10.3.373.
# 3. H. A. van der Vorst, Math. Comp. 39, 559-561 (1982),
#    https://doi.org/10.1090/s0025-5718-1982-0669648-0
# 4. M. Grüning, A. Marini, X. Gonze, Comput. Mater. Sci. 50, 2148-2156 (2011),
#    https://doi.org/10.1016/j.commatsci.2011.02.021.
function lanczos(mulA!, mulS!, v; niters, lhs=zeros(length(v), 0))
    αs = Float64[]
    βs = Float64[]

    vp  = zero(v)
    Sv  = zero(v)
    Svp = zero(v)
    w   = zero(v)
    Sw  = zero(v)
    Q_adj_lhs = zeros(niters, size(lhs, 2))

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
    mul!(view(Q_adj_lhs, 1:1, :), v', lhs)

    for i in 2:niters
        β = sqrt(real(dot(w, Sw)))
        iszero(β) && break
        @. vp = w / β
        @. Svp = Sw / β
        mulA!(w, Svp)
        α = real(dot(w, Svp))
        @. w = w - α * vp - β * v
        mulS!(Sw, w)
        @. v = vp
        push!(αs, α)
        push!(βs, β)
        mul!(view(Q_adj_lhs, i:i, :), v', lhs)
    end

    return SymTridiagonal(αs, βs), Q_adj_lhs
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
    T, _ = lanczos(mulA!, mulS!, v; niters)
    return eigmin(T), eigmax(T)
end
