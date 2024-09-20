# Reference implementation, without consideration of speed.
function lanczos_ref(A, S, v; niters)
    αs = Float64[]
    βs = Float64[]

    v /= sqrt(real(dot(v, S, v)))
    w = A * S * v
    α = real(dot(w, S, v))
    w = w - α * v
    push!(αs, α)

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
    end

    return SymTridiagonal(αs, βs)
end



# Use Lanczos iterations to find a truncated tridiagonal representation of A S,
# up to similarity transformation. Here, A is any Hermitian matrix, while S is
# both Hermitian and positive definite. Traditional Lanczos uses the identity
# matrix, S = I. The extension to non-identity matrices S is as follows: Each
# matrix-vector product A v becomes A S v, and each vector inner product w† v
# becomes w† S v. The implementation below follows the notation of Wikipedia,
# and is the most stable of the four variants considered by Paige [1]. Each
# Lanczos iteration requires only one matrix-vector multiplication for A and S,
# respectively.

# Similar generalizations of Lanczos have been considered in [2] and [3].
#
# 1. C. C. Paige, IMA J. Appl. Math., 373-381 (1972),
# https://doi.org/10.1093%2Fimamat%2F10.3.373.
# 2. H. A. van der Vorst, Math. Comp. 39, 559-561 (1982),
# https://doi.org/10.1090/s0025-5718-1982-0669648-0
# 3. M. Grüning, A. Marini, X. Gonze, Comput. Mater. Sci. 50, 2148-2156 (2011),
# https://doi.org/10.1016/j.commatsci.2011.02.021.
function lanczos(mulA!, mulS!, v; niters)
    αs = Float64[]
    βs = Float64[]

    vp  = zero(v)
    Sv  = zero(v)
    Svp = zero(v)
    w   = zero(v)
    Sw  = zero(v)

    mulS!(Sv, v)
    c = sqrt(real(dot(v, Sv)))
    @. v /= c
    @. Sv /= c

    mulA!(w, Sv)
    α = real(dot(w, Sv))
    @. w = w - α * v
    mulS!(Sw, w)
    push!(αs, α)

    for _ in 2:niters
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
    end

    return SymTridiagonal(αs, βs)
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
    tridiag = lanczos(mulA!, mulS!, v; niters)
    return eigmin(tridiag), eigmax(tridiag)
end
