#=

function lanczos(swt, q_reshaped, niters)
    L = nbands(swt)

    function inner_product(w, v)
        return @views real(dot(w[1:L], v[1:L]) - dot(w[L+1:2L], v[L+1:2L]))
    end

    function mat_vec_mul!(w, v, q_reshaped)
        mul_dynamical_matrix!(swt, reshape(w, 1, :), reshape(v, 1, :), [q_reshaped])
        view(w, L+1:2L) .*= -1
    end

    αs = zeros(Float64, niters)   # Diagonal part
    βs = zeros(Float64, niters-1) # Off-diagonal part

    vprev = randn(ComplexF64, 2L)
    v = zeros(ComplexF64, 2L)
    w = zeros(ComplexF64, 2L)

    mat_vec_mul!(w, vprev, q_reshaped)
    b0 = sqrt(inner_product(vprev, w))
    vprev = vprev / b0
    w = w / b0
    αs[1] = inner_product(w, w)
    @. w = w - αs[1]*vprev
    for j in 2:niters
        @. v = w
        mat_vec_mul!(w, v, q_reshaped)
        βs[j-1] = sqrt(inner_product(v, w)) # maybe v?
        if abs(βs[j-1]) < 1e-12
            # TODO: Terminate early if all β are small
            βs[j-1] = 0
            @. v = w = 0
        else
            @. v = v / βs[j-1]
            @. w = w / βs[j-1]
        end
        αs[j] = inner_product(w, w)
        @. w = w - αs[j]*v - βs[j-1]*vprev
        v, vprev = vprev, v
    end

    return SymTridiagonal(αs, βs)
end


"""
    eigbounds(A, niters; extend=0.0)

Returns estimates of the extremal eigenvalues of Hermitian matrix `A` using the
Lanczos algorithm. `niters` should be given a value smaller than the dimension
of `A`. `extend` specifies how much to shift the upper and lower bounds as a
percentage of the total scale of the estimated eigenvalues.
"""
function eigbounds(swt, q_reshaped, niters; extend)
    tridiag = lanczos(swt, Vec3(q_reshaped), niters)
    lo, hi = extrema(eigvals!(tridiag)) # TODO: pass irange=1 and irange=2L
    slack = extend*(hi-lo)
    return lo-slack, hi+slack
end

=#


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
function lanczos3(mulA!, mulS!, v; niters)
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
    eigbounds(swt, niters; extend=0.0)

Returns estimates of the extremal eigenvalues of the matrix (Ĩ D)^2 using the
Lanczos algorithm. Here Ĩ is the Bogoliubov identity and D is the dynamical
matrix from LSWT for the wavevector q_reshaped. `niters` should be given a value
smaller than the dimension of `A`. `extend` specifies how much to shift the
upper and lower bounds as a percentage of the total scale of the estimated
eigenvalues.
"""
function eigbounds(swt, q_reshaped, niters; extend)
    q_reshaped = Vec3(q_reshaped)
    L = nbands(swt)

    function mulA!(w, v)
        @views w[1:L]    = +v[1:L]
        @views w[L+1:2L] = -v[L+1:2L]
    end

    function mulS!(w, v)
        mul_dynamical_matrix!(swt, reshape(w, 1, :), reshape(v, 1, :), [q_reshaped])
    end

    v = randn(ComplexF64, 2L)
    tridiag = lanczos3(mulA!, mulS!, v; niters)
    lo, hi = extrema(eigvals!(tridiag)) # TODO: pass irange=1 and irange=2L
    slack = extend*(hi-lo)
    return lo-slack, hi+slack
end
