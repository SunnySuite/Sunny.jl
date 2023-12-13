"""
    eigbounds(A, niters; extend=0.0)

Returns estimates of the extremal eigenvalues of Hermitian matrix `A` using the
Lanczos algorithm. `niters` should be given a value smaller than the dimension
of `A`. `extend` specifies how much to shift the upper and lower bounds as a
percentage of the total scale of the estimated eigenvalues.
"""
function eigbounds(A, niters; extend=0.5)
    lo, hi = lanczos(A, niters) |> eigvals |> xs -> (first(xs), last(xs))
    slack = extend*(hi-lo)
    return lo-slack, hi+slack
end


"""
    lanczos(A, niters)

Takes an `N`x`N`` quasi-Hermitian matrix `A` and returns a real, symmetric,
tridiagonal matrix `T`. `niters` sets the number of iterations used by the
Lanczos algorithm.
"""

function lanczos(A, niters)
    N = size(A, 1)
    nmodes= Int(N/2)
    Ĩ = diagm([ones(nmodes); -ones(nmodes)])
    #@assert Ĩ * A == A' * Ĩ  "Matrix must be quasi-Hermitian"
    αs = zeros(Float64, niters)    # Main diagonal 
    βs = zeros(Float64, niters-1)  # Off diagonal
    buf = zeros(ComplexF64, N, 3)  # Vector buffer -- don't technically need 3 columns, but more legible this way
    lanczos_QH_aux!(αs, βs, buf, A, niters)  # Call non-allocating internal Lanczos function
    return SymTridiagonal(αs, βs)
end

function lanczos_aux!(αs, βs, buf, A, niters)
    v, vprev, w = view(buf,:,1), view(buf,:,2), view(buf,:,3)
    randn!(vprev)
    normalize!(vprev)
    mul!(w, A, vprev)
    αs[1] = real(w⋅vprev)
    @. w = w - αs[1]*vprev

    for j in 2:niters
        βs[j-1] = norm(w)
        iszero(βs[j-1]) && return
        @. v = w/βs[j-1]
        mul!(w, A, v)
        αs[j] = real(w⋅v)
        @. w = w - αs[j]*v - βs[j-1]*vprev
        v, vprev = vprev, v
    end

    return nothing
end
function Ĩ_dot_3(w, v, n)
    @views real(dot(w[1:n], v[1:n]) - dot(w[n+1:2n], v[n+1:2n]))
end
function lanczos_QH_aux!(αs, βs, buf, A, niters)
    v, vprev, w = view(buf,:,1), view(buf,:,2), view(buf,:,3)
    randn!(vprev)
    mat_size=Int(size(A,2)/2)
    mul!(w, A, vprev)
    b0 = sqrt(Ĩ_dot_3(vprev, w,mat_size))
    vprev = vprev/b0
    w = w/b0
    αs[1] = Ĩ_dot_3(w, w,mat_size)
    @. w = w - αs[1]*vprev
    for j in 2:niters
        @. v = w
        mul!(w, A, v)
        βs[j-1] = sqrt(Ĩ_dot_3(v, w,mat_size)) # maybe v?
        iszero(βs[j-1]) && return
        @. v = v/βs[j-1]
        @. w = w/βs[j-1]
        αs[j] = Ĩ_dot_3(w, w,mat_size)
        @. w = w - αs[j]*v - βs[j-1]*vprev
        v, vprev = vprev, v
    end
    return nothing
end




function lanczos_legacy(A, niters)
    @assert ishermitian(A) "Matrix must be Hermitian"
    N = size(A, 1)

    αs = zeros(Float64, niters)    # Main diagonal 
    βs = zeros(Float64, niters-1)  # Off diagonal
    buf = zeros(ComplexF64, N, 3)  # Vector buffer -- don't technically need 3 columns, but more legible this way 

    lanczos_aux!(αs, βs, buf, A, niters)  # Call non-allocating internal Lanczos function
    return SymTridiagonal(αs, βs)
end
