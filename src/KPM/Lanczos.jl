
"""
    parallel_lanczos!(buf1, buf2, buf3; niters, mat_vec_mul!, inner_product!)

Parallelized Lanczos. Takes functions for matrix-vector multiplication and inner
products. Returns a list of tridiagonal matrices. `niters` sets the number of
iterations used by the Lanczos algorithm.

References:
- [H. A. van der Vorst, Math. Comp. 39, 559-561
  (1982)](https://doi.org/10.1090/s0025-5718-1982-0669648-0).
- [M. Grüning, A. Marini, X. Gonze, Comput. Mater. Sci. 50, 2148-2156
  (2011)](https://doi.org/10.1016/j.commatsci.2011.02.021).
"""
function parallel_lanczos!(vprev, v, w;  niters, mat_vec_mul!, inner_product)
    # ...
end


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
