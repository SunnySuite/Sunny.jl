@inline δ(x, y) = (x==y)


# Set submatrix H21 to H12' without allocating
function set_H21!(H)
    L = round(Int, size(H, 1)/2)
    H12 = view(H, 1:L,L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    for i in CartesianIndices(H21)
        i, j = i.I
        H21[i,j] = conj(H12[j,i])
    end
end

# Calculating norm(H - H') without allocating
function hermiticity_norm(H)
    acc = 0.0
    for idx in CartesianIndices(H) 
        acc += abs2(H[idx] - H'[idx])
    end
    return sqrt(acc)
end

# Modified from LinearAlgebra.jl to not perform any conjugation
function dot_no_conj(x, A, y)
    (axes(x)..., axes(y)...) == axes(A) || throw(DimensionMismatch())
    T = typeof(dot(first(x), first(A), first(y)))
    s = zero(T)
    i₁ = first(eachindex(x))
    x₁ = first(x)
    @inbounds for j in eachindex(y)
        yj = y[j]
        if !iszero(yj)
            temp = zero(A[i₁,j] * x₁)
            @simd for i in eachindex(x)
                temp += A[i,j] * x[i]
            end
            s += temp * yj
        end
    end
    return s
end

# Scalar "metric"
function dot_no_conj(x, A::Float64, y)
    s = zero(eltype(x))
    axes(x) == axes(y) || throw(DimensionMismatch())
    for j in eachindex(x)
        s += x[j]*y[j]
    end
    return A*s
end

# Diagonal "metric"
function dot_no_conj(x, A::SVector{N, Float64}, y) where N
    s = zero(eltype(x))
    axes(x) == axes(y)  == axes(A) || throw(DimensionMismatch())
    for j in eachindex(x)
        s += A[j]*x[j]*y[j]
    end
    return s
end