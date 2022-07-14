"Mod functions for CartesianIndex"
@inline function modc(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod.(Tuple(i), Tuple(m)))
end
@inline function modc1(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod1.(Tuple(i), Tuple(m)))
end
@inline function offset(i::CartesianIndex{D}, n::SVector{D,Int}, m) :: CartesianIndex{D} where {D}
    CartesianIndex( Tuple(mod1.(Tuple(i) .+ n, Tuple(m))) )
end
"Splits a CartesianIndex into its first index, and the rest"
@inline function splitidx(i::CartesianIndex{D}) where {D}
    return (CartesianIndex(Tuple(i)[1:3]), i[4])
end

# Taken from:
# https://discourse.julialang.org/t/efficient-tuple-concatenation/5398/8
"Functions for joining tuples"
@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

# For efficiency, may need to look into Base.unsafe_wrap
#   and pointer trickery if we want to stick with Vec3.

"Reinterprets an array of Vec3 to an equivalent array of Float64"
@inline function _reinterpret_from_spin_array(A::Array{Vec3}) :: Array{Float64}
    Ar = reinterpret(reshape, Float64, A)
end

"Reinterprets an array of Floats with leading dimension 3 to an array of Vec3"
@inline function _reinterpret_to_spin_array(A::Array{Float64}) :: Array{Vec3}
    Ar = reinterpret(reshape, Vec3, A)
end

"Reinterprets an array of Mat3 to an equivalent array of Float64"
@inline function _reinterpret_dipole_tensor(A::OffsetArray{Mat3}) :: Array{Float64}
    Ar = reinterpret(reshape, Float64, parent(A))
    return reshape(Ar, 3, 3, size(A)...)    # make sure this doesn't mess up indexing
end

"Generalize dot product so works on vector of operators"
function LinearAlgebra.dot(a::Vec3, b::NTuple{3, Matrix{ComplexF64}}) where T
    a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
end

"Sparse tensor type for quartic anisotropies. (Could be used for quadratic as well.)"
struct SparseTensor{R,N}
    indices :: NTuple{N, NTuple{R, Int64}}    
    vals    :: SVector{N, Float64}
end

function SparseTensor(tens::T)  where T <: AbstractArray
    # Find non-zero elements
    rank = length(size(tens))
    indices = findall(!iszero, tens)
    vals = tens[indices]
    if length(indices) < 1
        SparseTensor{rank, 0}((), SVector{0, Float64}())
    end

    # Determine number of entries and rank of tensor
    num_entries = length(vals)

    # Convert to static types
    indices = map(i -> i.I, indices)
    indices = NTuple{num_entries, NTuple{rank, Int64}}(indices)
    vals = SVector{num_entries, Float64}(vals)

    SparseTensor{rank, num_entries}(indices, vals)
end


# Note that this approach requires a small allocation, since you can't
# setindex a static array. However, seems to be no slower than old approach.
"Calculate the gradient of the tensor contraction. Used for -∇E calculations in LL dynamics. Only for
rank-4 tensors."
@inline function grad_contract(tens::SparseTensor{4,N}, v::Vec3) where {R,N}
    result = zeros(Float64, 3)
    @inbounds for (indices, val) in zip(tens.indices, tens.vals)
        i1, i2, i3, i4 = indices
        result[i1] += val *         v[i2] * v[i3] * v[i4]
        result[i2] += val * v[i1] *         v[i3] * v[i4]
        result[i3] += val * v[i1] * v[i2] *         v[i4]
        result[i4] += val * v[i1] * v[i2] * v[i3]
    end
    Vec3(result)
end

"Non-generated version of tensor contraction using SparseTensor type. Only for rank-4 tensors."
@inline function contract(tens::SparseTensor{4,N}, v) where {R, N}
    result = zero(v[1])
    @inbounds for (indices, val) in zip(tens.indices, tens.vals) 
        i1, i2, i3, i4 = indices
        result += val*v[i1]*v[i2]*v[i3]*v[i4] 
    end
    result
end