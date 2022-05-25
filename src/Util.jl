"Mod functions for CartesianIndex"
@inline function modc(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod.(Tuple(i), Tuple(m)))
end
@inline function modc1(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod1.(Tuple(i), Tuple(m)))
end
@inline function offset(i::CartesianIndex{D}, n::SVector{D,Int}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(Tuple(mod1.(Tuple(i) .+ n, Tuple(m))))
end
"Splits a CartesianIndex into its first index, and the rest"
@inline function splitidx(i::CartesianIndex{D}) where {D}
    return (i[1], CartesianIndex(Tuple(i)[2:end]))
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
    return reshape(Ar, 3, 3, size(A)...)
end

"Generalize dot product so works on vector of operators"
function LinearAlgebra.dot(a::Vec3, b::NTuple{3, Matrix{ComplexF64}}) where T
    a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
end

"Sparse tensor type for quartic anisotropies. (Could be used for quadratic as well.)"
struct SparseTensor{R,N}
    indices :: NTuple{N, NTuple{R, Int64}}    
    counts  :: NTuple{N, NTuple{3, Int64}}
    vals    :: SVector{N, Float64}
end

numberof(i, v) = filter(x -> x == i, v) |> length
counts_from_indices(indices::NTuple{N, Int64}) where N = map(x -> numberof(x, indices), (1,2,3))
function SparseTensor(tens::T)  where T <: AbstractArray
    # Find non-zero elements
    indices = findall(!iszero, tens)
    vals = tens[indices]
    if length(indices) < 1
        @error "Cannot create a tensor with all zero entries"
    end

    # Determine number of entries and rank of tensor
    rank = length(size(tens))
    num_entries = length(vals)

    # Convert to static types
    indices = map(i -> i.I, indices)
    indices = NTuple{num_entries, NTuple{rank, Int64}}(indices)
    vals = SVector{num_entries, Float64}(vals)

    # Get counts of each index for gradient calculations
    counts = counts_from_indices.(indices)

    SparseTensor{rank, num_entries}(indices, counts, vals)
end


"Conracts a sparse tensor of arbitrary rank R that R vectors (or one vector R times)"
@generated function contract(T::SparseTensor{R,N}, v...) where {R, N}
    repeat_vec = length(v) < R ? true : false
    ex = nothing
    for n ∈ 1:N
        prod_ex = :(T.vals[$n])
        for r ∈ 1:R
            idx = repeat_vec ? 1 : r
            prod_ex = :($prod_ex * v[$idx][T.indices[$n][$r]])
        end
        ex = n == 1 ? prod_ex : :($ex + $prod_ex)
    end
    :(@inbounds $ex)
end

# Can write generated function to do this that would unroll loop and eliminate calls
# to pow. Note that it will not be able to the check for zero values, which is kind
# of an ugly optimization in any case. Will have to benchmark.
"Calculate the gradient of the tensor contraction, for -∇E calculations in LL dynamics"
@inline function grad_contract(T::SparseTensor{R,N}, v) where {R,N}
    result = zero(v)
    for (counts, val) ∈ zip(T.counts, T.vals)
        result += SA[counts[1] == zero(v[1]) ? zero(v[1]) : counts[1] * val * (v[1]^(counts[1]-1)) * (v[2]^counts[2])     * (v[3]^counts[3]),
                     counts[2] == zero(v[1]) ? zero(v[1]) : counts[2] * val * (v[1]^counts[1])     * (v[2]^(counts[2]-1)) * (v[3]^counts[3]),
                     counts[3] == zero(v[1]) ? zero(v[1]) : counts[3] * val * (v[1]^counts[1])     * (v[2]^counts[2])     * (v[3]^(counts[3]-1))]
    end
    result
end

# Note, this seems to be just as fast as generated function version in most cases.
# Both this and generated function substantially faster than full contraction on a Quad3.
"Non-generated version of tensor contraction using SparseTensor type. For comparison."
@inline function contract_sparse(T::SparseTensor{R,N}, v) where {R, N}
    result = zero(v[1])
    @inbounds for (indices, val) ∈ zip(T.indices, T.vals) 
        i1, i2, i3, i4 = indices
        result += val*v[i1]*v[i2]*v[i3]*v[i4] 
    end
    result
end