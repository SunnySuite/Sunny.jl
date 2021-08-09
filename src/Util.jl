const Vec3 = SVector{3, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}

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

# Taken from:
# https://discourse.julialang.org/t/efficient-tuple-concatenation/5398/8
"Functions for joining tuples"
@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

""" Hacky workaround of the puzzling default behavior of zeros, which is to share memory
     for each element initialized using zeros(T, dims...) for types which allocate their own mem.
"""
function zeros_sepmem(T, dims...)
    arr = Array{T, length(dims)}(undef, dims)
    for i in eachindex(arr)
        arr[i] = zero(T)
    end
    return arr
end

"Constructs 3D rotation matrices, rotating by θ around vector n (right-handed)"
function rot_matrix(n::Vector{Float64}, θ::Float64)
    @assert length(n) == 3

    nx, ny, nz = n / norm(n)
    sθ, cθ = sin(θ), cos(θ)
    R = [[(cθ + nx * nx * (1 - cθ)) (nx * ny * (1 - cθ) - nz * sθ) (nx * nz * (1 - cθ) + ny * sθ)]
         [(ny * nx * (1 - cθ) + nz * sθ) (cθ + ny * ny * (1 - cθ)) (ny * nz * (1 - cθ) - nx * sθ)]
         [(nz * nx * (1 - cθ) - ny * sθ) (nz * ny * (1 - cθ) + nx * sθ) (cθ + nz * nz * (1 - cθ))]]
    return R
end


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