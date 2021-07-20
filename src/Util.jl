const Vec3 = SVector{3, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}

"Mod functions for CartesianIndex"
@inline function modc(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod.(Tuple(i), Tuple(m)))
end
@inline function modc1(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod1.(Tuple(i), Tuple(m)))
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