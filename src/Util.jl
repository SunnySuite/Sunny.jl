"Mod functions for CartesianIndex"
@inline function modc(i::CartesianIndex{D}, m) :: CartesianIndex{D} where {D}
    CartesianIndex(mod.(Tuple(i), Tuple(m)))
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