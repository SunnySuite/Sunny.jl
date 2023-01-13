"Element-wise application of mod1(cell+off, latsize), returning CartesianIndex"
@inline offsetc(cell::CartesianIndex{3}, off, latsize) = CartesianIndex(mod1.(Tuple(cell).+Tuple(off), latsize))

"Split a Cartesian index (i,j,k,b) into its parts (i,j,k) and b."
@inline splitidx(idx::CartesianIndex{4}) = (CartesianIndex((idx[1],idx[2],idx[3])), idx[4])

@inline convert_idx(idx::CartesianIndex{4}) = idx
@inline convert_idx(idx::NTuple{4,Int}) = CartesianIndex(idx)

"Tensor product of 3-vectors"
(⊗)(a::Vec3,b::Vec3) = reshape(kron(a,b), 3, 3)

