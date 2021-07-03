using TOML
using StaticArrays
using LinearAlgebra
import Base.size

"""Defines lattice/basis vectors, and number of unit cells.
   Note this struct does not hold the degrees of freedom,
     but just simply sets the geometry of the underlying
     lattice.
"""
struct Lattice{D, L}
    lat_vecs     :: SMatrix{D, D, Float64, L}       # Columns are lattice vectors
    basis_vecs   :: Vector{SVector{D, Float64}}     # Each SVector gives a basis vector
    size         :: SVector{D, Int}                 # Number of cells along each dimension
end

"Calculate the total volume of the lattice (unit cell volume × num cells)"
@inline function volume(lat::Lattice) :: Float64
    abs(det(lat.lat_vecs)) * prod(lat.size)
end

"Produces an iterator over all (j, k, l) indexes for the lattice"
@inline function brav_indices(lat::Lattice{D}) :: CartesianIndices where {D}
    return CartesianIndices(Tuple(lat.size))
end

"Produces an iterator over all (j, k, l, basis) indexes for the lattice"
@inline function indices(lat::Lattice{D}) :: CartesianIndices where {D}
    nb = length(lat.basis_vecs)
    return CartesianIndices((lat.size..., nb))
end

@inline function Base.size(lat::Lattice)
    nb = length(lat.basis_vecs)
    return (lat.size..., nb)
end

#=== Indexing returns absolute coordinates of lattice points ===#

@inline function Base.getindex(lat::Lattice{D}, brav::CartesianIndex{D}, b::Int64) where {D}
    return lat.lat_vecs * convert(SVector{D, Float64}, brav) + lat.basis_vecs[b]
end

@inline function Base.getindex(lat::Lattice{D}, brav::NTuple{D, Int64}, b::Int64) where {D}
    return lat.lat_vecs * convert(SVector{D, Float64}, brav) + lat.basis_vecs[b]
end

struct ReciprocalLattice{D, L}
    lat_vecs     :: SMatrix{D, D, Float64, L}     # Columns of this are the reciprocal lattice vectors
    size         :: SVector{D, Int}
end

function indices(lat::ReciprocalLattice) :: CartesianIndices
    return CartesianIndices(Tuple(lat.size))
end

"Generates a reciprocal lattice from a real-space Lattice"
function gen_reciprocal(lat::Lattice) :: ReciprocalLattice
    recip_vecs = 2π * transpose(inv(lat.lat_vecs)) ./ lat.size
    return ReciprocalLattice(recip_vecs, lat.size)
end

struct ValidateError <: Exception end
Base.showerror(io::IO, e::ValidateError) = print(io, e)

"Checks if a config for a lattice is consistent / properly formatted"
function _parse_lattice(config::Dict{String, Any}) :: Lattice
    try
        dim = config["dimension"]
        lat_vecs = config["lattice_vectors"]
        basis_vecs = config["basis_vectors"]
        lattice_size = config["lattice_size"]

        @assert length(lat_vecs) == dim
        @assert all(v -> length(v) == dim, lat_vecs)
        @assert all(v -> length(v) == dim, basis_vecs)
        @assert length(lattice_size) == dim

        lat_vecs = SMatrix{dim, dim, Float64, dim*dim}(hcat(lat_vecs...))
        basis_vecs = [SVector{dim}(bv) for bv in basis_vecs]
        # basis_vecs = hcat(basis_vecs...
        lattice_size = SVector{dim}(lattice_size)

        return Lattice{dim, dim*dim}(lat_vecs, basis_vecs, lattice_size) 
    catch err
        if isa(err, KeyError)
            throw(ValidateError("lattice config missing mandatory key: $(err.key)"))
        else
            # Re-raise all other kinds of errors
            throw(err)
        end
    end
end
