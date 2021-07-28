import Base.size

"""Defines lattice/basis vectors, and number of unit cells.
   Note this struct does not hold the degrees of freedom,
     but just simply sets the geometry of the underlying
     lattice.
   Because of typing headaches, L = D^2, and Db = D + 1
    are needed in the definition, but can be automatically
    inferred by the constructor.
"""
struct Lattice{D, L, Db} <: AbstractArray{SVector{D, Float64}, Db}
    lat_vecs     :: SMatrix{D, D, Float64, L}       # Columns are lattice vectors
    basis_vecs   :: Vector{SVector{D, Float64}}     # Each SVector gives a basis vector
    size         :: SVector{D, Int}                 # Number of cells along each dimension
end

function Lattice(lat_vecs::SMatrix{D, D, Float64}, basis_vecs::Vector{SVector{D, Float64}}, size::SVector{D, Int}) where {D}
    return Lattice{D, D*D, D+1}(lat_vecs, basis_vecs, size)
end

function Lattice(lat_vecs::Array{Float64, 2}, basis_vecs::Vector{Vector{Float64}}, latsize::Vector{Int})
    D = size(lat_vecs, 1)
    @assert all(map(s->s==D, size(lat_vecs)))
    @assert all(map(v->length(v)==D, basis_vecs))
    @assert length(latsize) == D

    lat_vecs = SMatrix{D, D}(lat_vecs)
    basis_vecs = map(v->SVector{D}(v), basis_vecs)
    latsize = SVector{D, Int}(latsize)

    return Lattice{D, D*D, D+1}(lat_vecs, basis_vecs, latsize)
end

# Calculation taken from the _cell function of the TRI class of:
#  https://wiki.fysik.dtu.dk/ase/_modules/ase/lattice.html#BravaisLattice
"Specify a 3D Bravais Lattice by the lattice vector lengths and angles"
function Lattice(a::Float64, b::Float64, c::Float64, α::Float64, β::Float64, γ::Float64, size::Vector{Int})
    @assert all(map(x->0 < x ≤ π/2, (α, β, γ)))

    sγ, cγ = sin(γ), cos(γ)
    cβ, cα = cos(β), cos(α)
    v1 = [a, 0.0, 0.0]
    v2 = [b * cγ, b * sγ, 0.0]
    v3x = c * cγ
    v3y = c / sγ * (cα - cβ * cγ)
    v3z = c / sγ * √(sγ^2 - cα^2 - cβ^2 + 2 * cα * cβ * cγ)
    v3 = [v3x, v3y, v3z]

    lat_vecs = SMatrix{3, 3}([v1 v2 v3])
    basis_vecs = [SVector{3, Float64}([0.0, 0.0, 0.0])]
    size = SVector{3, Int}(size)
    return Lattice{3, 9, 4}(lat_vecs, basis_vecs, size)
end

@inline function nbasis(lat::Lattice) :: Int
    length(lat.basis_vecs)
end

"Calculate the total volume of the lattice (unit cell volume × num cells)"
@inline function volume(lat::Lattice) :: Float64
    abs(det(lat.lat_vecs)) * prod(lat.size)
end

"Produces an iterator over all (j, k, l) indexes for the lattice"
@inline function bravindexes(lat::Lattice{D}) :: CartesianIndices where {D}
    return CartesianIndices(Tuple(lat.size))
end

@inline function Base.size(lat::Lattice)
    nb = length(lat.basis_vecs)
    return (nb, lat.size...)
end

#=== Indexing returns absolute coordinates of lattice points ===#

@inline function Base.getindex(lat::Lattice{D}, b::Int64, brav::CartesianIndex{D}) where {D}
    return lat.lat_vecs * convert(SVector{D, Float64}, brav) + lat.basis_vecs[b]
end

@inline function Base.getindex(lat::Lattice{D}, b::Int64, brav::NTuple{D, Int64}) where {D}
    return lat.lat_vecs * convert(SVector{D, Float64}, brav) + lat.basis_vecs[b]
end

@inline function Base.getindex(lat::Lattice{D}, b::Int64, brav::Vararg{Int64, D}) where {D}
    return lat.lat_vecs * convert(SVector{D, Float64}, brav) + lat.basis_vecs[b]
end

# TODO: Should just be another Lattice
# TODO: Potentially confusing. These are the reciprocal lattice vectors
#        for the whole box, repeating in space. This is what Ewald wants.
#       However, we want the reciprocal vectors for the underlying (smaller)
#        bravais cell for structure factor calculations.
#       These differ just by a scale factor, but need to decide which
#        convention the struct uses. (Probably the latter.)
"Defines a reciprocal lattice structure"
struct ReciprocalLattice{D, L}
    lat_vecs  :: SMatrix{D, D, Float64, L}     # Columns of this are the reciprocal lattice vectors
    size      :: SVector{D, Int}
end

function Base.eachindex(lat::ReciprocalLattice) :: CartesianIndices
    return CartesianIndices(Tuple(lat.size))
end

"Generates a reciprocal lattice from a real-space Lattice"
function gen_reciprocal(lat::Lattice) :: ReciprocalLattice
    recip_vecs = 2π * transpose(inv(lat.lat_vecs))
    return ReciprocalLattice(recip_vecs, lat.size)
end

"Returns just the underlying Bravais lattice"
function brav_lattice(lat::Lattice{D}) :: Lattice{D} where {D}
    return Lattice(
        lat.lat_vecs,
        [@SVector zeros(D)],
        lat.size
    )
end

""" Some functions which construct common lattices
"""

function cubic_lattice(::Val{D}, a::Float64, latsize) :: Lattice{D} where {D}
    lat_vecs = SA[ a  0.0 0.0;
                  0.0  a  0.0;
                  0.0 0.0  a ]
    basis_vecs = [SA[0.0, 0.0, 0.0]]
    latsize = SVector{D, Int}(latsize)
    Lattice{D, D*D, D+1}(lat_vecs, basis_vecs, latsize)
end

function diamond_lattice(a::Float64, latsize) :: Lattice{3, 9, 4}
    lat_vecs = SA[ 4a 0.0 0.0;
                  0.0  4a 0.0;
                  0.0 0.0  4a]
    basis_vecs = [
        SA[0.0, 0.0, 0.0],
        SA[0.0,  2a,  2a],
        SA[ 2a, 0.0,  2a],
        SA[ 2a,  2a, 0.0],
        SA[ 3a,  3a,  3a],
        SA[ 3a,   a,   a],
        SA[  a,  3a,   a],
        SA[  a    a   3a]
    ]
    latsize = SVector{3, Int}(latsize)
    Lattice{3, 9, 4}(lat_vecs, basis_vecs, latsize)
end

function diamond_lattice2(a::Float64, latsize) :: Lattice{3, 9, 4}
    lat_vecs = SA[a/2 a/2 0.0;
                  a/2 0.0 a/2;
                  0.0 a/2 a/2]
    basis_vecs = [
        SA[0.0, 0.0, 0.0],
        SA[0.25a, 0.25a, 0.25a],
    ]
    latsize = SVector{3, Int}(latsize)
    Lattice{3, 9, 4}(lat_vecs, basis_vecs, latsize)
end

function kagome_lattice(a::Float64, latsize) :: Lattice{2, 4, 3}
    lat_vecs = SA[   a  0.0 ;
                  0.5a √3a/2]
    basis_vecs = [
        SA[0.0,     0.0],
        SA[0.5a,    0.0],
        SA[0.25a, √3a/4]
    ]
    latsize = SVector{2, Int}(latsize)
    Lattice{2, 4, 3}(lat_vecs, basis_vecs, latsize)
end
