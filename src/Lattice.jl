import Base.size

"""
    Lattice{D, L, Db}

A type holding geometry information about a lattice in a simulation box.
The type parameter `D` represents the dimensionality of the `Lattice`,
 while the others must satisfy `L = D^2, Db = D + 1`.

These other type parameters must be in the definition for technical reasons,
 but are inferred without needing them explicitly provided. For example,
 see the `Lattice{D}` constructor.
"""
struct Lattice{D, L, Db} <: AbstractArray{SVector{D, Float64}, Db}
    lat_vecs      :: SMatrix{D, D, Float64, L}       # Columns are lattice vectors
    basis_vecs    :: Vector{SVector{D, Float64}}     # Each SVector gives a basis vector
    species       :: Vector{String}                  # Indices labeling atom types
    size          :: SVector{D, Int}                 # Number of cells along each dimension

    @doc """
        Lattice{D}(lat_vecs, basis_vecs, species, latsize)

    Construct a `D`-dimensional `Lattice`.
    # Arguments
    - `lat_vecs`: A matrix where the lattice vectors form the columns.
                  Must be `convert`-able into `SMatrix{D, D, Float64, D^2}`.
    - `basis_vecs`: A `Vector` of basis positions of sites within the unit cell,
                     given in fractional coordinates.
                    Each element must be `convert`-able into `SVector{D, Float64}`.
    - `species::Vector{String}`: A list of atomic species identifiers for each site.
                                 Equivalent sites should have the same identifier.
    - `latsize`: Specifies the number of unit cells extending along each lattice vector.
                 Must be `convert`-able into `SVector{D, Int}`.
    """
    function Lattice{D}(lat_vecs, basis_vecs, species, latsize) where {D}
        @assert all(isequal(D), size(lat_vecs))          "All dims of lat_vecs should be equal size"
        @assert all(isequal(D), map(length, basis_vecs)) "All basis_vecs should be size $D to match lat_vecs"
        @assert length(basis_vecs) > 0                   "At least one basis atom required"
        @assert length(basis_vecs) == length(species)    "Length of basis_vecs and species should match"
        @assert all(v->all(0 .<= v .< 1), basis_vecs)    "All basis_vecs should be given in fractional coordinates [0, 1)"
        @assert length(latsize) == D                     "latsize should be size $D to match lat_vecs"
        lat_vecs = SMatrix{D, D}(lat_vecs)
        basis_vecs = map(v->lat_vecs * SVector{D, Float64}(v), basis_vecs)
        latsize = SVector{D, Int}(latsize)
        new{D, D*D, D+1}(lat_vecs, basis_vecs, species, latsize)
    end
end

function Lattice{D}(lat_vecs, basis_vecs, latsize) where {D}
    return Lattice{D}(
        lat_vecs, basis_vecs, fill("A", length(basis_vecs)), latsize
    )
end

"""
    Lattice(lat_vecs::SMatrix{D, D}, basis_vecs, latsize)

Construct a `Lattice`, with the dimension inferred by the shape of `lat_vecs`.
"""
function Lattice(lat_vecs::SMatrix{D,D}, basis_vecs, latsize) where {D}
    return Lattice{D}(
        lat_vecs, basis_vecs, fill("A", length(basis_vecs)), latsize
    )
end

"""
    lattice_params(lattice::Lattice{3})
"""

function Base.display(lattice::Lattice)
    D = length(size(lattice)) - 1
    println(join(size(lattice), "x"), " Lattice{$D}")
    # Print out lattice vectors, species, basis?
end

# Constructors converting Lattice -> Crystal, and Crystal -> Lattice

"""
    Crystal(lattice::Lattice)

Construct a `Crystal` using geometry information in `lattice`, inferring symmetry information.
"""
function Crystal(lattice::Lattice{3, 9, 4})
    L = lattice.lat_vecs
    # Convert absolute basis positions to fractional coordinates
    basis_coords = map(b -> inv(L) * b, lattice.basis_vecs)
    Crystal(L, basis_coords, lattice.species)
end

function Crystal(lattice::Lattice{2, 4, 3})
    L = lattice.lat_vecs

    # Expand the lattice to 3D, with a very long c axis
    L3 = @MMatrix zeros(3, 3)
    L3[1:2, 1:2] = L
    # Make the vertical axis 10x longer than the longest 2D lattice vector
    max_len = maximum(norm.(eachcol(lattice.lat_vecs)))
    L3[3, 3] = 10. * max_len
    L3 = SMatrix(L3)

    basis_coords = map(b -> SVector{3}((inv(L) * b)..., 0.), lattice.basis_vecs)

    Crystal(L3, basis_coords, lattice.species)
end

function Lattice(cryst::Crystal, latsize) :: Lattice{3, 9, 4}
    Lattice{3}(cryst.lat_vecs, cryst.positions, cryst.species, latsize)
end

"Return number of basis sites in Lattice"
nbasis(lat::Lattice) = length(lat.basis_vecs)
"Compute the volume of a unit cell."
cell_volume(lat::Lattice) = abs(det(lat.lat_vecs))
"Compute the volume of the full simulation box."
volume(lat::Lattice) = cell_volume(lat) * prod(lat.size)
lattice_vectors(lat::Lattice) = lat.lat_vecs
lattice_params(lattice::Lattice{3}) = lattice_params(lattice.lat_vecs)

"Produce an iterator over all unit cell indices."
@inline function eachcellindex(lat::Lattice{D}) :: CartesianIndices where {D}
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
# TODO: Potentially confusing. `gen_reciprocal` computes the reciprocal lattice
#        vectors of the unit cell, not the entire simulation box. Most use cases
#        want the former, but Ewald summations want the latter.
#       These differ just by a scale factor, which is currently explicitly handled
#        in Ewald summation calculations.
"Defines a reciprocal lattice structure"
struct ReciprocalLattice{D, L}
    lat_vecs  :: SMatrix{D, D, Float64, L}     # Columns of this are the reciprocal lattice vectors
    size      :: SVector{D, Int}
end

function Base.eachindex(lat::ReciprocalLattice) :: CartesianIndices
    return CartesianIndices(Tuple(lat.size))
end

"Generate a `ReciprocalLattice` for a given `Lattice`"
function gen_reciprocal(lat::Lattice) :: ReciprocalLattice
    recip_vecs = 2π * transpose(inv(lat.lat_vecs))
    return ReciprocalLattice(recip_vecs, lat.size)
end

"Returns just the underlying Bravais lattice"
function brav_lattice(lat::Lattice{D}) :: Lattice{D} where {D}
    return Lattice{D}(
        lat.lat_vecs,
        [@SVector zeros(D)],
        ["A"],
        lat.size
    )
end

cell_type(lattice::Lattice{3}) = cell_type(lattice.lat_vecs)

function distance(lat::Lattice{D}, b::Bond{D}) where {D}
    return norm(lat.lat_vecs * b.n + (lat.basis_vecs[b.j] - lat.basis_vecs[b.i]))
end

""" Some functions which construct common lattices
"""
function square_lattice(a::Float64, latsize) :: Lattice{2, 4, 3}
    lat_vecs = SA[ a  0.0;
                  0.0  a ]
    basis_vecs = [SA[0.0, 0.0]]
    basis_labels = ["A"]
    latsize = SVector{D, Int}(latsize)
    Lattice{2}(lat_vecs, basis_vecs, basis_labels, latsize)
end

function cubic_lattice(a::Float64, latsize) :: Lattice{3, 9, 4}
    lat_vecs = SA[ a  0.0 0.0;
                  0.0  a  0.0;
                  0.0 0.0  a ]
    basis_vecs = [SA[0.0, 0.0, 0.0]]
    basis_labels = ["A"]
    latsize = SVector{D, Int}(latsize)
    Lattice{3}(lat_vecs, basis_vecs, basis_labels, latsize)
end

function fcc_conventional(a::Float64, latsize) :: Lattice{3, 9, 4}
    lat_vecs = SA[ a  0.0 0.0;
                  0.0  a  0.0;
                  0.0 0.0  a ]
    basis_vecs = [SA[ 0.,  0.,  0.],
                  SA[1/2, 1/2, 0.0],
                  SA[1/2,   0, 1/2],
                  SA[  0, 1/2, 1/2]]
    basis_labels = fill("A", 4)
    latsize = SVector{3, Int}(latsize)
    Lattice{3}(lat_vecs, basis_vecs, basis_labels, latsize)
end

function diamond_conventional(a::Float64, latsize) :: Lattice{3, 9, 4}
    lat_vecs = SA[ 4a 0.0 0.0;
                  0.0  4a 0.0;
                  0.0 0.0  4a]
    basis_vecs = [
        SA[0.0, 0.0, 0.0],
        SA[0.0, 1/2, 1/2],
        SA[1/2, 0.0, 1/2],
        SA[1/2, 1/2, 0.0],
        SA[3/4, 3/4, 3/4],
        SA[3/4, 1/4, 1/4],
        SA[1/4, 3/4, 1/4],
        SA[1/4, 1/4, 3/4]
    ]
    basis_labels = fill("A", 8)
    latsize = SVector{3, Int}(latsize)
    Lattice{3}(lat_vecs, basis_vecs, basis_labels, latsize)
end

function diamond_primitive(a::Float64, latsize) :: Lattice{3, 9, 4}
    lat_vecs = SA[a/2 a/2 0.0;
                  a/2 0.0 a/2;
                  0.0 a/2 a/2]
    basis_vecs = [
        SA[0.0, 0.0, 0.0],
        SA[1/4, 1/4, 1/4],
    ]
    basis_labels = ["A", "A"]
    latsize = SVector{3, Int}(latsize)
    Lattice{3}(lat_vecs, basis_vecs, basis_labels, latsize)
end

function kagome_lattice(a::Float64, latsize) :: Lattice{2, 4, 3}
    lat_vecs = SA[  a    0.5a;
                   0.0  √3a/2]
    basis_vecs = [
        SA[0.0, 0.0],
        SA[1/2, 0.0],
        SA[0.0, 1/2]
    ]
    basis_labels = ["A", "A", "A"]
    latsize = SVector{2, Int}(latsize)
    Lattice{2}(lat_vecs, basis_vecs, basis_labels, latsize)
end
