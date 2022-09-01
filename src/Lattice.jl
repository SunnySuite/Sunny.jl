"""Definition and functions working with the Lattice type. Lattice is a lightweight
    wrapper that stores lattice geometry in absolute coordinates for fast site indexing,
    along with a box size for easy collection of all sites.
   Lattice is arbitrary-dimensional, but contains no symmetry information.
"""

"""
    Lattice

A type holding geometry information about a 3D lattice in a simulation box.
"""
struct Lattice <: AbstractArray{Vec3, 4}
    lat_vecs      :: Mat3                            # Columns are lattice vectors
    basis_vecs    :: Vector{Vec3}                    # Each SVector gives a basis vector
    types         :: Vector{String}                  # Indices labeling atom types
    size          :: SVector{3, Int}                 # Number of cells along each dimension

    @doc """
        Lattice(lat_vecs, basis_vecs, types, latsize)

    Construct a 3-dimensional `Lattice`.
    # Arguments
    - `lat_vecs`: A matrix where the lattice vectors form the columns.
                  Must be `convert`-able into `SMatrix{3, 3, Float64, 9}`.
    - `basis_vecs`: A `Vector` of basis positions of sites within the unit cell,
                     given in fractional coordinates.
                    Each element must be `convert`-able into `SVector{3, Float64}`.
    - `types::Vector{String}`: A list of atomic types identifiers for each site.
                                 Equivalent sites should have the same identifier.
    - `latsize`: Specifies the number of unit cells extending along each lattice vector.
                 Must be `convert`-able into `SVector{3, Int}`.
    """
    function Lattice(lat_vecs, basis_vecs, types, latsize)
        @assert all(isequal(3), size(lat_vecs))          "All `lat_vecs` dims must be length 3"
        @assert all(isequal(3), map(length, basis_vecs)) "All basis_vecs should be length 3 to match lat_vecs"
        @assert length(basis_vecs) > 0                   "At least one basis atom required"
        @assert length(basis_vecs) == length(types)      "Length of basis_vecs and types should match"
        @assert all(v->all(0 .<= v .< 1), basis_vecs)    "All basis_vecs should be given in fractional coordinates [0, 1)"
        @assert length(latsize) == 3                     "latsize should be length 3"
        lat_vecs = convert(Mat3, lat_vecs)
        basis_vecs = [lat_vecs * convert(Vec3, b) for b in basis_vecs]
        latsize = SVector{3, Int}(latsize)
        new(lat_vecs, basis_vecs, types, latsize)
    end
end

"""
    Lattice(lat_vecs, basis_vecs, latsize)

Construct a `Lattice`, with the dimension inferred by the shape of `lat_vecs`,
and all sites assumed to be the same type (labeled `"A"`).
"""
Lattice(lvecs, bvecs, latsize) = Lattice(lvecs, bvecs, fill("A", length(bvecs)), latsize)

"""
    lattice_params(lattice::Lattice)
"""

function Base.show(io::IO, ::MIME"text/plain", lattice::Lattice)
    println(io, join(size(lattice), "x"), " Lattice")
    # Print out lattice vectors, types, basis?
end

# Constructors converting Lattice -> Crystal, and Crystal -> Lattice

"""
    Crystal(lattice::Lattice)

Construct a `Crystal` using geometry information in `lattice`, inferring symmetry information.
"""
function Crystal(lattice::Lattice)
    L = lattice.lat_vecs
    # Convert absolute basis positions to fractional coordinates
    basis_coords = map(b -> inv(L) * b, lattice.basis_vecs)
    Crystal(L, basis_coords, types=lattice.types)
end

function Lattice(cryst::Crystal, latsize) :: Lattice
    Lattice(cryst.lat_vecs, cryst.positions, cryst.types, latsize)
end

# Number of sublattices.
nbasis(lat::Lattice) = length(lat.basis_vecs)
# Volume of a unit cell.
cell_volume(lat::Lattice) = abs(det(lat.lat_vecs))
# Volume of the full simulation box.
volume(lat::Lattice) = cell_volume(lat) * prod(lat.size)
lattice_vectors(lat::Lattice) = lat.lat_vecs
lattice_params(lat::Lattice) = lattice_params(lat.lat_vecs)

"Produce an iterator over all unit cell indices."
@inline function eachcellindex(lat::Lattice) :: CartesianIndices{3}
    return CartesianIndices(Tuple(lat.size))
end

# Perhaps it's confusing to have the Lattice.size field refer to just the Bravais
# lattice dimensions (3 numbers), but to have size(lat::Lattice) return dimensions 
# of the whole thing including sites within each unit cell (4 numbers).
@inline function Base.size(lat::Lattice)
    return (lat.size..., nbasis(lat))
end

#=== Indexing returns absolute coordinates of lattice points ===#

@inline function Base.getindex(lat::Lattice, brav::CartesianIndex{3}, b::Int64)
    return muladd(lat.lat_vecs, convert(Vec3, brav), lat.basis_vecs[b])
end

@inline function Base.getindex(lat::Lattice, brav::NTuple{3, Int64}, b::Int64)
    return muladd(lat.lat_vecs, convert(Vec3, brav), lat.basis_vecs[b])
end

@inline function Base.getindex(lat::Lattice, idx::Vararg{Int64, 4})
    brav, b = idx[1:3], idx[end]
    return muladd(lat.lat_vecs, convert(Vec3, brav), lat.basis_vecs[b])
end

# TODO: Should just be another Lattice
# TODO: Potentially confusing. `gen_reciprocal` computes the reciprocal lattice
#        vectors of the unit cell, not the entire simulation box. Most use cases
#        want the former, but Ewald summations want the latter.
#       These differ just by a scale factor, which is currently explicitly handled
#        in Ewald summation calculations.
"Defines a reciprocal lattice structure"
struct ReciprocalLattice
    lat_vecs  :: SMatrix{3, 3, Float64, 9}     # Columns of this are the reciprocal lattice vectors
    size      :: SVector{3, Int}
end

function ReciprocalLattice(lat::Lattice)
    recip_vecs = 2π * transpose(inv(lat.lat_vecs))
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
function brav_lattice(lat::Lattice) :: Lattice
    return Lattice(
        lat.lat_vecs,
        [zero(Vec3)],
        ["A"],
        lat.size
    )
end

cell_type(lattice::Lattice) = cell_type(lattice.lat_vecs)

function distance(lat::Lattice, b::Bond)
    return norm(lat.lat_vecs * b.n + (lat.basis_vecs[b.j] - lat.basis_vecs[b.i]))
end
