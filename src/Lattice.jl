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
    lattice_params(lat_vecs::Mat3)

Compute the lattice parameters ``(a, b, c, α, β, γ)`` from a set of lattice vectors,
 which form the columns of `lat_vecs`.
"""
function lattice_params(lat_vecs::Mat3) :: NTuple{6, Float64}
    v1, v2, v3 = eachcol(lat_vecs)
    a, b, c = norm(v1), norm(v2), norm(v3)
    α = acosd((v2 ⋅ v3) / (b * c))
    β = acosd((v1 ⋅ v3) / (a * c))
    γ = acosd((v1 ⋅ v2) / (a * b))
    return (a, b, c, α, β, γ)
end

"""
    lattice_params(lattice::Lattice{3})
"""
lattice_params(lattice::Lattice{3}) = lattice_params(lattice.lat_vecs)

function Base.display(lattice::Lattice)
    D = length(size(lattice)) - 1
    println(join(size(lattice), "x"), " Lattice{$D}")
    # Print out lattice vectors, species, basis?
end

"""
    lattice_vectors(a, b, c, α, β, γ) :: Mat3

Compute a set of lattice vectors (forming the columns of the result), specified by a given
 set of lattice parameters ``(a, b, c, α, β, γ)``.
"""
function lattice_vectors(a, b, c, α, β, γ) :: Mat3
    @assert all(0 < x < 180 for x in (α, β, γ))

    sγ, cγ = sind(γ), cosd(γ)
    cβ, cα = cosd(β), cosd(α)
    v1 = Vec3(a, 0, 0)
    v2 = Vec3(b * cγ, b * sγ, 0)
    v3x = c * cβ
    v3y = c / sγ * (cα - cβ * cγ)
    v3z = c / sγ * √(sγ^2 - cα^2 - cβ^2 + 2 * cα * cβ * cγ)
    v3 = Vec3(v3x, v3y, v3z)

    @assert norm(v1) ≈ a
    @assert norm(v2) ≈ b
    @assert norm(v3) ≈ c
    @assert acosd(v1⋅v2 / (a*b)) ≈ γ
    @assert acosd(v1⋅v3 / (a*c)) ≈ β
    @assert acosd(v2⋅v3 / (b*c)) ≈ α

    return [v1 v2 v3]
end

lattice_vectors(lattice::Lattice) = lattice.lat_vecs

"""
    Lattice(a, b, c, α, β, γ)

Specify a 3D Bravais `Lattice` by the lattice vector lengths and angles (in degrees)
"""
function Lattice(a::Float64, b::Float64, c::Float64, α::Float64, β::Float64, γ::Float64, size::Vector{Int})
    lat_vecs = lattice_vectors(a, b, c, α, β, γ)
    basis_vecs = [SVector{3, Float64}([0.0, 0.0, 0.0])]
    size = SVector{3, Int}(size)
    return Lattice{3, 9, 4}(lat_vecs, basis_vecs, ["A"], size)
end

@inline function nbasis(lat::Lattice) :: Int
    length(lat.basis_vecs)
end

"Compute the volume of a unit cell."
@inline cell_volume(lat::Lattice) = abs(det(lat.lat_vecs))
"Compute the volume of the full simulation box."
@inline volume(lat::Lattice) = cell_volume(lat) * prod(lat.size)

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

"""
    CellType

An enumeration over the different types of 3D Bravais unit cells.
"""
@enum CellType begin
    triclinic
    monoclinic
    orthorhombic
    tetragonal
    # Rhombohedral is a special case. It is a lattice type (a=b=c, α=β=γ) but
    # not a spacegroup type. Trigonal space groups are conventionally described
    # using either hexagonal or rhombohedral lattices.
    rhombohedral
    hexagonal
    cubic
end

"""
    cell_type(lat_vecs::Mat3)

Infer the `CellType` of a unit cell from its lattice vectors, i.e. the columns
of `lat_vecs`. Report an error if lattice vectors are not in conventional form.
"""
function cell_type(lat_vecs::Mat3)
    a, b, c, α, β, γ = lattice_params(lat_vecs)

    if !(lat_vecs ≈ lattice_vectors(a, b, c, α, β, γ))
        error("Lattice vectors are not in conventional form. Consider using `lattice_vectors(a, b, c, α, β, γ)`.")
    end

    if a ≈ b ≈ c
        if α ≈ β ≈ γ ≈ 90
            return cubic
        elseif α ≈ β ≈ γ
            return rhombohedral
        end
    end

    if α ≈ β ≈ γ ≈ 90
        if a ≈ b
            return tetragonal
        elseif b ≈ c || c ≈ a
            error("Found a nonconventional tetragonal unit cell. Use `lattice_vectors(a, a, c, 90, 90, 90)` instead.")
        else
            return orthorhombic
        end
    end

    if (a ≈ b && α ≈ β ≈ 90 && (γ ≈ 60 || γ ≈ 120)) ||
       (b ≈ c && β ≈ γ ≈ 90 && (α ≈ 60 || α ≈ 120)) ||
       (c ≈ a && γ ≈ α ≈ 90 && (β ≈ 60 || β ≈ 120))
        if γ ≈ 120
            return hexagonal
        else
            error("Found a nonconventional hexagonal unit cell. Use `lattice_vectors(a, a, c, 90, 90, 120)` instead.")
        end
    end

    # Accept any of three possible permutations for monoclinic unit cell
    if α ≈ β ≈ 90 || β ≈ γ ≈ 90 || α ≈ γ ≈ 90
        return monoclinic
    end
    
    return triclinic
end

cell_type(lattice::Lattice{3}) = cell_type(lattice.lat_vecs)

"Return the standard cell convention for a given Hall number"
# Using the convention of spglib, listed at
# http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37
function cell_type(hall_number::Int)
    if 1 <= hall_number <= 2
        triclinic
    elseif 3 <= hall_number <= 107
        monoclinic
    elseif 108 <= hall_number <= 348
        orthorhombic
    elseif 349 <= hall_number <= 429
        tetragonal
    elseif 430 <= hall_number <= 461 # trigonal spacegroups
        # These special Hall numbers have "setting" R
        hall_number in [434, 437, 445, 451, 453, 459, 461] ? rhombohedral : hexagonal
    elseif 462 <= hall_number <= 488
        hexagonal
    elseif 489 <= hall_number <= 530
        cubic
    else
        error("Invalid Hall number $hall_number. Allowed range is 1..530")
    end
end

function is_trigonal_symmetry(hall_number::Int)
    return 430 <= hall_number <= 461
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
