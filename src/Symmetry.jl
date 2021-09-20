module Symmetry

using Printf
using LinearAlgebra
using StaticArrays
using Parameters

import Spglib

export Crystal, Bond, canonical_bonds, print_bond_table
export lattice_params, lattice_vectors, nbasis, cell_volume, CellType, cell_type

const Vec3 = SVector{3, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}


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

function is_standard_form(lat_vecs::Mat3)
    lat_params = lattice_params(lat_vecs)
    conventional_lat_vecs = lattice_vectors(lat_params...)
    return lat_vecs ≈ conventional_lat_vecs
end

# Return true if lattice vectors are compatible with a monoclinic space group
# "setting" for a Hall number.
function is_compatible_monoclinic_cell(lat_vecs, hall_number)
    @assert cell_type(hall_number) == monoclinic

    # Special handling of monoclinic space groups. There are three possible
    # conventions for the unit cell, depending on which of α, β, or γ is
    # special.
    _, _, _, α, β, γ = lattice_params(lat_vecs)
    choice = Spglib.get_spacegroup_type(hall_number).choice
    x = first(replace(choice, "-" => ""))
    if x == 'a'
        return β≈90 && γ≈90
    elseif x == 'b'
        return α≈90 && γ≈90
    elseif x == 'c'
        return α≈90 && β≈90
    else
        error()
    end
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
    elseif 430 <= hall_number <= 461
        # The trigonal space groups require either rhombohedral or hexagonal
        # cells. The Hall numbers below have "setting" R.
        hall_number in [434, 437, 445, 451, 453, 459, 461] ? rhombohedral : hexagonal
    elseif 462 <= hall_number <= 488
        hexagonal
    elseif 489 <= hall_number <= 530
        cubic
    else
        error("Invalid Hall number $hall_number. Allowed range is 1..530")
    end
end

function all_compatible_cells(cell::CellType)
    if cell == triclinic
        [triclinic, monoclinic, orthorhombic, tetragonal, rhombohedral, hexagonal, cubic]
    elseif cell == monoclinic
        [monoclinic, orthorhombic, tetragonal, hexagonal, cubic]
    elseif cell == orthorhombic
        [orthorhombic, tetragonal, cubic]
    elseif cell == tetragonal
        [tetragonal, cubic]
    elseif cell == rhombohedral
        [rhombohedral]
    elseif cell == hexagonal
        [hexagonal]
    elseif cell == cubic
        [cubic]
    else
        error()
    end
end

function is_trigonal_symmetry(hall_number::Int)
    return 430 <= hall_number <= 461
end

"""
    SymOp

Defines a symmetry operation belonging to a 3D space group, operating on fractional coordinates.
"""
struct SymOp
    R::Mat3
    T::Vec3
end

"""
    Crystal

A type holding all geometry and symmetry information needed to represent
 a three-dimensional crystal.
"""
struct Crystal
    lat_vecs             :: Mat3             # Lattice vectors as columns
    positions            :: Vector{Vec3}     # Full set of atoms, fractional coords
    equiv_atoms          :: Vector{Int}      # Index to equivalent atom type
    species              :: Vector{String}   # Species for each atom
    symops               :: Vector{SymOp}    # Symmetry operations
    hall_number          :: Int              # Hall number
    symprec              :: Float64          # Tolerance to imperfections in symmetry
end

nbasis(cryst::Crystal) = length(cryst.positions)
cell_volume(cryst::Crystal) = abs(det(lat.lat_vecs))
lattice_params(cryst::Crystal) = lattice_params(cryst.lat_vecs)
lattice_vectors(cryst::Crystal) = cryst.lat_vecs

"""
    Bond{D}

Represents a class of bond between pairs of sites in a `D`-dimensional crystal,
 separated from each other by a certain vector.
"""
struct Bond{D}
    "The sublattice index of the first site."
    i :: Int
    "The sublattice index of the second site."
    j :: Int
    "The displacement vector between sites, in units of lattice vectors."
    n :: SVector{D, Int}
end

"Represents a bond expressed as two fractional coordinates"
struct BondRaw
    r1::SVector{3, Float64}
    r2::SVector{3, Float64}
end

"Convenience constructor for a 2D bond"
Bond2D(i, j, n::Vector{Int}) = Bond{2}(i, j, n)

"Convenience constructor for a 3D bond"
Bond3D(i, j, n::Vector{Int}) = Bond{3}(i, j, n)


function Bond(cryst::Crystal, b::BondRaw)
    i1 = position_to_index(cryst, b.r1)
    i2 = position_to_index(cryst, b.r2)
    r1 = cryst.positions[i1]
    r2 = cryst.positions[i2]
    n = round.(Int, (b.r2-b.r1) - (r2-r1))
    return Bond{3}(i1, i2, n)
end

function BondRaw(cryst::Crystal, b::Bond{3})
    return BondRaw(cryst.positions[b.i], cryst.positions[b.j]+b.n)
end

function is_same_position(x, y; symprec=1e-5)
    return norm(rem.(x-y, 1, RoundNearest)) < symprec
end

function position_to_index(cryst::Crystal, r::Vec3)
    return findfirst(r′ -> is_same_position(r, r′; symprec=cryst.symprec), cryst.positions)
end

# Wrap each coordinate of position r into the range [0,1). To account for finite
# precision, wrap 1-ϵ to -ϵ, where ϵ=symprec is a tolerance parameter.
function wrap_to_unit_cell(r::Vec3; symprec=1e-5)
    return @. mod(r+symprec, 1) - symprec
end

function distance(cryst::Crystal, b::BondRaw)
    return norm(cryst.lat_vecs * (b.r1 - b.r2))
end

function distance(cryst::Crystal, b::Bond{3})
    return distance(cryst, BondRaw(cryst, b))
end

function Base.isapprox(s1::SymOp, s2::SymOp; atol)
    isapprox(s1.R, s2.R; atol) && isapprox(s1.T, s2.T; atol)
end


"""
    Crystal(lat_vecs::Mat3, positions::Vector{Vec3}, species::Vector{String}; symprec=1e-5)

Construct a `Crystal` using explicit geometry information, with all symmetry information
automatically inferred. `positions` should be a list of site positions (in fractional
coordinates) within the unit cell defined by lattice vectors which are the columns of `lat_vecs`.
"""
function Crystal(lat_vecs::Mat3, positions::Vector{Vec3}, species::Vector{String}; symprec=1e-5)
    positions = wrap_to_unit_cell.(positions; symprec)

    cell = Spglib.Cell(lat_vecs, hcat(positions...), species)
    d = Spglib.get_dataset(cell, symprec)
    equiv_atoms = d.equivalent_atoms
    Rs = Mat3.(transpose.(eachslice(d.rotations, dims=3)))
    Ts = Vec3.(eachcol(d.translations))
    symops = map(SymOp, Rs, Ts)
    
    # Sort atoms so that they are contiguous in equivalence classes
    p = sortperm(equiv_atoms)
    positions .= positions[p]
    species .= species[p]

    # Base constructor
    ret = Crystal(lat_vecs, positions, equiv_atoms, species, symops, d.hall_number, symprec)
    validate(ret)
    return ret
end

# Build Crystal using the space group denoted by a unique Hall number. The complete
# list is given at http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
# TODO: Make hall_number a named parameter. But then we have to avoid a conflict
# with the constructor above. Maybe use unique names for every constructor?
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{String}, hall_number::Int; symprec=1e-5)
    cell = cell_type(lat_vecs)
    hall_cell = cell_type(hall_number)
    allowed_cells = all_compatible_cells(hall_cell)
    @assert cell in allowed_cells "Hall number $hall_number requires a $hall_cell cell, but found $cell."

    if hall_cell == monoclinic
        is_compatible = is_compatible_monoclinic_cell(lat_vecs, hall_number)
        @assert is_compatible "Lattice vectors define a monoclinic cell that is incompatible with Hall number $hall_number."
    end

    rotations, translations = Spglib.get_symmetry_from_database(hall_number)
    Rs = Mat3.(transpose.(eachslice(rotations, dims=3)))
    Ts = Vec3.(eachcol(translations))
    symops = map(SymOp, Rs, Ts)

    # Use symops to fill in symmetry-related atom positions
    return Crystal(lat_vecs, base_positions, base_species, symops; hall_number, symprec)
end

# Make best effort to build Crystal from symbolic representation of spacegroup
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{String}, symbol::String; symprec=1e-5)
    # See "Complete list of space groups" at Seto's home page:
    # http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
    n_space_groups = 530

    crysts = Crystal[]
    for hall_number in 1:n_space_groups
        sgt = Spglib.get_spacegroup_type(hall_number)

        if (replace(symbol, " "=>"") == sgt.international_short || 
            symbol in [sgt.hall_symbol, sgt.international, sgt.international_full, "I:"*string(sgt.number)])

            # Some Hall numbers may be incompatible with unit cell of provided
            # lattice vectors; skip them.
            is_compatible = true

            cell = cell_type(lat_vecs)
            hall_cell = cell_type(hall_number)
            allowed_cells = all_compatible_cells(hall_cell)

            # Special handling of trigonal space groups
            if is_trigonal_symmetry(hall_number)
                # Trigonal symmetry must have either hexagonal or rhombohedral
                # cell, according to the Hall number.
                is_latvecs_valid = cell in [rhombohedral, hexagonal]
                @assert is_latvecs_valid "Symbol $symbol requires a rhomobohedral or hexagonal cell, but found $cell."
                is_compatible = cell in allowed_cells
            else
                # For all other symmetry types, there is a unique cell for each Hall number
                @assert cell in allowed_cells "Symbol $symbol requires a $hall_cell cell, but found $cell."
            end

            if hall_cell == monoclinic
                is_compatible = is_compatible_monoclinic_cell(lat_vecs, hall_number)
            end

            if is_compatible
                # Build crystal using symops for hall_number from database
                c = Crystal(lat_vecs, base_positions, base_species, hall_number; symprec)
                push!(crysts, c)
            end
        end
    end

    if length(crysts) == 0
        error("Could not find symbol '$symbol' in database.")
    elseif length(crysts) == 1
        return first(crysts)
    else
        sort!(crysts, by=c->length(c.positions))

        println("Warning, the symbol '$symbol' is ambiguous. It could refer to:")
        for c in crysts
            hall_number = c.hall_number
            hm_symbol = Spglib.get_spacegroup_type(hall_number).international
            n_atoms = length(c.positions)
            println("   HM symbol '$hm_symbol' (Hall number $hall_number), which generates $n_atoms atoms")
        end
        println()
        println("Selecting Hall number $(first(crysts).hall_number). You may wish to specify")
        println("an alternative Hall number in place of the symbol '$symbol'.")
        return first(crysts)
    end
end

"Build Crystal from explicit set of symmetry operations and a minimal set of positions "
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{String}, symops::Vector{SymOp}; hall_number=nothing, symprec=1e-5)
    # Spglib can map a Hall number to symops (via `get_symmetry_from_database`)
    # and can map symops to a Hall number (via `get_hall_number_from_symmetry`).
    # Unfortunately the round trip is not the identity function:
    #   https://github.com/spglib/spglib/issues/132
    # As a workaround, allow the caller to pass an explicit `hall_number`. If
    # `nothing`, then Spglib can infer the Hall number.
    if isnothing(hall_number)
        rotation = zeros(3, 3, length(symops))
        translation = zeros(3, length(symops))
        for (i, s) = enumerate(symops)
            rotation[:, :, i] = s.R'
            translation[:, i] = s.T
        end
        hall_number = Int(Spglib.get_hall_number_from_symmetry(rotation, translation, length(symops)))
    end
    
    positions = Vec3[]
    species = String[]
    equiv_atoms = Int[]
    
    for i = eachindex(base_positions)
        for s = symops
            x = wrap_to_unit_cell(transform(s, base_positions[i]); symprec)

            idx = findfirst(y -> is_same_position(x, y; symprec), positions)
            if isnothing(idx)
                push!(positions, x)
                push!(species, base_species[i])
                push!(equiv_atoms, i)
            else
                j = equiv_atoms[idx]
                if i != j
                    error("Base positions $(base_positions[i]) and $(base_positions[j]) are symmetry equivalent.")
                end
            end
        end
    end

    # Check that symops are present in Spglib-inferred space group
    cryst′ = Crystal(lat_vecs, positions, species; symprec)
    for s in symops
        @assert any(cryst′.symops) do s′
            isapprox(s, s′; atol=symprec)
        end
    end

    # Call base constructor
    ret = Crystal(lat_vecs, positions, equiv_atoms, species, symops, hall_number, symprec)
    validate(ret)
    return ret
end

"""
    subcrystal(cryst, species) :: Crystal

Filter sublattices of a `Crystal` by species, keeping the symmetry group of the original `Crystal`.
"""
function subcrystal(cryst::Crystal, species::String) :: Crystal
    if !in(species, cryst.species)
        error("Species string '$species' is not present in crystal.")
    end
    subindexes = findall(isequal(species), cryst.species)
    new_positions = cryst.positions[subindexes]
    new_equiv_atoms = cryst.equiv_atoms[subindexes]
    # Reduce all of the equivalent atom indexes to count again from 1, 2,...
    unique_equiv = unique(new_equiv_atoms)
    for (i, unique) in enumerate(unique_equiv)
        equiv_sites = findall(isequal(unique), new_equiv_atoms)
        new_equiv_atoms[equiv_sites] .= i
    end
    new_species = fill(species, length(subindexes))

    return Crystal(cryst.lat_vecs, new_positions, new_equiv_atoms,
                   new_species, cryst.symops, cryst.hall_number, cryst.symprec)
end

"""
    subcrystal(cryst, equiv_idxs) :: Crystal

Filter sublattices of a `Crystal` by a list of indexes into `cryst.equiv_atoms`,
 keeping the symmetry group of the original `Crystal`.
"""
function subcrystal(cryst::Crystal, equiv_idxs::Vector{Int}) :: Crystal
    for equiv_idx in equiv_idxs
        if !in(equiv_idx, cryst.equiv_atoms)
            error("Equivalent index '$equiv_idx' is not present in crystal.")
        end
    end
    new_positions = empty(cryst.positions)
    new_species = empty(cryst.species)
    new_equiv_atoms = empty(cryst.equiv_atoms)

    for (i, equiv_idx) in enumerate(equiv_idxs)
        subindexes = findall(isequal(equiv_idx), cryst.equiv_atoms)
        append!(new_positions, cryst.positions[subindexes])
        append!(new_species, cryst.species[subindexes])
        append!(new_equiv_atoms, fill(i, length(subindexes)))
    end

    return Crystal(cryst.lat_vecs, new_positions, new_equiv_atoms,
                   new_species, cryst.symops, cryst.hall_number, cryst.symprec)
end

subcrystal(cryst::Crystal, equiv_idx::Int) = subcrystal(cryst, [equiv_idx])

function Base.display(cryst::Crystal)
    printstyled("Crystal info\n"; bold=true, color=:underline)
    sgt = Spglib.get_spacegroup_type(cryst.hall_number)
    # println("Hall group '$(sgt.hall_symbol)' (Hall number $(cryst.hall_number))")
    println("H-M space group '$(sgt.international)' ($(sgt.number))")

    if is_standard_form(cryst.lat_vecs)
        (a, b, c, α, β, γ) = lattice_params(cryst.lat_vecs)
        @printf "Lattice params a=%.4g, b=%.4g, c=%.4g, α=%.4g°, β=%.4g°, γ=%.4g°\n" a b c α β γ
    else
        println("Lattice vectors:")
        for a in eachcol(cryst.lat_vecs)
            @printf "   [%.4g %.4g %.4g]\n" a[1] a[2] a[3]
        end
    end

    println("Atoms:")
    for i in eachindex(cryst.positions)
        print("   $i. ")
        if length(unique(cryst.species)) > 1
            print("Species '$(cryst.species[i])', ")
        end
        if length(unique(cryst.equiv_atoms)) > length(unique(cryst.species))
            print("Class $(cryst.equiv_atoms[i]), ")
        end
        r = cryst.positions[i]
        @printf "Coords [%.4g, %.4g, %.4g]\n" r[1] r[2] r[3]
    end
end


function transform(s::SymOp, r::Vec3)
    return s.R*r + s.T
end

function transform(s::SymOp, b::BondRaw)
    return BondRaw(transform(s, b.r1), transform(s, b.r2))
end

function transform(cryst::Crystal, s::SymOp, b::Bond{3})
    return Bond(cryst, transform(s, BondRaw(cryst, b)))
end

function Base.reverse(b::BondRaw)
    return BondRaw(b.r2, b.r1)
end



function validate(cryst::Crystal)
    # Atoms must be sorted by equivalence class
    sortperm(cryst.equiv_atoms) == eachindex(cryst.equiv_atoms)

    # Equivalent atoms must have the same species
    for i in eachindex(cryst.positions)
        for j in eachindex(cryst.positions)
            if cryst.equiv_atoms[i] == cryst.equiv_atoms[j]
                @assert cryst.species[i] == cryst.species[j]
            end
        end
    end

    # Check symmetry rotations are orthogonal, up to periodicity
    for s in cryst.symops
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        # Due to possible imperfections in the lattice vectors, only require
        # that R is approximately orthogonal
        @assert norm(R*R' - I) < cryst.symprec "Lattice vectors and symmetry operations are incompatible."
    end

    # TODO: Check that space group is closed and that symops have inverse?
end


function is_equivalent_by_translation(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    # Displacements between the two bonds
    D1 = b2.r1 - b1.r1
    D2 = b2.r2 - b1.r2
    # Round components of D1 to nearest integers
    n = round.(D1, RoundNearest)
    # If both n ≈ D1 and n ≈ D2, then the bonds are equivalent by translation
    return norm(n - D1) < cryst.symprec && norm(n - D2) < cryst.symprec
end


function is_equivalent_by_symmetry(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    return !isempty(symmetries_between_bonds(cryst::Crystal, b1::BondRaw, b2::BondRaw))
end

function is_equivalent_by_symmetry(cryst::Crystal, b1::Bond{3}, b2::Bond{3})
    return is_equivalent_by_symmetry(cryst, BondRaw(cryst, b1), BondRaw(cryst, b2))
end

# Generate list of all symmetries that transform b2 into b1, along with parity
function symmetries_between_bonds(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    # Fail early if two bonds describe different real-space distances
    # (dimensionless error tolerance is measured relative to the minimum lattice
    # constant ℓ)
    if b1 != b2
        ℓ = minimum(norm, eachcol(cryst.lat_vecs))
        d1 = distance(cryst, b1) / ℓ
        d2 = distance(cryst, b2) / ℓ
        if abs(d1-d2) > cryst.symprec
            return Tuple{SymOp, Bool}[]
        end
    end

    function bonds_equiv(s)
        b2′ = transform(s, b2)
        if is_equivalent_by_translation(cryst, b1, b2′)
            (s, true)
        elseif is_equivalent_by_translation(cryst, b1, reverse(b2′))
            (s, false)
        else
            nothing
        end
    end

    ret = (bonds_equiv(s) for s in cryst.symops)
    return Iterators.filter(!isnothing, ret)
end



"Returns all bonds in `cryst` for which `bond.i == i`"
function all_bonds_for_atom(cryst::Crystal, i::Int, max_dist; min_dist=0.0)
    # be a little generous with the minimum and maximum distances
    ℓ = minimum(norm, eachcol(cryst.lat_vecs))
    max_dist += 4 * cryst.symprec * ℓ
    min_dist -= 4 * cryst.symprec * ℓ

    # columns are the reciprocal vectors
    recip_vecs = 2π * inv(cryst.lat_vecs)'

    # box_lengths[i] represents the perpendicular distance between two parallel
    # boundary planes spanned by lattice vectors a_j and a_k (where indices j
    # and k differ from i)
    box_lengths = [a⋅b/norm(b) for (a,b) = zip(eachcol(cryst.lat_vecs), eachcol(recip_vecs))]
    n_max = round.(Int, max_dist ./ box_lengths, RoundUp)

    bonds = Bond{3}[]

    # loop over neighboring cells
    for n1 in -n_max[1]:n_max[1]
        for n2 in -n_max[2]:n_max[2]
            for n3 in -n_max[3]:n_max[3]
                n = SVector(n1, n2, n3)
                
                # loop over all atoms within neighboring cell
                for j in eachindex(cryst.positions)
                    b = Bond{3}(i, j, n)
                    if min_dist <= distance(cryst, b) <= max_dist
                        push!(bonds, b)
                    end
                end
            end
        end
    end

    return bonds
end


# Calculate score for a bond. Lower would be preferred.
function _score_bond(cryst::Crystal, b)
    # Favor bonds with fewer nonzero elements in basis matrices J
    Js = basis_for_symmetry_allowed_couplings(cryst, b)
    nnz = [count(abs.(J) .> 1e-12) for J in Js]
    score = sum(nnz)

    # Favor bonds with smaller unit cell displacements. Positive
    # displacements are slightly favored over negative displacements.
    # Displacements in x are slightly favored over y, etc.
    score += norm((b.n .- 0.1) .* [0.07, 0.08, 0.09])

    # Favor smaller indices and indices where i < j
    score += 1e-2 * (b.i + b.j) + 1e-2 * (b.i < b.j ? -1 : +1)

    return score
end


function bond_multiplicity(cryst::Crystal, b::Bond{3})
    return length(all_symmetry_related_bonds_for_atom(cryst, b.i, b))
end


"Produces a list of 'canonical' bonds that belong to different symmetry equivalence classes."
function canonical_bonds(cryst::Crystal, max_dist)
    # List of distinct atom types
    atom_types = unique(cryst.equiv_atoms)
    
    # Atom indices, one for each equivalence class
    canon_atoms = [findfirst(isequal(t), cryst.equiv_atoms) for t in atom_types]
    
    # Bonds, one for each equivalent class
    cbonds = Bond[]
    for i in canon_atoms
        for b in all_bonds_for_atom(cryst, i, max_dist)
            if !any(is_equivalent_by_symmetry(cryst, b, b′) for b′ in cbonds)
                push!(cbonds, b)
            end
        end
    end

    # Sort by distance
    sort!(cbonds, by=b->distance(cryst, b))

    # Replace each canonical bond by the "best" equivalent bond
    return map(cbonds) do cb
        # Find full set of symmetry equivalent bonds
        equiv_bonds = unique([transform(cryst, s, cb) for s in cryst.symops])
        # Take the bond with lowest score
        scores = [_score_bond(cryst, b) for b in equiv_bonds]
        return equiv_bonds[findmin(scores)[2]]
    end
end

# For each element of `bonds`, return an index into `canonical_bonds` that gives
# the equivalent bond
function equivalent_bond_indices(cryst::Crystal, canonical_bonds, bonds)
    map(bonds) do b
        findfirst(canonical_bonds) do cb
             is_equivalent_by_symmetry(cryst, b, cb)
        end
    end
end


function print_allowed_exchange(name::String, allowed_J_basis::Vector{Mat3})
    J_strings = _coupling_basis_strings(allowed_J_basis)
    max_len = maximum(length, J_strings)
    for i in 1:3
        print(i == 1 ? name : repeat(' ', length(name)))
        print('|')
        for j in 1:3
            elem = J_strings[i, j]
            elem = repeat(' ', max_len-length(elem)) * elem
            print(elem * " ")
        end
        println('|')
    end
end

function print_bond(cryst::Crystal, b::Bond{3})
    ri = cryst.positions[b.i]
    rj = cryst.positions[b.j] + b.n

    @printf "Bond{3}(%d, %d, [%d, %d, %d])\n" b.i b.j b.n[1] b.n[2] b.n[3]

    if b.i == b.j && iszero(b.n)
        # On site interaction
        if length(unique(cryst.species)) == 1
            @printf "At fractional coordinates [%.4g, %.4g, %.4g]\n" ri[1] ri[2] ri[3]
        else
            @printf "'%s' at fractional coordinates [%.4g, %.4g, %.4g]\n" cryst.species[b.i] ri[1] ri[2] ri[3]
        end
        allowed_J_basis = basis_for_symmetry_allowed_couplings(cryst, b)
        allowed_J_basis_sym = filter(J -> J ≈ J', allowed_J_basis)
        print_allowed_exchange("Allowed on-site anisotropy: ", allowed_J_basis_sym)
        # Antisymmetric terms are relevant to g-tensor. Report these separately.
        if length(allowed_J_basis) > length(allowed_J_basis_sym)
            print_allowed_exchange("Allowed effective g-tensor: ", allowed_J_basis)
        end
    else
        @printf "Distance %.4g, multiplicity %i\n" distance(cryst, b) bond_multiplicity(cryst, b)
        if length(unique(cryst.species)) == 1
            @printf "Connects [%.4g, %.4g, %.4g] to [%.4g, %.4g, %.4g]\n" ri[1] ri[2] ri[3] rj[1] rj[2] rj[3]
        else
            @printf "Connects '%s' at [%.4g, %.4g, %.4g] to '%s' at [%.4g, %.4g, %.4g]\n" cryst.species[b.i] ri[1] ri[2] ri[3] cryst.species[b.j] rj[1] rj[2] rj[3]
        end
        allowed_J_basis = basis_for_symmetry_allowed_couplings(cryst, b)
        print_allowed_exchange("Allowed exchange matrix: ", allowed_J_basis)
    end
    println()
end


"""
    print_bond_table(cryst::Crystal, max_dist)

Pretty-prints a table of all symmetry classes of bonds present in `cryst`, up
 to a maximum bond length of `max_dist` in absolute coordinates. For each class,
 a single "canonical" bond is shown, with `i`, `j` being the two sublattices
 connected, and `n` being the displacement vector in units of the lattice vectors.
"""
function print_bond_table(cryst::Crystal, max_dist)
    for b in canonical_bonds(cryst, max_dist)
        print_bond(cryst, b)
    end
end


# Interaction matrix J should be invariant under symmetry operations. For a
# symop involving the orthogonal matrix R, we require either `R J Rᵀ = J` or 
# `R J Rᵀ = Jᵀ`, depending on the parity of the bond symmetry.
#
# Let F denote the linear operator such that `F J = R J Rᵀ - J` or
# `F J = R J Rᵀ - Jᵀ`. We can view F as a 9x9 matrix. In the first case,
# F = R⊗R - I. In the second case, we should replace I by the suitable
# transpose operation, 'transpose_op_3x3'.
#
# Any linear combination "vectors" in the null space of F, i.e. `F J = 0`,
# satisfy the original constraint . The singular value decomposition
#
#     F = U Σ Vᵀ
#
# can be used to produce an orthogonal basis for this null space. It is spanned
# by columns v of V corresponding to the zero singular values. The space spanned
# by v equivalently represented as a projection matrix,
#
#     P_R = v vᵀ
#
# When there are multiple constraints (R1, R2, ... Rn) we should take the
# intersection of the spanned spaces of P_R1, ... P_Rn. To calculate this
# intersection, form the product of the projectors:
#
#     P = P_R1 P_R2 ... P_Rn
# 
# The allowable J values correspond to the eigenvectors of P with eigenvalue 1.
# An orthogonal basis for this space can again be calculated with an SVD.

const transpose_op_3x3 = [
    1 0 0  0 0 0  0 0 0
    0 0 0  1 0 0  0 0 0
    0 0 0  0 0 0  1 0 0

    0 1 0  0 0 0  0 0 0
    0 0 0  0 1 0  0 0 0
    0 0 0  0 0 0  0 1 0

    0 0 1  0 0 0  0 0 0
    0 0 0  0 0 1  0 0 0
    0 0 0  0 0 0  0 0 1
]


# Returns a projection operator P that maps to zero any symmetry-unallowed
# coupling matrix J. The space spanned by the eigenvectors of P with eigenvalue
# 1 represents the allowed coupling matrices J.
function projector_for_symop(cryst::Crystal, s::SymOp, parity::Bool)
    # Cartesian-space rotation operator corresponding to `s`
    R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)

    # Constraint is modeled as `F J = 0`
    F = kron(R, R) - (parity ? I : transpose_op_3x3)

    # Orthogonal column vectors that span the null space of F
    v = nullspace(F; atol=1e-12)

    # Projector onto the null space of F
    P = v * v'
    return P
end


# Return an operator P that implicitly gives the space of symmetry allowed
# coupling matrices for bond b. Specifically, x is an allowed coupling if and
# only if it is an eigenvector of P with eigenvalue 1, i.e., `P x = x`.
function symmetry_allowed_couplings_operator(cryst::Crystal, b::BondRaw)
    P = I
    for (s, parity) in symmetries_between_bonds(cryst, b, b)
        P = P * projector_for_symop(cryst, s, parity)
    end
    # Allowed coupling matrices J are simultaneously eigenvectors for all
    # projectors above, with eigenvalue 1.
    return P
end

# Check that a coupling matrix J is consistent with symmetries of a bond
function verify_coupling_matrix(cryst::Crystal, b::BondRaw, J::Mat3)
    for (s, parity) in symmetries_between_bonds(cryst, b, b)
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        @assert norm(R*J*R' - (parity ? J : J')) < 1e-12 "Specified J matrix not in allowed space!"
    end
end

function verify_coupling_matrix(cryst::Crystal, b::Bond{3}, J::Mat3)
    verify_coupling_matrix(cryst, BondRaw(cryst, b), J)
end


# Orthonormal basis of 3x3 antisymmetric matrices
const asym_basis = begin
    b = [[ 0  1  0
          -1  0  0
           0  0  0]/√2,
         [ 0  0  1
           0  0  0
          -1  0  0]/√2,
         [ 0  0  0
           0  0  1
           0 -1  0]/√2]
    SMatrix{9, 3, Float64}(hcat(reshape.(b, 9)...))
end

# Orthonormal basis of 3x3 symmetric matrices
const sym_basis = begin
    b = [[0 1 0
          1 0 0
          0 0 0]/√2,
         [0 0 1
          0 0 0
          1 0 0]/√2,
         [0 0 0
          0 0 1
          0 1 0]/√2,
         diagm([1, 0, 0]),
         diagm([0, 1, 0]),
         diagm([0, 0, 1])]
    SMatrix{9, 6, Float64}(hcat(reshape.(b, 9)...))
end

@assert sym_basis * sym_basis' + asym_basis * asym_basis' ≈ I


# Given an m×n matrix A with empty nullspace, linearly combine the n columns to
# make them sparser.
function sparsify_columns(A; atol)
    if size(A, 2) <= 1
        return A
    else
        # By assumption, the n columns of A are linearly independent
        @assert isempty(nullspace(A; atol))
        # Since row rank equals column rank, it should be possible to find n
        # linearly independent rows in A
        indep_rows = Vector{Float64}[]
        for row in eachrow(A)
            # Add `row` to list if linearly independent with existing ones
            if isempty(nullspace(hcat(indep_rows..., row); atol))
                push!(indep_rows, row)
            end
        end
        return A * inv(hcat(indep_rows...))'
    end
end

# Find a convenient basis for the symmetry allowed couplings on bond b
function basis_for_symmetry_allowed_couplings(cryst::Crystal, b::BondRaw)
    P = symmetry_allowed_couplings_operator(cryst, b)

    acc = SVector{9, Float64}[]

    # If any "reference" basis vectors are eigenvalues of P with eigenvalue 1,
    # use them as outputs, and remove them from P
    for x in eachcol(hcat(asym_basis, sym_basis))
        if isapprox(P*x, x; atol=1e-12)
            push!(acc, x)
            P = P * (I - x*x')
        end
    end

    # Any solution to the original symmetry constraints `R J Rᵀ = J` or `R J Rᵀ
    # = Jᵀ` decomposes into purely symmetric/antisymmetric solutions. Therefore
    # we can pick a basis that separates into symmetric and antisymmetric parts.
    # We will do so by decomposing P. By construction, P = P_sym+P_asym.
    P_sym  = P *  sym_basis *  sym_basis'
    P_asym = P * asym_basis * asym_basis'

    # Search for eigenvectors of P_sym with eigenvalue 1. These provide an
    # orthonormal basis for symmetric couplings.
    v = nullspace(P_sym-I; atol=1e-12)
    v = sparsify_columns(v; atol=1e-12)
    append!(acc, eachcol(v))

    # Similarly for antisymmetric couplings.
    v = nullspace(P_asym-I; atol=1e-12)
    v = sparsify_columns(v; atol=1e-12)
    append!(acc, eachcol(v))

    return map(acc) do x
        # Normalize each basis vector so that its maximum component is 1. The
        # shift by ϵ avoids unnecessary sign change in the case where the
        # maximum magnitude values of x appear symmetrically as ±c for some c.
        ϵ = 1e-12
        _, i = findmax(abs.(x.+ϵ))
        x = x / x[i]

        # Reinterpret as 3x3 matrix
        x = Mat3(reshape(x, 3, 3))
        
        # Double check that x indeed satifies the necessary symmetries
        verify_coupling_matrix(cryst, b, x)

        return x
    end
end

function basis_for_symmetry_allowed_couplings(cryst::Crystal, b::Bond{3})
    return basis_for_symmetry_allowed_couplings(cryst, BondRaw(cryst, b))
end

"Removes trailing zeros after a decimal place, and returns empty for 1.0"
function _strip_decimal_string(str)
    decimal_idx = findfirst('.', str)
    chop_idx = length(str)
    if !isnothing(decimal_idx)
        for i in length(str):-1:decimal_idx-1
            chop_idx = i
            cur_char = str[chop_idx]

            if (cur_char != '0' && cur_char != '.')
                break
            end
        end
    end
    # If we're left with just a 1 or -1, remove the 1
    if (chop_idx == decimal_idx-1 && str[chop_idx] == '1')
        chop_idx -= 1
    end
    return str[1:chop_idx]
end

"Converts a list of basis elements for a J matrix into a nice string summary"
function _coupling_basis_strings(coup_basis; digits=2, tol=1e-4) :: Matrix{String}
    J = fill("", size(coup_basis[1])...)
    for (letter, basis_mat) in zip('A':'Z', coup_basis)
        for idx in eachindex(basis_mat)
            coeff = basis_mat[idx]
            if abs(coeff) > tol
                coeff = round(coeff; digits=digits)
                if J[idx] == ""
                    float_str = @sprintf "%.4f" coeff
                else
                    float_str = @sprintf "%+.4f" coeff
                end
                float_str = _strip_decimal_string(float_str)
                J[idx] *= float_str * letter
            end
        end
    end
    for idx in eachindex(J)
        if J[idx] == ""
            J[idx] = "0"
        end
    end
    return J
end

"""
    allowed_J(cryst::Crystal, b::Bond{3}; digits=2, tol=1e-4)

Given a bond `b`, returns a `Matrix{String}` representing the allowed
 form of a bilinear exchange interaction matrix on this bond, given the
 symmetry constraints of `cryst`.
"""
function allowed_J(cryst::Crystal, b::Bond{3}; digits=2, tol=1e-4)
    J_basis = basis_for_symmetry_allowed_couplings(cryst, b)
    _coupling_basis_strings(J_basis; digits=digits, tol=tol)
end


function all_symmetry_related_bonds_for_atom(cryst::Crystal, i::Int, b_ref::Bond{3})
    bs = Bond{3}[]
    dist = distance(cryst, b_ref)
    for b in all_bonds_for_atom(cryst, i, dist; min_dist=dist)
        if is_equivalent_by_symmetry(cryst, b_ref, b)
            push!(bs, b)
        end
    end
    return bs
end

"""
    all_symmetry_related_bonds(cryst::Crystal, b_ref::Bond{3}) :: Vector{Bond{3}}

Construct a list of all bonds which are symmetry-equivalent to the reference bond `b_ref`
 within `cryst.
"""
function all_symmetry_related_bonds(cryst::Crystal, b_ref::Bond{3})
    bs = Bond{3}[]
    for i in eachindex(cryst.positions)
        append!(bs, all_symmetry_related_bonds_for_atom(cryst, i, b_ref))
    end
    bs
end

function all_symmetry_related_interactions_for_atom(cryst::Crystal, i::Int, b_ref::Bond{3}, J_ref::Mat3)
    verify_coupling_matrix(cryst, BondRaw(cryst, b_ref), J_ref)

    bs = Bond{3}[]
    Js = Mat3[]

    for b in all_symmetry_related_bonds_for_atom(cryst, i, b_ref)
        push!(bs, b)
        (s, parity) = first(symmetries_between_bonds(cryst, BondRaw(cryst, b_ref), BondRaw(cryst, b)))
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        push!(Js, R * (parity ? J_ref : J_ref') * R')
    end

    return (bs, Js)
end

"""
    all_symmetry_related_interactions(cryst::Crystal, b_ref::Bond{3}, J_ref::Mat3) :: Tuple{Vector{Bond3}, Vector{Mat3}}

Given a reference bond `b_ref` and exchange matrix `J_ref` on that bond, construct lists of all
 symmetry-equivalent bonds and their respective transformed exchange matrices.
"""
function all_symmetry_related_interactions(cryst::Crystal, b_ref::Bond{3}, J_ref::Mat3)
    bs = Bond{3}[]
    Js = Mat3[]

    for i in eachindex(cryst.positions)
        (bs_i, Js_i) = all_symmetry_related_interactions_for_atom(cryst, i, b_ref, J_ref)
        append!(bs, bs_i)
        append!(Js, Js_i)
    end

    return (bs, Js)
end

end