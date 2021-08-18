module Symmetry

using Printf
using LinearAlgebra
using StaticArrays
using Parameters
import Spglib

import FastDipole: Vec3, Mat3, Lattice, lattice_params, lattice_vectors
export Crystal, Bond, canonical_bonds, print_bond_table

@enum CellType begin
    triclinic
    monoclinic
    orthorhombic
    tetragonal
    trigonal
    hexagonal
    cubic
end
const rhombohedral = trigonal

"Infer a 3D Bravais lattice cell type from its lattice vectors"
function cell_type(lat_vecs::Mat3)
    (a, b, c, α, β, γ) = lattice_params(lat_vecs)
    p = sortperm([a, b, c])
    a, b, c = (a, b, c)[p]
    α, β, γ = (α, β, γ)[p]
    if a ≈ b ≈ c
        if α ≈ β ≈ γ ≈ 90.
            cubic
        elseif α ≈ β ≈ γ
            trigonal
        end
    elseif a ≈ b
        if α ≈ β ≈ 90.
            if γ ≈ 90.
                tetragonal
            elseif γ ≈ 120.
                hexagonal
            end
        end
    elseif b ≈ c
        if β ≈ γ ≈ 90.
            if α ≈ 90.
                tetragonal
            elseif α ≈ 120.
                hexagonal
            end
        end
    elseif α ≈ β ≈ γ ≈ 90.
        orthorhombic
    elseif α ≈ β ≈ 90. || β ≈ γ ≈ 90. || α ≈ γ ≈ 90.
        monoclinic
    else
        triclinic
    end
end

"Return the standard cell convention for a given Hall number"
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
        trigonal
    elseif 462 <= hall_number <= 488
        hexagonal
    elseif 489 <= hall_number <= 530
        cubic
    else
        error("Invalid Hall number $hall_number. Allowed range is 1..530")
    end
end

function is_standard_form(lat_vecs::Mat3)
    lat_params = lattice_params(lat_vecs)
    conventional_lat_vecs = lattice_vectors(lat_params...)
    return lat_vecs ≈ conventional_lat_vecs
end

# A SymOp is composed of a rotation matrix and a translation vector. These are
# represented in fractional coordinates.
struct SymOp
    R::Mat3
    T::Vec3
end

"Holds all of the symmetry information about a crystal's unit cell"
struct Crystal{T}
    lat_vecs             :: Mat3             # Lattice vectors as columns
    positions            :: Vector{Vec3}     # Full set of atoms, fractional coords
    equiv_atoms          :: Vector{Int}      # Index to equivalent atom type
    species              :: Vector{T}        # Species for each atom
    symops               :: Vector{SymOp}    # Symmetry operations
    hall_number          :: Int              # Hall number
    symprec              :: Float64          # Tolerance to imperfections in symmetry
end

"Represents a bond between two sublattice sites, and displaced by some number of unit cells"
struct Bond{D}
    i :: Int
    j :: Int
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

function distance(lat::Lattice{D}, b::Bond{D}) where {D}
    return norm(lat.lat_vecs * b.n + (lat.basis_vecs[b.j] - lat.basis_vecs[b.i]))
end

function midpoint(b::BondRaw)
    # fractional coordinates
    return (b.r1 + b.r2) / 2
end

function Base.isapprox(s1::SymOp, s2::SymOp; atol)
    isapprox(s1.R, s2.R; atol) && isapprox(s1.T, s2.T; atol)
end


# Let Spglib infer the space group
function Crystal(lat_vecs::Mat3, positions::Vector{Vec3}, species::Vector{T}; symprec=1e-5) where {T}
    positions = wrap_to_unit_cell.(positions; symprec)

    cell = Spglib.Cell(lat_vecs, hcat(positions...), species)
    d = Spglib.get_dataset(cell, symprec)

    Rs = Mat3.(transpose.(eachslice(d.rotations, dims=3)))
    Ts = Vec3.(eachcol(d.translations))
    symops = map(SymOp, Rs, Ts)
    
    Crystal{T}(lat_vecs, positions, d.equivalent_atoms, species, symops, d.hall_number, symprec)
end

# Build Crystal using the space group denoted by a unique Hall number. The complete
# list is given at http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{T}, hall_number::Int; symprec=1e-5) where {T}
    latvec_cell_type = cell_type(lat_vecs)
    hall_cell_type = cell_type(hall_number)
    @assert latvec_cell_type == hall_cell_type "Hall number $hall_number requires cell type $hall_cell_type; received instead $latvec_cell_type."
    @assert is_standard_form(lat_vecs) "Lattice vectors must be in standard form. Consider using `lattice_vectors(a, b, c, α, β, γ)`."

    rotations, translations = Spglib.get_symmetry_from_database(hall_number)
    Rs = Mat3.(transpose.(eachslice(rotations, dims=3)))
    Ts = Vec3.(eachcol(translations))
    symops = map(SymOp, Rs, Ts)

    return Crystal(lat_vecs, base_positions, base_species, symops; hall_number, symprec)
end

# Make best effort to build Crystal from symbolic representation of spacegroup
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{T}, symbol::String; symprec=1e-5) where {T}
    # See "Complete list of space groups" at Seto's home page:
    # http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
    n_space_groups = 530

    crysts = Crystal{T}[]
    for hall_number in 1:n_space_groups
        sgt = Spglib.get_spacegroup_type(hall_number)
        if symbol in [sgt.hall_symbol, sgt.international, sgt.international_short, sgt.international_full]
            c = Crystal(lat_vecs, base_positions, base_species, hall_number; symprec)
            push!(crysts, c)
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
            hall_symbol = Spglib.get_spacegroup_type(hall_number).hall_symbol
            n_atoms     = length(c.positions)
            println("   Hall group '$hall_symbol' (number $hall_number), which generates $n_atoms atoms")
        end
        println()
        println("I will select Hall group number $(first(crysts).hall_number). You may wish to specify")
        println("an alternative Hall group in place of the symbol '$symbol'.")
        return first(crysts)
    end
end

"Build Crystal from explicit set of symmetry operations and a minimal set of positions "
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{T}, symops::Vector{SymOp}; hall_number=nothing, symprec=1e-5) where {T}
    # Spglib can map a Hall number to symops (via `get_symmetry_from_database`)
    # and can map symops to a Hall number (via `get_hall_number_from_symmetry`).
    # Unfortunately the round trip is not the identity function. For diamond
    # cubic (Hall numbers 525 and 526), I observed 526 -> symops -> 525. This
    # might be a bug in Spglib. For now, allow the caller to pass an explicit
    # `hall_number`. If `nothing`, then Spglib can infer the Hall number.
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
    species = T[]
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

    ret = Crystal{T}(lat_vecs, positions, equiv_atoms, species, symops, hall_number, symprec)
    validate(ret)
    return ret
end

"Filter sublattices of a Crystal by species, keeping the symmetry group of the original Crystal"
function subcrystal(cryst::Crystal{T}, species::T) :: Crystal{T} where {T}
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

    return Crystal{T}(cryst.lat_vecs, new_positions, new_equiv_atoms,
                      new_species, cryst.symops, cryst.hall_number, cryst.symprec)
end

"Filter sublattices by equivalency indices, keeping the symmetry group of the original Crystal"
function subcrystal(cryst::Crystal{T}, equiv_idxs::Vector{Int}) :: Crystal{T} where {T}
    new_positions = empty(cryst.positions)
    new_species = empty(cryst.species)
    new_equiv_atoms = empty(cryst.equiv_atoms)

    for (i, equiv_idx) in enumerate(equiv_idxs)
        subindexes = findall(isequal(equiv_idx), cryst.equiv_atoms)
        append!(new_positions, cryst.positions[subindexes])
        append!(new_species, cryst.species[subindexes])
        append!(new_equiv_atoms, fill(i, length(subindexes)))
    end

    return Crystal{T}(cryst.lat_vecs, new_positions, new_equiv_atoms,
                      new_species, cryst.symops, cryst.hall_number, cryst.symprec)
end

subcrystal(cryst::Crystal{T}, equiv_idx::Int) where {T} = subcrystal(cryst, [equiv_idx])

function Base.display(cryst::Crystal)
    printstyled("Crystal info\n"; bold=true, color=:underline)
    sgt = Spglib.get_spacegroup_type(cryst.hall_number)
    println("Hall group '$(sgt.hall_symbol)' (Hall number $(cryst.hall_number))")
    println("International symbol '$(sgt.international)' (number $(sgt.number))")
    println("Unit cell contains:")
    for s in unique(cryst.species)
        idxs = findall(==(s), cryst.species)
        n = length(idxs)
        uniq = length(unique(cryst.equiv_atoms[idxs]))
        println("    $n atoms of species '$s' ($uniq equivalence classes)")
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

# Constructors converting Lattice -> Crystal, and Crystal -> Lattice

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
    # Make the vertical axis 100x longer than the longest 2D lattice vector
    max_len = maximum(norm.(eachcol(lattice.lat_vecs)))
    L3[3, 3] = 100. * max_len
    L3 = SMatrix(L3)

    basis_coords = map(b -> SVector{3}((inv(L) * b)..., 0.), lattice.basis_vecs)

    Crystal(L3, basis_coords, lattice.species)
end

function Lattice(cryst::Crystal{String}, latsize) :: Lattice{3, 9, 4}
    basis_vecs = map(b -> cryst.lat_vecs * b, cryst.positions)
    Lattice{3}(cryst.lat_vecs, basis_vecs, cryst.species, latsize)
end


function validate(cryst::Crystal)
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
        @assert norm(R*R' - I) < 1e-12
    end

    # TODO: Check that space group is closed and that symops have inverse?

    # Check that symmetry elements are present in Spglib-inferred space group
    cryst′ = Crystal(cryst.lat_vecs, cryst.positions, cryst.species; cryst.symprec)
    for s in cryst.symops
        @assert any(cryst′.symops) do s′
            isapprox(s, s′; atol=cryst.symprec)
        end
    end
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


function find_symmetry_between_bonds(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    # Fail early if two bonds describe different real-space distances
    # (dimensionless error tolerance is measured relative to the minimum lattice
    # constant ℓ)
    ℓ = minimum(norm, eachcol(cryst.lat_vecs))
    d1 = distance(cryst, b1) / ℓ
    d2 = distance(cryst, b2) / ℓ
    if abs(d1-d2) > cryst.symprec
        return nothing
    end

    # Check whether any symmetry operation s maps b2 to b1
    for s = cryst.symops
        b2′ = transform(s, b2)

        # Check equivalence of b1 and b2′
        if (is_equivalent_by_translation(cryst, b1, b2′) ||
            is_equivalent_by_translation(cryst, b1, reverse(b2′)))
            return s
        end
    end

    # No symmetry was found
    return nothing
end

function is_equivalent_by_symmetry(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    return !isnothing(find_symmetry_between_bonds(cryst::Crystal, b1::BondRaw, b2::BondRaw))
end

function is_equivalent_by_symmetry(cryst::Crystal, b1::Bond{3}, b2::Bond{3})
    return is_equivalent_by_symmetry(cryst, BondRaw(cryst, b1), BondRaw(cryst, b2))
end

# List of all symmetries for bond, along with parity
function symmetries(cryst::Crystal, b::BondRaw)
    acc = Tuple{SymOp, Bool}[]
    for s = cryst.symops
        b′ = transform(s, b)
        # Check equivalence of b and b′
        if is_equivalent_by_translation(cryst, b, b′)
            push!(acc, (s, true))
        elseif is_equivalent_by_translation(cryst, b, reverse(b′))
            push!(acc, (s, false))
        end
    end
    return acc
end

function symmetries(cryst::Crystal, b::Bond{3})
    return symmetries(cryst, BondRaw(cryst, b))
end


function all_bonds_for_atom(cryst::Crystal, i::Int, max_dist::Float64)
    # be a little generous with the maximum distance
    max_dist += 4 * cryst.symprec

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
                
                # loop over all atoms within neighboring cryst
                for j in eachindex(cryst.positions)
                    b = Bond{3}(i, j, n)
                    if distance(cryst, b) <= max_dist
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

    # Favor indices where i < j
    score += 1e-2 * (b.i < b.j ? -1 : +1)

    return score
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

function _print_bond_table_header()
    header = "  i   j              n    dist     allowed J\n"
    printstyled(header; bold=true, color=:underline)
end

# Print a nice table of all possible bond classes to a maximum distance, and
#   the allowed interactions on the canonical bonds for each class.
function print_bond_table(cryst::Crystal, max_dist)
    _print_bond_table_header()

    canon_bonds = canonical_bonds(cryst, max_dist)
    dists = map(b->distance(cryst, b), canon_bonds)

    for (bond, dist) in zip(canon_bonds, dists)
        @unpack i, j, n = bond
        d = distance(cryst, bond)
        line = @sprintf "%3i %3i  [%3i,%3i,%3i]    %.2f     " i j n[1] n[2] n[3] d
        print(line)

        allowed_J_basis = basis_for_symmetry_allowed_couplings(cryst, bond)
        J_strings = _coupling_basis_strings(allowed_J_basis)
        max_len = maximum(length, J_strings)
        for i in 1:3
            print('|')
            for j in 1:3
                elem = J_strings[i, j]
                elem = repeat(' ', max_len-length(elem)) * elem
                print(elem * " ")
            end
            println('|')
            if i != 3
                print(repeat(' ', 35))
            end
        end
        println()
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
    for (s, parity) in symmetries(cryst, b)
        P = P * projector_for_symop(cryst, s, parity)
    end
    # Allowed coupling matrices J are simultaneously eigenvectors for all
    # projectors above, with eigenvalue 1.
    return P
end

# Check that a coupling matrix J is consistent with symmetries of a bond
function verify_coupling_matrix(cryst::Crystal, b::BondRaw, J::Mat3)
    for (s, parity) in symmetries(cryst, b)
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        @assert norm(R*J*R' - (parity ? J : J')) < 1e-12 "Specified J matrix not in allowed space!"
    end
end

function verify_coupling_matrix(cryst::Crystal, b::Bond, J::Mat3)
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
        for row = eachrow(A)
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

function allowed_J(cryst::Crystal, b::Bond{3}; digits=2, tol=1e-4)
    J_basis = basis_for_symmetry_allowed_couplings(cryst, b)
    _coupling_basis_strings(J_basis; digits=digits, tol=tol)
end


function all_symmetry_related_bonds(cryst::Crystal, b_ref::Bond{3})
    bs = Bond{3}[]
    for i in eachindex(cryst.positions)
        for b in all_bonds_for_atom(cryst, i, distance(cryst, b_ref))
            if is_equivalent_by_symmetry(cryst, b_ref, b)
                push!(bs, b)
            end
        end
    end
    return bs
end


function all_symmetry_related_interactions(cryst::Crystal, b_ref::Bond{3}, J_ref::Mat3)
    verify_coupling_matrix(cryst, BondRaw(cryst, b_ref), J_ref)

    bs = Bond{3}[]
    Js = Mat3[]

    for i in eachindex(cryst.positions)
        for b in all_bonds_for_atom(cryst, i, distance(cryst, b_ref))
            s = find_symmetry_between_bonds(cryst, BondRaw(cryst, b_ref), BondRaw(cryst, b))
            if !isnothing(s)
                R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
                J = R * J_ref * R'
                push!(bs, b)
                push!(Js, J)
            end
        end
    end

    return (bs, Js)
end

end