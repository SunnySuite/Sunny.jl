module Symmetry

using Printf
using LinearAlgebra
using StaticArrays
using Parameters
import Spglib

import FastDipole: Vec3, Mat3, Lattice
export Crystal, Bond, canonical_bonds, print_bond_table

# Rotational and translational parts in fractional coordinates
struct SymOp
    R::Mat3
    T::Vec3
end

"Holds all of the symmetry information about a crystal's unit cryst"
struct Crystal{T}
    lat_vecs             :: Mat3             # Lattice vectors as columns
    positions            :: Vector{Vec3}     # Full set of atoms, fractional coords
    equiv_atoms          :: Vector{Int}      # Index to equivalent atom type
    species              :: Vector{T}        # Species for each atom
    symops               :: Vector{SymOp}    # Symmetry operations
    hall_number          :: Int              # Hall number
    symprec              :: Float64          # Tolerance to imperfections in symmetry
end

struct Bond{D}
    i :: Int
    j :: Int
    n :: SVector{D, Int}
end

# Warning: This constructor is most likely type-unstable. Don't use
#  in production code.
function Bond(i, j, n)
    D = length(n)
    n = SVector{D, Float64}(n)
    Bond{D}(i, j, n)
end

# Like Bond3, but with fractional positions instead of atom indices
struct BondRaw
    r1::SVector{3, Float64}
    r2::SVector{3, Float64}
end


function BondRaw(cryst::Crystal, b::Bond{3})
    return BondRaw(cryst.positions[b.i], cryst.positions[b.j]+b.n)
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

function Base.isapprox(s1::SymOp, s2::SymOp; atol)
    isapprox(s1.R, s2.R; atol) && isapprox(s1.T, s2.T; atol)
end


# Let Spglib infer the space group
function Crystal(lat_vecs::Mat3, positions::Vector{Vec3}, species::Vector{T}; symprec=1e-5) where {T}
    cryst = Spglib.Cell(lat_vecs, hcat(positions...), species)
    d = Spglib.get_dataset(cryst, symprec)

    Rs = Mat3.(transpose.(eachslice(d.rotations, dims=3)))
    Ts = Vec3.(eachcol(d.translations))
    symops = map(SymOp, Rs, Ts)
    
    Crystal{T}(lat_vecs, positions, d.equivalent_atoms, species, symops, d.hall_number, symprec)
end

# Build cryst using the space group denoted by a unique Hall number. The complete
# list is given at http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{T}, hall_number::Int; symprec=1e-5) where {T}
    is_same_position(x, y) = norm(lat_vecs * rem.(x-y, 1, RoundNearest)) < symprec

    rotations, translations = Spglib.get_symmetry_from_database(hall_number)
    Rs = Mat3.(transpose.(eachslice(rotations, dims=3)))
    Ts = Vec3.(eachcol(translations))
    symops = map(SymOp, Rs, Ts)

    positions = Vec3[]
    species = T[]
    equiv_atoms = Int[]
    
    for i = eachindex(base_positions)
        for s = symops
            x = transform(s, base_positions[i])

            idx = findfirst(y -> is_same_position(x, y), positions)
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

# Make best effort to build cryst from symbolic representation of spacegroup
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{T}, symbol::String; symprec=1e-5) where {T}
    # See "Complete list of space groups" at Seto's home page:
    # http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
    n_space_groups = 530

    cells = Crystal{T}[]
    for hall_number in 1:n_space_groups
        sgt = Spglib.get_spacegroup_type(hall_number)
        if symbol in [sgt.hall_symbol, sgt.international, sgt.international_short, sgt.international_full]
            c = Crystal(lat_vecs, base_positions, base_species, hall_number; symprec)
            push!(cells, c)
        end
    end

    if length(cells) == 0
        error("Could not find symbol '$symbol' in database.")
    elseif length(cells) == 1
        return first(cells)
    else
        sort!(cells, by = c->length(c.positions))

        println("Warning, the symbol '$symbol' is ambiguous. It could refer to:")
        for c in cells
            hall_number = c.hall_number
            hall_symbol = Spglib.get_spacegroup_type(hall_number).hall_symbol
            n_atoms     = length(c.positions)
            println("   Hall group '$hall_symbol' (number $hall_number), which generates $n_atoms atoms")
        end
        println()
        println("I will select Hall group number $(first(cells).hall_number). You may wish to specify")
        println("an alternative Hall group in place of the symbol '$symbol'.")
        return first(cells)
    end
end

# Build Crystal from explicit set of symmetry operations
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_species::Vector{T}, symops::Vector{SymOp}; symprec=1e-5) where {T}
    rotation = zeros(3, 3, length(symops))
    translation = zeros(3, length(symops))
    for (i, s) = enumerate(symops)
        rotation[:, :, i] = s.R'
        translation[:, i] = s.T
    end
    hall_number = Int(Spglib.get_hall_number_from_symmetry(rotation, translation, length(symops)))
    return Crystal(lat_vecs, base_positions, base_species, hall_number; symprec)
end

function display(cryst::Crystal)
    printstyled("Crystal info\n"; bold=true, color=:underline)
    sgt = Spglib.get_spacegroup_type(cryst.hall_number)
    println("Hall group '$(sgt.hall_symbol)' (number $(cryst.hall_number))")
    println("International symbol '$(sgt.international)'")
    println("Unit cryst contains:")
    for s in unique(cryst.species)
        n = count(==(s), cryst.species)
        println("    $n atoms of species '$s'")
    end
end


function transform(s::SymOp, r::Vec3)
    return s.R*r + s.T
end

function transform(s::SymOp, b::BondRaw)
    return BondRaw(transform(s, b.r1), transform(s, b.r2))
end

function reverse(b::BondRaw)
    return BondRaw(b.r2, b.r1)
end

# Constructors converting Lattice -> Crystal, and Crystal -> Lattice

function Crystal(lattice::Lattice{3, 9, 4})
    L = lattice.lat_vecs
    # Convert absolute basis positions to fractional coordinates
    basis_coords = map(b -> inv(L) * b, lattice.basis_vecs)
    Crystal(L, basis_coords, lattice.basis_species)
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

    Crystal(L3, basis_coords, lattice.basis_species)
end

function Lattice(cryst::Crystal{String}, latsize) :: Lattice{3, 9, 4}
    basis_vecs = map(b -> cryst.lat_vecs * b, cryst.positions)
    Lattice{3, 9, 4}(cryst.lat_vecs, basis_vecs, cryst.species, latsize)
end


function validate(cryst::Crystal)
    # Equivalent atoms must have the same species
    for i = eachindex(cryst.positions)
        for j = eachindex(cryst.positions)
            if cryst.equiv_atoms[i] == cryst.equiv_atoms[j]
                @assert cryst.species[i] == cryst.species[j]
            end
        end
    end

    # Check symmetry rotations are orthogonal
    for s = cryst.symops
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        @assert norm(R*R' - I) < 1e-12
    end

    # TODO: Check that space group is closed and that symops have inverse?

    # Check that symmetry elements are present in Spglib-inferred space group
    cryst′ = Crystal(cryst.lat_vecs, cryst.positions, cryst.species; cryst.symprec)
    for s = cryst.symops
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
    d1 = norm(cryst.lat_vecs*(b1.r2 - b1.r1))
    d2 = norm(cryst.lat_vecs*(b2.r2 - b2.r1))
    if abs(d1 - d2) > cryst.symprec
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

# function symmetries(cryst::Crystal, b::Bond{3})
#     return symmetries(cryst, BondRaw(cryst, b))
# end


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

"Produces a list of symmetry equivalence classes, specified by an arbitrary bond."
function canonical_bonds(cryst::Crystal, max_dist)
    # list of distinct atom types
    atom_types = unique(cryst.equiv_atoms)
    
    # canonical indices for atoms of each type
    canon_atoms = [findfirst(isequal(t), cryst.equiv_atoms) for t = atom_types]
    
    cbonds = Bond[]

    for i in canon_atoms
        for b in all_bonds_for_atom(cryst, i, max_dist)
            if !any(is_equivalent_by_symmetry(cryst, b, b′) for b′=cbonds)
                push!(cbonds, b)
            end
        end
    end

    sort!(cbonds, by=b->distance(cryst, b))
    return cbonds
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
    header = "dist  class   i   j   n               allowed J\n"
    printstyled(header; bold=true, color=:underline)
end


# Print a nice table of all possible bonds up to a maximum distance, and
#   the allowed interactions on these bonds.
function print_bond_table(cryst::Crystal, max_dist)
    _print_bond_table_header()

    canon_bonds = canonical_bonds(cryst, max_dist)
    lengths = map(b->distance(cryst, b), canon_bonds)
    # Integer orderings of the length of all canon bonds
    dists = indexin(lengths, unique(lengths)) .- 1

    # Enumerate symmetry classes of all canon bonds by seeing
    #   where dists changes. Kind of ugly.
    classes = [1]
    for (i, dist) in enumerate(dists[2:end])
        if dists[i] == dist
            push!(classes, classes[i] + 1)
        else
            push!(classes, 1)
        end
    end

    for (bond, dist, class) in zip(canon_bonds, dists, classes)
        @unpack i, j, n = bond
        line = @sprintf "%4i  %5i %3i %3i  [%3i,%3i,%3i]    " dist class i j n[1] n[2] n[3]
        print(line)

        allowed_J_basis = basis_for_symmetry_allowed_couplings(cryst, bond)
        J_strings = coupling_basis_strings(allowed_J_basis)
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
                print(repeat(' ', 38))
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
    θ = svd(F)
    v = θ.Vt[findall(x -> abs(x) < 1e-12, θ.S), :]'

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
    θ = svd(P_sym-I)
    v = θ.Vt[findall(λ -> abs(λ) < 1e-12, θ.S), :]'
    append!(acc, eachcol(v))

    # Similarly for antisymmetric couplings.
    θ = svd(P_asym-I)
    v = θ.Vt[findall(λ -> abs(λ) < 1e-12, θ.S), :]'
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
function coupling_basis_strings(coup_basis; digits=2, tol=1e-4) :: Matrix{String}
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