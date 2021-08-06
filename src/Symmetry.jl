module Symmetry

using Printf
using LinearAlgebra
using StaticArrays
using Parameters
import Spglib

import FastDipole: Vec3, Mat3, Lattice
export Cell, Bond, canonical_bonds, print_bond_table

# Rotational and translational parts in fractional coordinates
struct SymOp
    R::Mat3
    T::Vec3
end

"Holds all of the symmetry information about a crystal's unit cell"
struct Cell{T}
    lattice              :: Mat3             # Lattice vectors in columns
    positions            :: Vector{Vec3}     # Full set of atoms, fractional coords
    equiv_atoms          :: Vector{Int}      # Index to equivalent atom type
    species              :: Vector{T}        # Species for each atom
    symops               :: Vector{SymOp}    # Symmetry operations
    international_number :: Int              # International number for space group
    international_symbol :: String           # International symbol (H-M)
    hall_number          :: Int              # Hall number (-1 if unknown)
    hall_symbol          :: String           # Hall symbol ("" if unknown)
    symprec              :: Float64          # Tolerance to imperfections in symmetry
end

struct Bond{D}
    i :: Int
    j :: Int
    n :: SVector{D, Int}
end

# Like Bond3, but with fractional positions instead of atom indices
struct BondRaw
    r1::SVector{3, Float64}
    r2::SVector{3, Float64}
end


function BondRaw(cell::Cell, b::Bond{3})
    return BondRaw(cell.positions[b.i], cell.positions[b.j]+b.n)
end

function distance(cell::Cell, b::BondRaw)
    return norm(cell.lattice * (b.r1 - b.r2))
end

function distance(cell::Cell, b::Bond{3})
    return distance(cell, BondRaw(cell, b))
end

function Base.isapprox(s1::SymOp, s2::SymOp; atol)
    isapprox(s1.R, s2.R; atol) && isapprox(s1.T, s2.T; atol)
end


# Let Spglib infer the space group
function Cell(lattice::Mat3, positions::Vector{Vec3}, species::Vector{T}; symprec=1e-5) where {T}
    cell = Spglib.Cell(lattice, hcat(positions...), species)
    d = Spglib.get_dataset(cell, symprec)

    Rs = Mat3.(eachslice(d.rotations, dims=3))'
    Ts = Vec3.(eachcol(d.translations))
    symops = map(SymOp, Rs, Ts)
    
    international_number = d.spacegroup_number
    Cell{T}(Mat3(lattice), Vec3.(positions), d.equivalent_atoms, species, symops, international_number, d.international_symbol, d.hall_number, d.hall_symbol, symprec)
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


# Build Cell from explicit set of symmetry operations
function Cell(lattice, base_positions, species::Vector{T}, symops::Vector{SymOp}; symprec=1e-5) where {T}
    is_same_position(x, y) = norm(lattice * rem.(x-y, 1, RoundNearest)) < symprec

    lattice = Mat3(lattice)
    base_positions = Vec3.(base_positions)

    positions = Vec3[]
    equiv_atoms = Int[]
    
    for i = eachindex(base_positions)
        for s = symops
            x = transform(s, base_positions[i])

            idx = findfirst(y -> is_same_position(x, y), positions)
            if isnothing(idx)
                push!(positions, x)
                push!(equiv_atoms, i)
            else
                j = equiv_atoms[idx]
                if i != j
                    error("Base positions $(base_positions[i]) and $(base_positions[j]) are symmetry equivalent.")
                end
            end
        end
    end

    species = [species[i] for i=equiv_atoms]

    ret = Cell{T}(lattice, positions, equiv_atoms, species, symops, -1, "", -1, "", symprec)
    validate(ret)
    return ret
end

# Constructors converting Lattice -> Cell, and Cell -> Lattice

function Cell(lattice::Lattice{3, 9, 4})
    L = lattice.lat_vecs
    # Convert absolute basis positions to fractional coordinates
    basis_coords = map(b -> inv(L) * b, lattice.basis_vecs)
    Cell(L, basis_coords, lattice.basis_species)
end

function Cell(lattice::Lattice{2, 4, 3})
    L = lattice.lat_vecs

    # Expand the lattice to 3D, with a very long c axis
    L3 = @MMatrix zeros(3, 3)
    L3[1:2, 1:2] = L
    # Make the vertical axis 10x longer than the longest 2D lattice vector
    max_len = maximum(norm.(eachcol(lattice.lat_vecs)))
    L3[3, 3] = 10. * max_len
    L3 = SMatrix(L3)

    basis_coords = map(b -> SVector{3}((inv(L) * b)..., 0.), lattice.basis_vecs)

    Cell(L3, basis_coords, lattice.basis_species)
end

function Lattice(cell::Cell{String}, latsize) :: Lattice{3, 9, 4}
    basis_vecs = map(b -> cell.lattice * b, cell.positions)
    Lattice{3, 9, 4}(cell.lattice, basis_vecs, cell.species, latsize)
end


function validate(cell::Cell)
    # Equivalent atoms must have the same species
    for i = eachindex(cell.positions)
        for j = eachindex(cell.positions)
            if cell.equiv_atoms[i] == cell.equiv_atoms[j]
                @assert cell.species[i] == cell.species[j]
            end
        end
    end

    # Check symmetry rotations are orthogonal
    for s = cell.symops
        R = cell.lattice * s.R * inv(cell.lattice)
        @assert norm(R*R' - I) < 1e-12
    end

    # TODO: Check that space group is closed and that symops have inverse?

    # Check that symmetry elements are present in Spglib-inferred space group
    cell′ = Cell(cell.lattice, cell.positions, cell.species; cell.symprec)
    for s = cell.symops
        @assert any(cell′.symops) do s′
            isapprox(s, s′; atol=cell.symprec)
        end
    end
end


function is_equivalent_by_translation(cell::Cell, b1::BondRaw, b2::BondRaw)
    # Displacements between the two bonds
    D1 = b2.r1 - b1.r1
    D2 = b2.r2 - b1.r2
    # Round components of D1 to nearest integers
    n = round.(D1, RoundNearest)
    # If both n ≈ D1 and n ≈ D2, then the bonds are equivalent by translation
    return norm(n - D1) < cell.symprec && norm(n - D2) < cell.symprec
end


function find_symmetry_between_bonds(cell::Cell, b1::BondRaw, b2::BondRaw)
    # Fail early if two bonds describe different real-space distances
    d1 = norm(cell.lattice*(b1.r2 - b1.r1))
    d2 = norm(cell.lattice*(b2.r2 - b2.r1))
    if abs(d1 - d2) > cell.symprec
        return nothing
    end

    # Check whether any symmetry operation s maps b2 to b1
    for s = cell.symops
        b2′ = transform(s, b2)

        # Check equivalence of l1 and l2′
        if (is_equivalent_by_translation(cell, b1, b2′) ||
            is_equivalent_by_translation(cell, b1, reverse(b2′)))
            return s
        end
    end

    # No symmetry was found
    return nothing
end

function is_equivalent_by_symmetry(cell::Cell, b1::BondRaw, b2::BondRaw)
    return !isnothing(find_symmetry_between_bonds(cell::Cell, b1::BondRaw, b2::BondRaw))
end

function is_equivalent_by_symmetry(cell::Cell, b1::Bond{3}, b2::Bond{3})
    return is_equivalent_by_symmetry(cell, BondRaw(cell, b1), BondRaw(cell, b2))
end


# List of all symmetries for bond, along with parity
function symmetries(cell::Cell, b::BondRaw)
    acc = Tuple{SymOp, Bool}[]
    for s = cell.symops
        b′ = transform(s, b)
        # Check equivalence of l and l′
        if is_equivalent_by_translation(cell, b, b′)
            push!(acc, (s, true))
        elseif is_equivalent_by_translation(cell, b, reverse(b′))
            push!(acc, (s, false))
        end
    end
    return acc
end

# function symmetries(cell::Cell, b::Bond{3})
#     return symmetries(cell, BondRaw(cell, b))
# end


function all_bonds_for_atom(cell::Cell, i::Int, max_dist::Float64)
    # be a little generous with the maximum distance
    max_dist += 4cell.symprec

    # columns are the reciprocal vectors
    recip_lattice = 2π*inv(cell.lattice)'

    # box_lengths[i] represents the perpendicular distance between two parallel
    # boundary planes spanned by lattice vectors a_j and a_k (where indices j
    # and k differ from i)
    box_lengths = [a⋅b/norm(b) for (a,b) = zip(eachcol(cell.lattice), eachcol(recip_lattice))]
    n_max = round.(Int, max_dist ./ box_lengths, RoundUp)

    bonds = Bond{3}[]

    # loop over neighboring cells
    for n1 in -n_max[1]:n_max[1]
        for n2 in -n_max[2]:n_max[2]
            for n3 in -n_max[3]:n_max[3]
                n = SVector(n1, n2, n3)
                
                # loop over all atoms within neighboring cell
                for j in eachindex(cell.positions)
                    b = Bond{3}(i, j, n)
                    if distance(cell, b) <= max_dist
                        push!(bonds, b)
                    end
                end
            end
        end
    end

    return bonds
end


function canonical_bonds(cell::Cell, max_dist)
    # list of distinct atom types
    atom_types = unique(cell.equiv_atoms)
    
    # canonical indices for atoms of each type
    canon_atoms = [findfirst(isequal(t), cell.equiv_atoms) for t = atom_types]
    
    cbonds = BondRaw[]

    for i in canon_atoms
        for b in all_bonds_for_atom(cell, i, max_dist)
            b = BondRaw(cell, b)
            if !any(is_equivalent_by_symmetry(cell, b, b′) for b′=cbonds)
                push!(cbonds, b)
            end
        end
    end

    sort!(cbonds, by=b->distance(cell, b))
    return cbonds
end

# For each element of `bonds`, return an index into `canonical_bonds` that gives
# the equivalent bond
function equivalent_bond_indices(cell::Cell, canonical_bonds, bonds)
    map(bonds) do b
        findfirst(canonical_bonds) do cb
             is_equivalent_by_symmetry(cell, b, cb)
        end
    end
end

function _print_bond_table_header()
    header = "dist  class   i   j   n               allowed J\n"
    printstyled(header; bold=true, color=:underline)
end


# Print a nice table of all possible bonds up to a maximum distance, and
#   the allowed interactions on these bonds.
function print_bond_table(cell::Cell, max_dist)
    _print_bond_table_header()

    canon_links = canonical_bonds(cell, max_dist)
    lengths = map(l->distance(cell, l), canon_links)
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

    for (link, dist, class) in zip(canon_links, dists, classes)
        all_bonds = all_symmetry_related_bonds(cell, link)
        canon_bond = first(all_bonds)
        @unpack i, j, n = canon_bond
        line = @sprintf "%4i  %5i %3i %3i  [%3i,%3i,%3i]    " dist class i j n[1] n[2] n[3]
        println(line)

        allowed_J_basis = basis_for_symmetry_allowed_couplings(cell, link)
        # 
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
function projector_for_symop(cell::Cell, s::SymOp, parity::Bool)
    # Cartesian-space rotation operator corresponding to `s`
    R = cell.lattice * s.R * inv(cell.lattice)

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
function symmetry_allowed_couplings_operator(cell::Cell, b::BondRaw)
    P = I
    for (s, parity) in symmetries(cell, b)
        P = P * projector_for_symop(cell, s, parity)
    end
    # Allowed coupling matrices J are simultaneously eigenvectors for all
    # projectors above, with eigenvalue 1.
    return P
end

# Check that a coupling matrix J is consistent with symmetries of a bond
function verify_coupling_matrix(cell::Cell, b::BondRaw, J::Mat3)
    for (s, parity) in symmetries(cell, b)
        R = cell.lattice * s.R * inv(cell.lattice)
        @assert norm(R*J*R' - (parity ? J : J')) < 1e-12
    end
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
function basis_for_symmetry_allowed_couplings(cell::Cell, b::BondRaw)
    P = symmetry_allowed_couplings_operator(cell, b)

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
        verify_coupling_matrix(cell, b, x)

        return x
    end
end

"Converts a list of basis elements for a J matrix into a nice string summary"
function pretty_print_coupling_basis(coup_basis)
    nothing
end


function all_symmetry_related_bonds(cell::Cell, b_ref::BondRaw)
    bs = Bond{3}[]
    for i in eachindex(cell.positions)
        for b in all_bonds_for_atom(cell, i, distance(cell, b_ref))
            if is_equivalent_by_symmetry(cell, b_ref, BondRaw(cell, b))
                push!(bs, b)
            end
        end
    end
    return bs
end


function all_symmetry_related_interactions(cell::Cell, b_ref::BondRaw, J_ref::Mat3)
    verify_coupling_matrix(cell, b_ref, J_ref)

    bs = Bond{3}[]
    Js = Mat3[]

    for i in eachindex(cell.positions)
        for b in all_bonds_for_atom(cell, i, distance(cell, b_ref))
            s = find_symmetry_between_bonds(cell, b_ref, BondRaw(cell, b))
            if !isnothing(s)
                R = cell.lattice * s.R * inv(cell.lattice)
                J = R * J_ref * R'
                push!(bs, b)
                push!(Js, J)
            end
        end
    end

    return (bs, Js)
end

end