"""
Structs and functions for implementing various pair interaction energies
 and fields, given a specific lattice to operate on.
Upon creation of a SpinSystem, all pair interactions get converted into their
 corresponding type here.
"""

"""
    BondTable{T}

Stores all of the bonds and associated data (e.g. exchange scalars/matrices) of type T
contained within a single pair interaction densely in memory, with quick access to bonds/Ts
on any sublattice.
"""
struct BondTable{T}
    bonds         :: Vector{Bond}      # All bonds, sorted on first basis index
    culled_bonds  :: Vector{Bond}      # Culled to minimal 1/2 bonds, sorted on first basis index
    data          :: Vector{T}         # Data T associated with each bond in `bonds`
    culled_data   :: Vector{T}         # Culled data T associated with each bond in `bonds`
    basis_indices :: Vector{Int}       # Indices to first bond of each basis index in `bonds`
end

"""Construct a `BondTable` given a `Crystal`, and a `Bond` along with its associated
    exchange matrix. Performs the symmetry transformations to obtain all other bonds
    and exchanges matrices.
"""
function BondTable(crystal::Crystal, bond::Bond, J::Mat3)
    bonds = Vector{Bond}()
    Js = Vector{Mat3}()
    basis_indices = [1]

    for i in 1:nbasis(crystal)
        (symbonds, symJs) = all_symmetry_related_couplings_for_atom(crystal, i, bond, J)
        append!(bonds, symbonds)
        append!(Js, symJs)
        push!(basis_indices, last(basis_indices) + length(symbonds))
    end

    # Compute the minimal set of the bonds/Js needed to specify all bonds
    # Bonds with i ↔ j and n ↔ -n are equivalent: only keep one of each pair
    culled_bonds = Vector{Bond}()
    culled_Js = Vector{Mat3}()
    for (bond, J) in zip(bonds, Js)
        neg_bond = Bond(bond.j, bond.i, -bond.n)
        if !(neg_bond in culled_bonds)
            push!(culled_bonds, bond)
            push!(culled_Js, J)
        end
    end

    BondTable{Mat3}(bonds, culled_bonds, Js, culled_Js, basis_indices)
end

function Base.show(io::IO, ::MIME"text/plain", bondtable::BondTable{T}) where {T}
    print(io, "$(length(bondtable))-element, $(length(bondtable.basis_indices)-1)-basis BondTable{$T}")
end

Base.length(bondtable::BondTable) = length(bondtable.bonds)
nbasis(bondtable::BondTable) = length(bondtable.basis_indices) - 1
## These functions generate iterators producing (bond, data) for various subsets of bonds ##
all_bonds(bondtable::BondTable) = zip(bondtable.bonds, bondtable.data)
culled_bonds(bondtable::BondTable) = zip(bondtable.culled_bonds, bondtable.culled_data)
function sublat_bonds(bondtable::BondTable, i::Int)
    @unpack bonds, data, basis_indices = bondtable
    @boundscheck checkindex(Bool, 1:length(basis_indices)-1, i) ? nothing : throw(BoundsError(bondtable, i))
    zip(bonds[basis_indices[i]:basis_indices[i+1]-1], data[basis_indices[i]:basis_indices[i+1]-1])
end

# This method `map`'s only the data within a BondTable, not the bonds themselves.
function Base.map(f, bondtable::BondTable)
    @unpack bonds, culled_bonds, data, culled_data, basis_indices = bondtable
    new_data = map(f, data)
    new_culled_data = map(f, culled_data)
    return BondTable{eltype(new_data)}(
        bonds, culled_bonds, new_data, new_culled_data, basis_indices
    )
end

abstract type AbstractPairIntCPU <: AbstractInteractionCPU end

"""
    HeisenbergCPU <: AbstractPairIntCPU

Implements an exchange interaction which is proportional to
the identity matrix.
"""
struct HeisenbergCPU <: AbstractPairIntCPU
    bondtable    :: BondTable{Float64}    # Bonds store effective J's = S_i J S_j (all identical)
    label        :: String
end

"""
    DiagonalCouplingCPU <: AbstractPairIntCPU

Implements an exchange interaction where matrices on all bonds
are diagonal.
"""
struct DiagonalCouplingCPU <: AbstractPairIntCPU
    bondtable     :: BondTable{Vec3}    # Bonds store diagonal of effective J's = S_i J S_j
    label         :: String
end

"""
    GeneralCouplingCPU <: AbstractPairIntCPU

Implements the most generalized interaction, where matrices on
all bonds are full 3x3 matrices which vary bond-to-bond.
"""
struct GeneralCouplingCPU <: AbstractPairIntCPU
    bondtable    :: BondTable{Mat3}    # Bonds store effective J's = S_i J S_j
    label        :: String
end

# Helper functions producing predicates checking if a matrix is approximately
# Heisenberg or diagonal
isheisen(J, tol) = isapprox(diagm(fill(J[1,1], 3)), J; atol=tol)
isdiag(J, tol) = isapprox(diagm(diag(J)), J; atol=tol)
isheisen(tol) = Base.Fix2(isheisen, tol)
isdiag(tol) = Base.Fix2(isdiag, tol)

# Figures out the correct maximally-efficient backend type for a quadratic interaction, and perform
#  all of the symmetry propagation.
function convert_quadratic(int::QuadraticInteraction, cryst::Crystal, site_infos::Vector{SiteInfo}; tol=1e-6)
    @unpack J, bond, label = int

    bondtable = BondTable(cryst, bond, J)

    if all(isheisen(tol), bondtable.culled_data)
        # Pull out the Heisenberg scalar from each exchange matrix
        bondtable = map(J->J[1,1], bondtable)
        return HeisenbergCPU(bondtable, label)
    elseif all(isdiag(tol), bondtable.culled_data)
        # Pull out the diagonal of each exchange matrix
        bondtable = map(J->diag(J), bondtable)
        return DiagonalCouplingCPU(bondtable, label)
    else
        return GeneralCouplingCPU(bondtable, label)
    end
end

function energy(dipoles::Array{Vec3}, heisenberg::HeisenbergCPU)
    bondtable = heisenberg.bondtable
    effJ = first(bondtable.data)
    E = 0.0

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, _) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = dipoles[cell, i]
            sⱼ = dipoles[offset(cell, n, latsize), j]
            E += sᵢ ⋅ sⱼ
        end
    end
    return effJ * E
end

function energy(dipoles::Array{Vec3}, diag_coup::DiagonalCouplingCPU)
    bondtable = diag_coup.bondtable
    E = 0.0

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = dipoles[cell, i]
            sⱼ = dipoles[offset(cell, n, latsize), j]
            E += (J .* sᵢ) ⋅ sⱼ
        end
    end
    return E
end

function energy(dipoles::Array{Vec3}, gen_coup::GeneralCouplingCPU)
    bondtable = gen_coup.bondtable
    E = 0.0

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = dipoles[cell, i]
            sⱼ = dipoles[offset(cell, n, latsize), j]
            E += dot(sᵢ, J, sⱼ)
        end
    end
    return E
end

"Accumulates the local -∇ℋ coming from Heisenberg couplings into `B`"
function _accum_neggrad!(B::Array{Vec3}, dipoles::Array{Vec3}, heisen::HeisenbergCPU)
    bondtable = heisen.bondtable
    effJ = first(bondtable.data)

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, _) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            offsetcell = offset(cell, n, latsize)
            B[cell, i] = B[cell, i] - effJ * dipoles[offsetcell, j]
            B[offsetcell, j] = B[offsetcell, j] - effJ * dipoles[cell, i]
        end
    end
end

function _neggrad(dipoles::Array{Vec3}, heisen::HeisenbergCPU, i)
    bondtable = heisen.bondtable
    B = Vec3(0,0,0)
    J = first(bondtable.data)
    _, site = splitidx(i)

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, _) in sublat_bonds(bondtable, site)
        (; j, n) = bond
        offsetcell = offset(cell, n, latsize)
        B -= J * dipoles[offsetcell, j]
    end

    return B
end

"Accumulates the local -∇ℋ coming from diagonal couplings into `B`"
function _accum_neggrad!(B::Array{Vec3}, dipoles::Array{Vec3}, diag_coup::DiagonalCouplingCPU)
    bondtable = diag_coup.bondtable

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            offsetcell = offset(cell, n, latsize)
            B[cell, i] = B[cell, i] - J .* dipoles[offsetcell, j]
            B[offsetcell, j] = B[offsetcell, j] - J .* dipoles[cell, i]
        end
    end
end

function _neggrad(dipoles::Array{Vec3}, diag_coup::DiagonalCouplingCPU, i)
    bondtable = diag_coup.bondtable
    B = Vec3(0,0,0)
    cell, site = splitidx(i)

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, J) ∈ sublat_bonds(bondtable, site)
        (; j, n) = bond
        offsetcell = offset(cell, n, latsize)
        B -= J .* dipoles[offsetcell, j]
    end

    return B
end

"Accumulates the local -∇ℋ coming from general couplings into `B`"
function _accum_neggrad!(B::Array{Vec3}, dipoles::Array{Vec3}, gen_coup::GeneralCouplingCPU)
    bondtable = gen_coup.bondtable

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            offsetcell = offset(cell, n, latsize)
            B[cell, i] = B[cell, i] - J * dipoles[offsetcell, j]
            B[offsetcell, j] = B[offsetcell, j] - J * dipoles[cell, i]
        end
    end
end


function _neggrad(dipoles::Array{Vec3}, gen_coup::GeneralCouplingCPU, idx)
    bondtable = gen_coup.bondtable
    B = Vec3(0,0,0)
    cell, site = splitidx(idx)

    latsize = size(dipoles)[1:3]
    @inbounds for (bond, J) ∈ sublat_bonds(bondtable, site)
        (; j, n) = bond
        offsetcell = offset(cell, n, latsize)
        B -= J * dipoles[offsetcell, j]
    end

    return B
end