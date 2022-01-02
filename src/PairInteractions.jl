"""
Structs and functions for implementing various pair interaction energies
 and fields, given a specific lattice to operate on.
Upon creation of a SpinSystem, all pair interactions get converted into their
 corresponding type here.
"""

"""
    cull_bonds(bonds::Vector{Vector{Bond{D}}}, Js::Vector{Vector{S}}) where {D, S}

Culls presorted lists of bonds to keep only those with `bond.i` <= `bond.j`.
Bonds with `bond.i == bond.j` need to be further culled to only keep half,
  removing one from each pair with inverted `bond.n`.
"""
function cull_bonds(bonds::Vector{Vector{Bond{D}}}, Js::Vector{Vector{S}}) where {D, S}
    nb = length(bonds)

    culled_bonds = [Bond{D}[] for _ in 1:nb]
    culled_Js = [S[] for _ in 1:nb]
    for i in 1:nb
        for (bond, J) in zip(bonds[i], Js[i])
            if bond.i < bond.j
                push!(culled_bonds[i], bond)
                push!(culled_Js[i], J)
            end
            if bond.i == bond.j
                neg_bond = Bond{D}(bond.i, bond.j, -bond.n)
                if !(neg_bond in culled_bonds[i])
                    push!(culled_bonds[i], bond)
                    push!(culled_Js[i], J)
                end
            end
        end
    end
    (culled_bonds, culled_Js)
end


"""
    BondTable{D, T}

Stores all of the bonds of dimension D and associated data (e.g. interaction matrix) T
contained within a single pair interaction densely in memory, with quick access to bonds/Ts
on any sublattice.
"""
struct BondTable{D, T}
    bonds         :: Vector{Bond{D}}   # All bonds, sorted on first basis index
    culled_bonds  :: Vector{Bond{D}}   # Culled to minimal 1/2 bonds, sorted on first basis index
    data          :: Vector{T}         # Data T associated with each bond in `bonds`
    culled_data   :: Vector{T}         # Culled data T associated with each bond in `bonds`
    basis_indices :: Vector{Int}       # Indices to first bond of each i in `bonds`
    function BondTable(bonds::Vector{Vector{Bond{D}}}, data::Vector{Vector{T}}) where {D, T}
        (culled_bonds, culled_data) = cull_bonds(bonds, data)

        num_bonds_per_basis = map(length, bonds)
        basis_indices = cumsum(num_bonds_per_basis) .- length(bonds[1]) .+ 1
        push!(basis_indices, sum(num_bonds_per_basis) + 1)

        # Flatten everything from Vector{Vector{...}} to Vector{...}
        bonds = reduce(vcat, bonds)
        culled_bonds = reduce(vcat, culled_bonds)
        data = reduce(vcat, data)
        culled_data = reduce(vcat, culled_data)

        new{D, T}(bonds, culled_bonds, data, culled_data, basis_indices)
    end
end


function Base.show(io::IO, ::MIME"text/plain", bondtable::BondTable{D, T}) where {D, T}
    print(io, "$(length(bondtable))-element, $(length(bondtable.basis_indices)-1)-basis BondTable{$D, $T}")
end

Base.length(bondtable::BondTable) = length(bondtable.bonds)
## These functions generate iterators producing (bond, data) for various subsets of bonds ##
all_bonds(bondtable::BondTable) = zip(bondtable.bonds, bondtable.data)
culled_bonds(bondtable::BondTable) = zip(bondtable.culled_bonds, bondtable.culled_data)
function sublat_bonds(bondtable::BondTable, i::Int)
    @unpack bonds, data, basis_indices = bondtable
    @boundscheck checkindex(Bool, 1:length(basis_indices)-1, i) ? nothing : throw(BoundsError(bondtable, i))
    zip(bonds[basis_indices[i]:basis_indices[i+1]-1], data[basis_indices[i]:basis_indices[i+1]-1])
end

abstract type AbstractPairIntCPU{D} <: AbstractInteractionCPU end

"""
    HeisenbergCPU{D} <: AbstractPairIntCPU{D}

Implements an exchange interaction which is proportional to
the identity matrix.
"""
struct HeisenbergCPU{D} <: AbstractPairIntCPU{D}
    effJ         :: Float64                  # In Heisenberg interaction, all bonds have identical J
    bondtable    :: BondTable{D, Float64}    # Bonds store effective J's = S_i J S_j (all identical)
    label        :: String
end

"""
    DiagonalCouplingCPU{D} <: AbstractPairIntCPU{D}

Implements an exchange interaction where matrices on all bonds
are diagonal.
"""
struct DiagonalCouplingCPU{D} <: AbstractPairIntCPU{D}
    bondtable     :: BondTable{D, Vec3}       # Bonds store diagonal of effective J's = S_i J S_j
    label         :: String
end

"""
    GeneralCouplingCPU{D} <: AbstractPairIntCPU{D}

Implements the most generalized interaction, where matrices on
all bonds are full 3x3 matrices which vary bond-to-bond.
"""
struct GeneralCouplingCPU{D} <: AbstractPairIntCPU{D}
    bondtable    :: BondTable{D, Mat3}       # Bonds store effective J's = S_i J S_j
    label        :: String
end

# Helper functions producing predicates checking if a matrix is approximately
# Heisenberg or diagonal
isheisen(J, tol) = isapprox(diagm(fill(J[1,1], 3)), J; atol=tol)
isdiag(J, tol) = isapprox(diagm(diag(J)), J; atol=tol)
isheisen(tol) = Base.Fix2(isheisen, tol)
isdiag(tol) = Base.Fix2(isdiag, tol)

# Figures out the correct maximally-efficient backend type for a quadratic interaction
function convert_quadratic(int::QuadraticInteraction{D}, cryst::Crystal, sites_info::Vector{SiteInfo};
                           tol=1e-6) where {D}
    @unpack J, bond, label = int
    # Bonds and Js on each sublattice
    sorted_bonds = Vector{Vector{Bond{D}}}()
    sorted_Js = Vector{Vector{Mat3}}()
    for i in 1:nbasis(cryst)
        (bs, Js) = all_symmetry_related_couplings_for_atom(cryst, i, bond, J)
        push!(sorted_bonds, bs)
        push!(sorted_Js, Js)
    end

    # Product S_i S_j between sites connected by the bond -- same for all bonds in a class
    # To get interactions quadratic in the unit vectors stored by `SpinSystem`, we store
    #  (Si J Sj) internally.
    SiSj = sites_info[bond.i].S * sites_info[bond.j].S
    sorted_Js .*= SiSj

    if all(isheisen(tol), Base.Iterators.flatten(sorted_Js))
        scalar_sorted_Js = [
            [J[1,1] for J in Js]
            for Js in sorted_Js
        ]
        # In the Heisenberg case, all J values are identical -- store that scalar
        if !iszero(length(scalar_sorted_Js)) && !iszero(length(scalar_sorted_Js[1]))
            effJ = scalar_sorted_Js[1][1][1,1]
        else
            effJ = 0.0
        end
        bondtable = BondTable(sorted_bonds, scalar_sorted_Js)
        return HeisenbergCPU{D}(effJ, bondtable, label)
    elseif all(isdiag(tol), Base.Iterators.flatten(sorted_Js))
        vec_sorted_Js = [
            [diag(M) for M in Js]
            for Js in sorted_Js
        ]
        bondtable = BondTable(sorted_bonds, vec_sorted_Js)
        return DiagonalCouplingCPU{D}(bondtable, label)
    else
        bondtable = BondTable(sorted_bonds, sorted_Js)
        return GeneralCouplingCPU{D}(bondtable, label)
    end
end

function energy(spins::Array{Vec3}, heisenberg::HeisenbergCPU)
    @unpack effJ, bondtable = heisenberg
    E = 0.0

    latsize = size(spins)[2:end]
    for (bond, _) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = spins[i, cell]
            sⱼ = spins[j, offset(cell, n, latsize)]
            E += sᵢ ⋅ sⱼ
        end
    end
    return effJ * E
end

function energy(spins::Array{Vec3}, diag_coup::DiagonalCouplingCPU)
    bondtable = diag_coup.bondtable
    E = 0.0

    latsize = size(spins)[2:end]
    for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = spins[i, cell]
            sⱼ = spins[j, offset(cell, n, latsize)]
            E += (J .* sᵢ) ⋅ sⱼ
        end
    end
    return E
end

function energy(spins::Array{Vec3}, gen_coup::GeneralCouplingCPU)
    bondtable = gen_coup.bondtable
    E = 0.0

    latsize = size(spins)[2:end]
    for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            sᵢ = spins[i, cell]
            sⱼ = spins[j, offset(cell, n, latsize)]
            E += dot(sᵢ, J, sⱼ)
        end
    end
    return E
end

"Accumulates the local -∇ℋ coming from Heisenberg couplings into `B`"
@inline function _accum_neggrad!(B::Array{Vec3}, spins::Array{Vec3}, heisen::HeisenbergCPU)
    @unpack effJ, bondtable = heisen

    latsize = size(spins)[2:end]
    for (bond, _) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            offsetcell = offset(cell, n, latsize)
            B[i, cell] = B[i, cell] - effJ * spins[j, offsetcell]
            B[j, offsetcell] = B[j, offsetcell] - effJ * spins[i, cell]
        end
    end
end

"Accumulates the local -∇ℋ coming from diagonal couplings into `B`"
@inline function _accum_neggrad!(B::Array{Vec3}, spins::Array{Vec3}, diag_coup::DiagonalCouplingCPU)
    bondtable = diag_coup.bondtable

    latsize = size(spins)[2:end]
    for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            offsetcell = offset(cell, n, latsize)
            B[i, cell] = B[i, cell] - J .* spins[j, offsetcell]
            B[j, offsetcell] = B[j, offsetcell] - J .* spins[i, cell]
        end
    end
end

"Accumulates the local -∇ℋ coming from general couplings into `B`"
@inline function _accum_neggrad!(B::Array{Vec3}, spins::Array{Vec3}, gen_coup::GeneralCouplingCPU)
    bondtable = gen_coup.bondtable

    latsize = size(spins)[2:end]
    for (bond, J) in culled_bonds(bondtable)
        @unpack i, j, n = bond
        for cell in CartesianIndices(latsize)
            offsetcell = offset(cell, n, latsize)
            B[i, cell] = B[i, cell] - J * spins[j, offsetcell]
            B[j, offsetcell] = B[j, offsetcell] - J * spins[i, cell]
        end
    end
end
