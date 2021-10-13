"""
Structs and functions for implementing various pair interaction energies
 and fields, given a specific lattice to operate on.
Upon creation of a SpinSystem, all pair interactions get converted into their
 corresponding type here.
"""

"Presorts a flat list of bonds into a nested list, with the outer list corresponding to `bond.i`"
function presort_bonds(bonds::Vector{Bond{D}}) :: Vector{Vector{Bond{D}}} where {D}
    nb = maximum(b->b.i, bonds)

    sorted_bonds = [Bond{D}[] for _ in 1:nb]
    for bond in bonds
        push!(sorted_bonds[bond.i], bond)
    end

    sorted_bonds
end

"""Culls presorted lists of bonds to keep only those with `bond.i` <= `bond.j`.
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
    HeisenbergCPU{D}

Implements an exchange interaction which is proportional to
the identity matrix.
"""
struct HeisenbergCPU{D} <: InteractionCPU
    J            :: Float64
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

"""
    DiagonalCouplingCPU{D}

Implements an exchange interaction where matrices on all bonds
are diagonal.
"""
struct DiagonalCouplingCPU{D} <: InteractionCPU
    Js            :: Vector{Vector{Vec3}}
    culled_Js     :: Vector{Vector{Vec3}}
    bonds         :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds  :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label         :: String
end

"""
    GeneralCouplingCPU{D}

Implements the most generalized interaction, where matrices on
all bonds are full 3x3 matrices which vary bond-to-bond.
"""
struct GeneralCouplingCPU{D} <: InteractionCPU
    Js           :: Vector{Vector{Mat3}}
    culled_Js    :: Vector{Vector{Mat3}}
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

const PairInt{D} = Union{HeisenbergCPU{D}, DiagonalCouplingCPU{D}, GeneralCouplingCPU{D}}

# Helper functions producing predicates checking if a matrix is approximately
# Heisenberg or diagonal
isheisen(J, tol) = isapprox(diagm(fill(J[1,1], 3)), J; atol=tol)
isdiag(J, tol) = isapprox(diagm(diag(J)), J; atol=tol)
isheisen(tol) = Base.Fix2(isheisen, tol)
isdiag(tol) = Base.Fix2(isdiag, tol)

# Figures out the correct maximally-efficient backend type for a quadratic interaction
function convert_quadratic(int::QuadraticInteraction{D}, cryst; tol=1e-6) where {D}
    @unpack J, bond, label = int
    sorted_bonds = Vector{Vector{Bond{D}}}()
    sorted_Js = Vector{Vector{Mat3}}()
    for i in 1:nbasis(cryst)
        (bs, Js) = all_symmetry_related_couplings_for_atom(cryst, i, bond, J)
        push!(sorted_bonds, bs)
        push!(sorted_Js, Js)
    end
    (culled_bonds, culled_Js) = cull_bonds(sorted_bonds, sorted_Js)

    if all(isheisen(tol), Base.Iterators.flatten(culled_Js))
        return HeisenbergCPU{D}(J[1,1], sorted_bonds, culled_bonds, label)
    elseif all(isdiag(tol), Base.Iterators.flatten(culled_Js))
        vec_sorted_Js = [
            [diag(M) for M in Js]
            for Js in sorted_Js
        ]
        vec_culled_Js = [
            [diag(M) for M in Js]
            for Js in culled_Js
        ]
        return DiagonalCouplingCPU{D}(
            vec_sorted_Js, vec_culled_Js,
            sorted_bonds, culled_bonds, label
        )
    else
        return GeneralCouplingCPU{D}(
            sorted_Js, culled_Js,
            sorted_bonds, culled_bonds, label
        )
    end
end

function energy(spins::Array{Vec3}, heisenberg::HeisenbergCPU)
    @unpack J, culled_bonds = heisenberg
    E = 0.0
    latsize = size(spins)[2:end]
    for (i, bonds) in enumerate(culled_bonds)
        for bond in bonds
            @unpack j, n = bond
            for cell in CartesianIndices(latsize)
                Sᵢ = spins[i, cell]
                Sⱼ = spins[j, offset(cell, n, latsize)]
                E += Sᵢ ⋅ Sⱼ
            end
        end
    end
    return J * E
end

function energy(spins::Array{Vec3}, diag_coup::DiagonalCouplingCPU)
    @unpack culled_Js, culled_bonds = diag_coup
    E = 0.0
    latsize = size(spins)[2:end]
    for (i, (Js, bonds)) in enumerate(zip(culled_Js, culled_bonds))
        for (J, bond) in zip(Js, bonds)
            @unpack j, n = bond
            for cell in CartesianIndices(latsize)
                Sᵢ = spins[i, cell]
                Sⱼ = spins[j, offset(cell, n, latsize)]
                E += (J .* Sᵢ) ⋅ Sⱼ
            end
        end
    end
    return E
end

function energy(spins::Array{Vec3}, gen_coup::GeneralCouplingCPU)
    @unpack culled_Js, culled_bonds = gen_coup
    E = 0.0
    latsize = size(spins)[2:end]
    for (i, (Js, bonds)) in enumerate(zip(culled_Js, culled_bonds))
        for (J, bond) in zip(Js, bonds)
            @unpack j, n = bond
            for cell in CartesianIndices(latsize)
                Sᵢ = spins[i, cell]
                Sⱼ = spins[j, offset(cell, n, latsize)]
                E += dot(Sᵢ, J, Sⱼ)
            end
        end
    end
    return E
end

"Accumulates the local field coming from Heisenberg couplings"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, heisen::HeisenbergCPU)
    latsize = size(spins)[2:end]
    @unpack J, culled_bonds = heisen
    for (i, bonds) in enumerate(culled_bonds)
        for bond in bonds
            @unpack j, n = bond
            for cell in CartesianIndices(latsize)
                offsetcell = offset(cell, n, latsize)
                B[i, cell] = B[i, cell] - J * spins[j, offsetcell]
                B[j, offsetcell] = B[j, offsetcell] - J * spins[i, cell]
            end
        end
    end
end

"Accumulates the local field coming from diagonal couplings"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, diag_coup::DiagonalCouplingCPU)
    latsize = size(spins)[2:end]

    @unpack culled_Js, culled_bonds = diag_coup
    for (i, (Js, bonds)) in enumerate(zip(culled_Js, culled_bonds))
        for (J, bond) in zip(Js, bonds)
            @unpack j, n = bond
            for cell in CartesianIndices(latsize)
                offsetcell = offset(cell, n, latsize)
                B[i, cell] = B[i, cell] - J .* spins[j, offsetcell]
                B[j, offsetcell] = B[j, offsetcell] - J .* spins[i, cell]
            end
        end
    end
end

"Accumulates the local field coming from general couplings"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, gen_coup::GeneralCouplingCPU)
    latsize = size(spins)[2:end]

    @unpack culled_Js, culled_bonds = gen_coup
    for (i, (Js, bonds)) in enumerate(zip(culled_Js, culled_bonds))
        for (J, bond) in zip(Js, bonds)
            @unpack j, n = bond
            for cell in CartesianIndices(latsize)
                offsetcell = offset(cell, n, latsize)
                B[i, cell] = B[i, cell] - J * spins[j, offsetcell]
                B[j, offsetcell] = B[j, offsetcell] - J * spins[i, cell]
            end
        end
    end
end
