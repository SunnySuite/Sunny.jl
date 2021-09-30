"""Structs and functions for implementing various pair interaction energies and fields,
    given a specific lattice to operate on.
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
function cull_bonds(bonds::Vector{Vector{Bond{D}}}) :: Vector{Vector{Bond{D}}} where {D}
    nb = length(bonds)

    culled_bonds = [Bond{D}[] for _ in 1:nb]
    for i in 1:nb
        for bond in bonds[i]
            if bond.i < bond.j
                push!(culled_bonds[i], bond)
            end
            if bond.i == bond.j
                neg_bond = Bond{D}(bond.i, bond.j, -bond.n)
                if !(neg_bond in culled_bonds[i])
                    push!(culled_bonds[i], bond)
                end
            end
        end
    end

    culled_bonds
end

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

The type implementing `Heisenberg` on CPU.
"""
struct HeisenbergCPU{D} <: InteractionCPU
    J            :: Float64
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

"""
    HeisenbergCPU(heisen::Heisenberg{D}, cryst::Crystal)

Construct a `Heisenberg{D}` interaction of strength `J` acting on all bonds
symmetry-equivalent to `bond` in `cryst`.
"""
function HeisenbergCPU(heisen::Heisenberg{D}, cryst::Crystal) where {D}
    @unpack J, bond, label = heisen
    sorted_bonds = [
        all_symmetry_related_bonds_for_atom(cryst, i, bond)
        for i in 1:nbasis(cryst)
    ]
    culled_bonds = cull_bonds(sorted_bonds)

    return HeisenbergCPU{D}(J, sorted_bonds, culled_bonds, label)
end

"""
    DiagonalCouplingCPU{D}

The type implementing `DiagonalCoupling` on CPU.
"""
struct DiagonalCouplingCPU{D} <: InteractionCPU
    Js            :: Vector{Vector{Vec3}}
    culled_Js     :: Vector{Vector{Vec3}}
    bonds         :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds  :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label         :: String
end

function DiagonalCouplingCPU(diag_coup::DiagonalCoupling{D}, cryst::Crystal; tol=1e-8) where {D}
    @unpack J, bond, label = diag_coup
    sorted_bonds = Vector{Vector{Bond{D}}}()
    sorted_Js = Vector{Vector{Vec3}}()
    for i in 1:nbasis(cryst)
        (bs, Js) = all_symmetry_related_interactions_for_atom(cryst, i, bond, diagm(J))
        # Check that all J matrices are actually diagonal
        diagJs = Vector{Vec3}()
        for J in Js
            if max(abs(J[1,2]), abs(J[1,3]), abs(J[2,3])) > tol
                error("A symmetry-transformed J in interaction '$(label)' is non-diagonal. Please use GeneralCoupling.")
            end
            push!(diagJs, diag(J))
        end
        push!(sorted_bonds, bs)
        push!(sorted_Js, diagJs)
    end
    (culled_bonds, culled_Js) = cull_bonds(sorted_bonds, sorted_Js)

    DiagonalCouplingCPU{D}(sorted_Js, culled_Js, sorted_bonds, culled_bonds, label)
end

"""
    GeneralCouplingCPU{D}

The type implementing `GeneralCoupling` on CPU.
"""
struct GeneralCouplingCPU{D} <: InteractionCPU
    Js           :: Vector{Vector{Mat3}}
    culled_Js    :: Vector{Vector{Mat3}}
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

function GeneralCouplingCPU(gen_coup::GeneralCoupling{D}, cryst::Crystal) where {D}
    @unpack J, bond, label = gen_coup
    sorted_bonds = Vector{Vector{Bond{D}}}()
    sorted_Js = Vector{Vector{Mat3}}()
    for i in 1:nbasis(cryst)
        (bs, Js) = all_symmetry_related_interactions_for_atom(cryst, i, bond, J)
        push!(sorted_bonds, bs)
        push!(sorted_Js, Js)
    end
    (culled_bonds, culled_Js) = cull_bonds(sorted_bonds, sorted_Js)

    return GeneralCouplingCPU{D}(sorted_Js, culled_Js, sorted_bonds, culled_bonds, label)
end

const PairIntCPU{D} = Union{HeisenbergCPU{D}, DiagonalCouplingCPU{D}, GeneralCouplingCPU{D}}

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
