abstract type Interaction end

struct ExternalField <: Interaction
    B :: Vec3
end

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

function cull_bonds(bonds::Vector{Vector{Bond{D}}}, Js::Vector{Vector{Mat3}}) where {D}
    nb = length(bonds)

    culled_bonds = [Bond{D}[] for _ in 1:nb]
    culled_Js = [Mat3[] for _ in 1:nb]
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

struct Heisenberg{D} <: Interaction
    J            :: Float64
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

"Create a Heisenberg interaction of strength `J` on all bonds symmetry-equivalent to `bond`"
function Heisenberg(J::Float64, cryst::Crystal, bond::Bond{D}, label::String="Heisen") where {D}
    sorted_bonds = [
        Symmetry.all_symmetry_related_bonds_for_atom(cryst, i, bond)
        for i in 1:nbasis(cryst)
    ]
    culled_bonds = cull_bonds(sorted_bonds)

    return Heisenberg{D}(J, sorted_bonds, culled_bonds, label)
end

"Equivalent to DiagonalCoupling with bonds = [(0,0,0)], but faster."
struct OnSite <: Interaction
    J     :: Vec3
    label :: String
end

function OnSite(J::Vec3)
    OnSite(J, "OnSite")
end

struct DiagonalCoupling{D} <: Interaction
    J            :: Vec3
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

function DiagonalCoupling(J::Vec3, cryst::Crystal, bond::Bond{D}, label::String="DiagJ") where {D}
    sorted_bonds = [
        Symmetry.all_symmetry_related_bonds_for_atom(cryst, i, bond)
        for i in 1:nbasis(cryst)
    ]
    culled_bonds = cull_bonds(sorted_bonds)

    DiagonalCoupling{D}(J, sorted_bonds, culled_bonds, label)
end

struct GeneralCoupling{D} <: Interaction
    Js           :: Vector{Vector{Mat3}}
    culled_Js    :: Vector{Vector{Mat3}}
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

function GeneralCoupling(J::Mat3, cryst::Crystal, bond::Bond{D}, label::String="GenJ") where {D}
    sorted_bonds = Vector{Vector{Bond{D}}}()
    sorted_Js = Vector{Vector{Mat3}}()
    for i in 1:nbasis(cryst)
        (bs, Js) = Symmetry.all_symmetry_related_interactions_for_atom(cryst, i, bond, J)
        push!(sorted_bonds, bs)
        push!(sorted_Js, Js)
    end
    (culled_bonds, culled_Js) = cull_bonds(sorted_bonds, sorted_Js)

    return GeneralCoupling{D}(sorted_Js, culled_Js, sorted_bonds, culled_bonds, label)
end

const PairInt{D} = Union{Heisenberg{D}, DiagonalCoupling{D}, GeneralCoupling{D}}

# Dipole-dipole interactions computed in real 3D space,
#   using a pre-computed interaction tensor.
struct DipoleReal <: Interaction
    int_mat :: OffsetArray{Mat3, 5, Array{Mat3, 5}}
end

# FFTW types for various relevant Fourier transform plans
const FTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const BFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const IFTPlan = AbstractFFTs.ScaledPlan{ComplexF64, BFTPlan, Float64}

# Dipole-dipole interactions computed in Fourier-space
struct DipoleFourier <: Interaction
    int_mat     :: Array{ComplexF64, 7}
    _spins_ft   :: Array{ComplexF64, 5}  # Space for Fourier-transforming spins
    _field_ft   :: Array{ComplexF64, 5}  # Space for holding Fourier-transformed fields
    _field_real :: Array{Float64, 5}     # Space for holding IFT-transformed fields
    _plan       :: FTPlan
    _ift_plan   :: IFTPlan
end

struct Hamiltonian{D}
    ext_field   :: Union{Nothing, ExternalField}
    heisenbergs :: Vector{Heisenberg{D}}
    on_sites    :: Vector{OnSite}
    diag_coups  :: Vector{DiagonalCoupling{D}}
    gen_coups   :: Vector{GeneralCoupling{D}}
    dipole_int  :: Union{Nothing, DipoleFourier}
end

function Hamiltonian{D}() where {D}
    return Hamiltonian{D}(nothing, [], [], [], nothing)
end

function Hamiltonian{D}(ints) where {D}
    ext_field   = nothing
    heisenbergs = Vector{Heisenberg{D}}()
    on_sites    = Vector{OnSite}()
    diag_coups  = Vector{DiagonalCoupling{D}}()
    gen_coups   = Vector{GeneralCoupling{D}}()
    dipole_int  = nothing
    for int in ints
        if isa(int, ExternalField)
            if !isnothing(ext_field)
                @warn "Provided multiple external fields. Only using last one."
            end
            ext_field = int
        elseif isa(int, Heisenberg)
            push!(heisenbergs, int)
        elseif isa(int, OnSite)
            push!(on_sites, int)
        elseif isa(int, DiagonalCoupling)
            push!(diag_coups, int)
        elseif isa(int, GeneralCoupling)
            push!(gen_coups, int)
        elseif isa(int, DipoleFourier)
            if !isnothing(dipole_int)
                @warn "Provided multiple dipole interactions. Only using last one."
            end
            dipole_int = int
        end
    end
    return Hamiltonian{D}(ext_field, heisenbergs, on_sites, diag_coups, gen_coups, dipole_int)
end

function Hamiltonian{D}(ints::Vararg{I}) where {D, I <: Interaction}
    Hamiltonian{D}(ints)
end