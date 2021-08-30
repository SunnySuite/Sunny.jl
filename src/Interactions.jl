abstract type Interaction end

"""
    ExternalField(B::Vec3)

Defines an external field acting on each spin, specifically the term

```math
    -∑_i 𝐁 ⋅ 𝐒_i
```
"""
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

"""
    Heisenberg{D}

Defines a Heisenberg interaction on a `D`-dimensional lattice. Specifically, the term

```math
    J ∑_{⟨ij⟩} 𝐒_i ⋅ 𝐒_j 
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once.
"""
struct Heisenberg{D} <: Interaction
    J            :: Float64
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

"""
    Heisenberg(J::Float64, cryst::Crystal, bond::Bond{D}, label="Heisenberg")

Construct a `Heisenberg{D}` interaction of strength `J` acting on all bonds
symmetry-equivalent to `bond` in `cryst`.
"""
function Heisenberg(J::Float64, cryst::Crystal, bond::Bond{D}, label::String="Heisen") where {D}
    sorted_bonds = [
        Symmetry.all_symmetry_related_bonds_for_atom(cryst, i, bond)
        for i in 1:nbasis(cryst)
    ]
    culled_bonds = cull_bonds(sorted_bonds)

    return Heisenberg{D}(J, sorted_bonds, culled_bonds, label)
end

"""
    OnSite(J::Vec3, label="OnSite")

Construct an on-site anisotropy term specified by the vector J. Specifically, the term

```math
    ∑_i ∑_α J_α S_{i, α}^2
```
Equivalent to `DiagonalCoupling` with `bonds = [(0,0,0)]`, but faster.
"""
struct OnSite <: Interaction
    J     :: Vec3
    label :: String
end

function OnSite(J::Vec3)
    OnSite(J, "OnSite")
end

"""
    DiagonalCoupling{D}

Defines a diagonal exchange interaction on a `D`-dimensional lattice. Specifically, the term

```math
    ∑_{⟨ij⟩} ∑_α J_α S_{i,α} ⋅ S_{j,α} 
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once.

**Warning**: This type currently (and incorrectly) assumes that ``𝐉`` takes
the same form on every bond within the equivalency class.
"""
struct DiagonalCoupling{D} <: Interaction
    J            :: Vec3
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

"""
    DiagonalCoupling(J::Vec3, cryst::Crystal, bond::Bond{D}, label="Heisenberg")

Construct a `DiagonalCoupling{D}` interaction with diagonal elements specified
by `J`, acting on all bonds symmetry-equivalent to `bond` in `cryst`.
"""
function DiagonalCoupling(J::Vec3, cryst::Crystal, bond::Bond{D}, label::String="DiagJ") where {D}
    sorted_bonds = [
        Symmetry.all_symmetry_related_bonds_for_atom(cryst, i, bond)
        for i in 1:nbasis(cryst)
    ]
    culled_bonds = cull_bonds(sorted_bonds)

    DiagonalCoupling{D}(J, sorted_bonds, culled_bonds, label)
end

"""
    GeneralCoupling{D}

Defines a general exchange interaction on a `D`-dimensional lattice. Specifically, the term

```math
    ∑_{⟨ij⟩} 𝐒_i^⊤ J 𝐒_j
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once. ``J`` is a
``3 × 3`` matrix which may vary from bond to bond, under symmetry constraints.
"""
struct GeneralCoupling{D} <: Interaction
    Js           :: Vector{Vector{Mat3}}
    culled_Js    :: Vector{Vector{Mat3}}
    bonds        :: Vector{Vector{Bond{D}}}  # Each outer Vector is bonds on one sublattice
    culled_bonds :: Vector{Vector{Bond{D}}}  # Like `bonds`, but culled to the minimal 1/2 bonds
    label        :: String
end

"""
    GeneralCoupling(J::Mat3, cryst::Crystal, bond::Bond{D}, label="GenJ")

Construct a `GeneralCoupling{D}` interaction specified by the `J` matrix
 on the given `bond` in `cryst`. Automatically propagates the correctly
 transformed exchange matrices to all symmetry-equivalent bonds.
"""
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

"""
Dipole-dipole interactions computed in real-space. `DipoleFourier` should
be preferred in actual simulations, but this type persists as a cross-check
to test the Fourier-space calculations.
"""
struct DipoleReal <: Interaction
    int_mat :: OffsetArray{Mat3, 5, Array{Mat3, 5}}
end

# FFTW types for various relevant Fourier transform plans
const FTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const BFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const IFTPlan = AbstractFFTs.ScaledPlan{ComplexF64, BFTPlan, Float64}

"""
Dipole-dipole interactions computed in Fourier-space. Should produce
identical results (up to numerical precision) as `DipoleReal`, but
is asymptotically faster.
"""
struct DipoleFourier <: Interaction
    int_mat     :: Array{ComplexF64, 7}
    _spins_ft   :: Array{ComplexF64, 5}  # Space for Fourier-transforming spins
    _field_ft   :: Array{ComplexF64, 5}  # Space for holding Fourier-transformed fields
    _field_real :: Array{Float64, 5}     # Space for holding IFT-transformed fields
    _plan       :: FTPlan
    _ift_plan   :: IFTPlan
end

"""
    Hamiltonian{D}

Defines a Hamiltonian for a `D`-dimensional spin system.
"""
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

"""
    Hamiltonian{D}(ints)
    Hamiltonian{D}(ints...)

Construct a `Hamiltonian{D}` from a list of `Interaction`'s.
"""
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