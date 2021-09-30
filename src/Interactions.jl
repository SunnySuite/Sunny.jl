"""Structs for defining various terms in a spin Hamiltonian.
"""

abstract type Interaction end      # Subtype this for user-facing interfaces
abstract type InteractionCPU end   # Subtype this for actual internal CPU implementations
abstract type InteractionGPU end   # Subtype this for actual internal GPU implementations

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

"""
    Heisenberg(J, bond::Bond{D}, label::String="Heisen")

Defines a Heisenberg interaction on a `D`-dimensional lattice. Specifically, the term

```math
    J ∑_{⟨ij⟩} 𝐒_i ⋅ 𝐒_j
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once.
"""
struct Heisenberg{D} <: Interaction
    J     :: Float64
    bond  :: Bond{D}
    label :: String
end

Heisenberg(J, bond::Bond{D}, label::String="Heisen") where {D} = Heisenberg{D}(J, bond, label)

"""
    OnSite(J, label="OnSite")

Construct an on-site anisotropy term specified by the vector J. Specifically, the term

```math
    ∑_i ∑_α J_α S_{i, α}^2
```
Equivalent to `DiagonalCoupling` with `bonds = [(0,0,0)]`, but faster.
"""
struct OnSite <: Interaction
    J     :: Vec3
    label :: String
    function OnSite(J, label::String="OnSite")
        new(J, label)
    end
end

"""
    DiagonalCoupling(J, bond::Bond{D}, label::String="DiagCoup")

Defines a diagonal exchange interaction on a `D`-dimensional lattice. Specifically, the term

```math
    ∑_{⟨ij⟩} ∑_α J_α S_{i,α} ⋅ S_{j,α} 
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once. ``J`` is a
length-3 vector which may vary from bond to bond, under symmetry constraints.
"""
struct DiagonalCoupling{D} <: Interaction
    J     :: Vec3
    bond  :: Bond{D}
    label :: String
end

DiagonalCoupling(J, bond::Bond{D}, label::String="DiagCoup") where {D} = DiagonalCoupling{D}(J, bond, label)

"""
    GeneralCoupling(J, bond::Bond{D}, label::String="GenCoup")

Defines a general exchange interaction on a `D`-dimensional lattice. Specifically, the term

```math
    ∑_{⟨ij⟩} 𝐒_i^⊤ J 𝐒_j
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once. ``J`` is a
``3 × 3`` matrix which may vary from bond to bond, under symmetry constraints.
"""
struct GeneralCoupling{D} <: Interaction
    J     :: Mat3
    bond  :: Bond{D}
    label :: String
end

GeneralCoupling(J, bond::Bond{D}, label::String="GenCoup") where {D} = GeneralCoupling{D}(J, bond, label)

const PairInt{D} = Union{Heisenberg{D}, DiagonalCoupling{D}, GeneralCoupling{D}}

"""
    DipoleDipole(strength, extent::Int=4, η::Float64=0.5)

Defines long-range dipole-dipole interactions under the Ewald summation convention,
assumed to be on a 3-dimensional lattice. Specifically, the term

```math
    ∑_{⟨ij⟩}
```
evaluated under the Ewald summation convention.

`extent` controls the number of periodic copies of the unit cell summed over in the
Ewald summation (higher is more accurate, but higher creation-time cost), while `η`
controls the direct/reciprocal-space tradeoff in the Ewald summation.
"""
struct DipoleDipole <: Interaction
    strength :: Float64
    extent   :: Int
    η        :: Float64
    function DipoleDipole(strength; extent=4, η=0.5)
        new(strength, extent, η)
    end
end

#= Energy and field functions for "simple" interactions that aren't geometry-dependent =#

function energy(spins::Array{Vec3}, field::ExternalField)
    B = field.B
    E = 0.0
    for S in spins
        E += S ⋅ B
    end
    return -E
end

function energy(spins::Array{Vec3}, on_site::OnSite)
    J = on_site.J
    E = 0.0
    for S in spins
        E += S ⋅ (J .* S)
    end
    return E
end

"Accumulates the local field coming from the external field"
@inline function _accum_field!(B::Array{Vec3}, field::ExternalField)
    for idx in eachindex(B)
        B[idx] = B[idx] + field.B
    end
end

"Accumulates the local field coming from on-site terms"
@inline function _accum_field!(B::Array{Vec3}, spins::Array{Vec3}, on_site::OnSite)
    J = on_site.J
    for idx in eachindex(spins)
        S = spins[idx]
        B[idx] = B[idx] - 2 * J .* S
    end
end