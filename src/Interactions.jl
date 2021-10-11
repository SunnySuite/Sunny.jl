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
    QuadraticInteraction{D}

Defines a general quadratic interaction. Specifically, the term

```math
    ∑_{⟨ij⟩} 𝐒_i^⊤ J_{ij} 𝐒_j
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once. ``J`` is a
``3 × 3`` matrix which may vary from bond to bond, under symmetry constraints.
"""
struct QuadraticInteraction{D} <: Interaction
    J     :: Mat3
    bond  :: Bond{D}
    label :: String
    function QuadraticInteraction(J, bond::Bond{D}, label::String) where {D}
        if bond.i == bond.j && all(isequal(0), bond.n)
            error("This interaction looks on-site. Please use `OnSiteQuadratic`.")
        end
        new{D}(J, bond, label)
    end
end

"""
    OnSiteQuadratic

Defines a general on-site quadratic anisotropy. Specifically, the term

```math
    ∑_i ∑_α J_α S_{(b, i), α}^2
```

for a given vector ``𝐉`` and sublattice index ``b``.
"""
struct OnSiteQuadratic <: Interaction
    J     :: Vec3
    site  :: Int
    label :: String
end

"""
    heisenberg(J, bond::Bond, label::String="Heisen")

Creates a Heisenberg interaction term
```math
    J ∑_{⟨ij⟩} 𝐒_i ⋅ 𝐒_j
```
where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once.
"""
heisenberg(J, bond::Bond, label::String="Heisen") = QuadraticInteraction(diagm(fill(J, 3)), bond, label)


"""
    dm_interaction(DMvec, bond::Bond, label::String="DMInt")

Creates a DM Interaction term
```math
    ∑_{⟨ij⟩} 𝐃_{ij} ⋅ (𝐒_i × 𝐒_j)
```
where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once.
"""
function dm_interaction(DMvec, bond::Bond, label::String="DMInt")
    J = SA[     0.0   DMvec[3] -DMvec[2];
           -DMvec[3]       0.0  DMvec[1];
            DMvec[2] -DMvec[1]      0.0]
    QuadraticInteraction(J, bond, label)
end


"""
    onsite_anisotropy(J, site, label="OnSiteAniso")

Creates an on-site anisotropy term specified by the vector `J`,
applied to all sites of `site` type in the crystal.
Specifically, the term
```math
    ∑_i ∑_α J_α S_{(b, i), α}^2
```
where ``b`` is the `site` index.
"""
function onsite_anisotropy(J, site::Int, label::String="OnSiteAniso")
    OnSiteQuadratic(J, site, label)
end

"""
    exchange(J, bond::Bond, label="Exchange")

Defines a general quadratic interaction. Specifically, the term

```math
    ∑_{⟨ij⟩} 𝐒_i^⊤ J_{ij} 𝐒_j
```

where ``⟨ij⟩`` indicates a sum over all bonds in the lattice of a certain
symmetry equivalence class, with each bond appearing only once. ``J`` is a
``3 × 3`` matrix which may vary from bond to bond, under symmetry constraints.
`i` and `j` can be the same index, in which case this is an on-site
anisotropy.
"""
function exchange(J, bond::Bond, label::String="Exchange")
    QuadraticInteraction(J, bond, label)
end

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

"Accumulates the local field coming from the external field"
@inline function _accum_field!(B::Array{Vec3}, field::ExternalField)
    for idx in eachindex(B)
        B[idx] = B[idx] + field.B
    end
end
