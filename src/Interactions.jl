"""Structs for defining various terms in a spin Hamiltonian.
"""

abstract type Interaction end      # Subtype this for user-facing interfaces
abstract type InteractionCPU end   # Subtype this for actual internal CPU implementations
abstract type InteractionGPU end   # Subtype this for actual internal GPU implementations


struct QuadraticInteraction{D} <: Interaction
    J     :: Mat3
    bond  :: Bond{D}
    label :: String
end

struct OnSiteQuadratic <: Interaction
    J     :: Mat3
    site  :: Int
    label :: String
end

struct ExternalField <: Interaction
    B :: Vec3
end

struct DipoleDipole <: Interaction
    strength :: Float64
    extent   :: Int
    η        :: Float64
end


"""
    exchange(J, bond::Bond, label="Exchange")

Creates a quadratic interaction,

```math
    ∑_{⟨ij⟩} 𝐒_i^T J^{(ij)} 𝐒_j
```

where ``⟨ij⟩`` runs over all bonds (not doubly counted) that are symmetry
equivalent to `bond`. The ``3 × 3`` interaction matrix ``J^{(ij)}`` is the
covariant transformation of `J` appropriate for the bond ``⟨ij⟩``.
"""
function exchange(J, bond::Bond, label::String="Exchange")
    QuadraticInteraction(Mat3(J), bond, label)
end


"""
    heisenberg(J, bond::Bond, label::String="Heisen")

Creates a Heisenberg interaction
```math
    J ∑_{⟨ij⟩} 𝐒_i ⋅ 𝐒_j
```
where ``⟨ij⟩`` runs over all bonds symmetry equivalent to `bond`.
"""
heisenberg(J, bond::Bond, label::String="Heisen") = QuadraticInteraction(J*Mat3(I), bond, label)


"""
    dm_interaction(DMvec, bond::Bond, label::String="DMInt")

Creates a DM Interaction
```math
    ∑_{⟨ij⟩} 𝐃^{(ij)} ⋅ (𝐒_i × 𝐒_j)
```
where ``⟨ij⟩`` runs over all bonds symmetry equivalent to `bond`, and
``𝐃^{(ij)}`` is the covariant transformation of the DM pseudo-vector `DMvec`
appropriate for the bond ``⟨ij⟩``.
"""
function dm_interaction(DMvec, bond::Bond, label::String="DMInt")
    J = SA[     0.0   DMvec[3] -DMvec[2]
           -DMvec[3]       0.0  DMvec[1]
            DMvec[2] -DMvec[1]      0.0]
    QuadraticInteraction(J, bond, label)
end


"""
    single_ion_anisotropy(J, site, label="Anisotropy")

Creates a quadratic single-ion anisotropy,
```math
    ∑_i 𝐒_i^T J^{(i)} 𝐒_i
```
where ``i`` runs over all sublattices that are symmetry equivalent to `site`,
and ``J^{(i)}`` is the covariant transformation of the ``3 × 3`` anisotropy
matrix `J` appropriate for ``i``. Without loss of generality, we require that
`J` is symmetric.
"""
function single_ion_anisotropy(J, site::Int, label::String="Anisotropy")
    OnSiteQuadratic(Mat3(J), site, label)
end


"""
    easy_axis(D, n, site, label="EasyAxis")

Creates an easy axis anisotropy,
```math
    - D ∑_i (𝐧̂^{(i)}⋅𝐒_i)^2
```
where ``i`` runs over all sublattices that are symmetry equivalent to `site`,
``𝐧̂^{(i)}`` is the covariant transformation of the unit vector `n`, and ``D > 0``
is the interaction strength.
"""
function easy_axis(D, n, site::Int, label::String="EasyAxis")
    if D <= 0
        error("Parameter `D` must be nonnegative.")
    end
    if !(norm(n) ≈ 1)
        error("Parameter `n` must be a unit vector. Consider using `normalize(n)`.")
    end
    OnSiteQuadratic(-D*Mat3(n*n'), site, label)
end


"""
    easy_plane(D, n, site, label="EasyPlane")

Creates an easy plane anisotropy,
```math
    + D ∑_i (𝐧̂^{(i)}⋅𝐒_i)^2
```
where ``i`` runs over all sublattices that are symmetry equivalent to `site`,
``𝐧̂^{(i)}`` is the covariant transformation of the unit vector `n`, and ``D > 0``
is the interaction strength.
"""
function easy_plane(D, n, site::Int, label::String="EasyAxis")
    if D <= 0
        error("Parameter `D` must be nonnegative.")
    end
    if !(norm(n) ≈ 1)
        error("Parameter `n` must be a unit vector. Consider using `normalize(n)`.")
    end
    OnSiteQuadratic(+D*Mat3(n*n'), site, label)
end


"""
    external_field(B::Vec3)

Adds an external field ``𝐁`` and the energy term

```math
    -∑_i 𝐁 ⋅ 𝐦_i.
```

The magnetic moments are ``𝐦_i = g 𝐬_i`` where ``g`` is in general a tensor and
the spin magnitude ``|𝐬_i|`` is typically a multiple of 1/2.
"""
external_field(B) = ExternalField(Vec3(B))


"""
    dipole_dipole(; extent::Int=4, η::Float64=0.5)

Adds long-range dipole-dipole interactions,

```math
    -(μ₀/4π) ∑_{ij}  (3 (𝐦_j⋅𝐫̂_{ij})(𝐦_i⋅𝐫̂_{ij}) - 𝐦_i⋅𝐦_j) / |𝐫_{ij}|^3
```

where the sum is over all pairs ``i \neq j``, singly counted, including
periodic images, regularized using the Ewald summation convention. The magnetic
moments are ``𝐦_i = g 𝐬_i`` where ``g`` is in general a tensor and the spin
magnitude ``|𝐬_i|`` is typically a multiple of 1/2.

A three-dimensional system is required.

`extent` controls the number of periodic copies of the unit cell summed over in
the Ewald summation (higher is more accurate, but higher creation-time cost),
while `η` controls the direct/reciprocal-space tradeoff in the Ewald summation.
"""
dipole_dipole(; extent=4, η=0.5) = DipoleDipole(1.0, extent, η)


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
