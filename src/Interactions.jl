"""Structs for defining various terms in a spin Hamiltonian.
"""

abstract type AbstractInteraction end      # Subtype this for user-facing interfaces
abstract type AbstractInteractionCPU end   # Subtype this for actual internal CPU implementations
abstract type AbstractInteractionGPU end   # Subtype this for actual internal GPU implementations
# abstract type AbstractAnisotropy <: AbstractInteraction end


struct QuadraticInteraction <: AbstractInteraction
    J     :: Mat3
    bond  :: Bond
    label :: String
end

function Base.show(io::IO, ::MIME"text/plain", int::QuadraticInteraction)
    b = repr("text/plain", int.bond)
    J = int.J
    if J ≈ -J'                             # Catch purely DM interactions
        x = J[2, 3]
        y = J[3, 1]
        z = J[1, 2]
        @printf io "dm_interaction([%.4g, %.4g, %.4g], %s)" x y z b
    elseif diagm(fill(J[1,1], 3)) ≈ J      # Catch Heisenberg interactions
        @printf io "heisenberg(%.4g, %s)" J[1,1] b
    elseif diagm(diag(J)) ≈ J              # Catch diagonal interactions
        @printf io "exchange(diagm([%.4g, %.4g, %.4g]), %s)" J[1,1] J[2,2] J[3,3] b
    else                                   # Rest -- general exchange interactions
        @printf io "exchange([%.4g %.4g %.4g; %.4g %.4g %.4g; %.4g %.4g %.4g], %s)" J[1,1] J[1,2] J[1,3] J[2,1] J[2,2] J[2,3] J[3,1] J[3,2] J[3,3] b
        # TODO: Figure out how to reenable this depending on context:
        # @printf io "exchange([%.4f %.4f %.4f\n"   J[1,1] J[1,2] J[1,3]
        # @printf io "          %.4f %.4f %.4f\n"   J[2,1] J[2,2] J[2,3]
        # @printf io "          %.4f %.4f %.4f],\n" J[3,1] J[3,2] J[3,3]
        # @printf io "    %s)" b
    end
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
    J = SA[      0.0  DMvec[3] -DMvec[2]
           -DMvec[3]       0.0  DMvec[1]
            DMvec[2] -DMvec[1]      0.0]
    QuadraticInteraction(J, bond, label)
end

struct OperatorAnisotropy <: AbstractInteraction
    op    :: DP.AbstractPolynomialLike
    site  :: Int
    label :: String # Maybe remove
end


"""
    anisotropy(op, site)

Creates a general anisotropy specified as a polynomial of spin operators `𝒮` or
Stevens operators `𝒪`.
"""
function anisotropy(op::DP.AbstractPolynomialLike, site, label="OperatorAniso")
    OperatorAnisotropy(op, site, label)
end


"""
    quadratic_anisotropy(J, site, label="Anisotropy")

Creates a quadratic single-ion anisotropy,
```math
    ∑_i 𝐒_i^T J^{(i)} 𝐒_i
```
where ``i`` runs over all sublattices that are symmetry equivalent to `site`,
and ``J^{(i)}`` is the covariant transformation of the ``3 × 3`` anisotropy
matrix `J` appropriate for ``i``. Without loss of generality, we require that
`J` is symmetric.
"""
function quadratic_anisotropy(J, site::Int, label::String="Anisotropy")
    if !(J ≈ J')
        error("Single-ion anisotropy must be symmetric.")
    end
    OperatorAnisotropy(𝒮'*Mat3(J)*𝒮, site, label)
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
    OperatorAnisotropy(-D*(𝒮⋅n)^2, site, label)
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
    OperatorAnisotropy(+D*(𝒮⋅n)^2, site, label)
end

struct ExternalField <: AbstractInteraction
    B :: Vec3
end

"""
    external_field(B::Vec3)

Adds an external field ``𝐁`` with Zeeman coupling,

```math
    -∑_i 𝐁 ⋅ 𝐌_i.
```

The magnetic moments are ``𝐌_i = μ_B g 𝐒_i`` where ``g`` is the g-factor or
g-tensor, and the spin magnitude ``|𝐒_i|`` is typically a multiple of 1/2. The
Bohr magneton ``μ_B`` is a physical constant, with numerical value determined by
the unit system.
"""
external_field(B) = ExternalField(Vec3(B))

function Base.show(io::IO, ::MIME"text/plain", int::ExternalField)
    B = int.B
    @printf io "external_field([%.4g, %.4g, %.4g])" B[1] B[2] B[3]
end


#= Energy and field functions for "simple" interactions that aren't geometry-dependent.
   See Hamiltonian.jl for expectations on `_accum_neggrad!` functions.
=#

struct ExternalFieldCPU
    effBs :: Vector{Vec3}  # |S_b|gᵀB for each basis index b
end

function ExternalFieldCPU(ext_field::ExternalField, site_infos::Vector{SiteInfo}; μB)
    # As E = -∑_i 𝐁^T g 𝐒_i, we can precompute effB = g^T S B, so that
    #  we can compute E = -∑_i effB ⋅ 𝐬_i during simulation.
    # However, S_i may be basis-dependent, so we need to store an effB
    #  per sublattice.
    effBs = [μB * site.g' * ext_field.B for site in site_infos]
    ExternalFieldCPU(effBs)
end

function energy(dipoles::Array{Vec3, 4}, field::ExternalFieldCPU)
    E = 0.0
    @inbounds for site in 1:size(dipoles)[end]
        effB = field.effBs[site]
        for s in selectdim(dipoles, 4, site)
            E += effB ⋅ s
        end
    end
    return -E
end

"Accumulates the negative local Hamiltonian gradient coming from the external field"
@inline function _accum_neggrad!(B::Array{Vec3, 4}, field::ExternalFieldCPU)
    @inbounds for site in 1:size(B)[end]
        effB = field.effBs[site]
        for cell in CartesianIndices(size(B)[1:3])
            B[cell, site] = B[cell, site] + effB
        end
    end
end
