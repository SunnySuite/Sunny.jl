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

struct BiQuadraticInteraction <: AbstractInteraction
    J::Mat3
    bond::Bond
    label::String
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
    biquadratic(B, bond::Bond, label::String="BHeisen")

Creates a Biquadratic interaction
```math
    B ∑_{⟨ij⟩} (𝐒_i ⋅ 𝐒_j)^2
```
where ``⟨ij⟩`` runs over all bonds symmetry equivalent to `bond`.
"""
function biquadratic(B, bond::Bond, label::String="BHeisen")
    BiQuadraticInteraction(B * Mat3(I), bond, label)
end


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
