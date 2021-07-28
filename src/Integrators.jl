""" Structs and methods for dynamics integrators -- attempted to stay as close as possible
     to the types of David's, but didn't make it parametric over the Float type (T)
"""

""" Standard (non-damped) Landau-Lifshitz dynamics integrates:
      dS/dt = -S × H, where H = B + ∑_neigh S_neigh
"""

abstract type Integrator end

"Integrator for a 2nd-order energy-conserving Heun + projection scheme"
mutable struct HeunP{D, L, Db} <: Integrator
    sys :: SpinSystem{D, L, Db}
    _S₁ :: Array{Vec3, Db}
    _S₂ :: Array{Vec3, Db}
    _B  :: Array{Vec3, Db}
    _f₁ :: Array{Vec3, Db}
end

"Integrator implementing Langevin dynamics using a 2nd-order Heun + projection scheme"
mutable struct LangevinHeunP{D, L, Db} <: Integrator
    α   :: Float64                # Normalized damping, α/S
    D   :: Float64                # Langevin stochastic coefficient
    sys :: SpinSystem{D, L, Db}
    _S₁ :: Array{Vec3, Db}        # Intermediate integration variable space
    _S₂ :: Array{Vec3, Db}
    _B  :: Array{Vec3, Db}
    _f₁ :: Array{Vec3, Db}
    _r₁ :: Array{Vec3, Db}
    _ξ  :: Array{Vec3, Db}
end

"""
    HeunP(sys)
"""
function HeunP(sys::SpinSystem)
    return HeunP(
        sys, zero(sys.sites), zero(sys.sites),
        zero(sys.sites), zero(sys.sites),
    )
end

"""
    LangevinHeunP(sys, kT, α)

Implements Langevin dynamics on `sys` targetting a temperature `kT`,
 with a damping coefficient `α`. Provided `α` should not be normalized
 by the spin magnitude -- this is done internally.
"""
function LangevinHeunP(sys::SpinSystem, kT::Float64, α::Float64)
    # D should also have a factor of 1/S, the magnitude of the spin. Assumed here 1.
    return LangevinHeunP(
        α/sys.S, α*kT/((1+α*α)*sys.S), sys,
        zero(sys.sites), zero(sys.sites), zero(sys.sites),
        zero(sys.sites), zero(sys.sites), zero(sys.sites)
    )
end


@inline f(S, B) = -S × B
@inline f(S, B, α) = -S × (B + α * (S × B))

"""
    evolve!(integrator, Δt)

Performs a single integrator timestep of size Δt.
"""
function evolve!(integrator::HeunP, Δt::Float64)
    @unpack sys, _S₁, _S₂, _B, _f₁ = integrator
    S = sys.sites
    
    # Euler step
    field!(_B, S, sys.interactions)
    @. _f₁ = f(S, _B)
    @. _S₁ = S + Δt * _f₁

    # Corrector step
    field!(_B, _S₁, sys.interactions)
    @. _S₂ = S + 0.5 * Δt * (_f₁ + f(_S₁, _B))
    @. _S₂ /= norm(_S₂)

    # Swap buffers
    sys.sites, integrator._S₂ = integrator._S₂, sys.sites
    nothing
end

function evolve!(integrator::LangevinHeunP, Δt::Float64)
    @unpack α, D, sys, _S₁, _S₂, _B, _f₁, _r₁, _ξ = integrator
    S = sys.sites

    Random.randn!(_ξ)
    _ξ .*= √(2D)

    # Euler step
    field!(_B, S, sys.interactions)
    @. _f₁ = f(S, _B, α)
    @. _r₁ = f(S, _ξ, α)
    @. _S₁ = S + Δt * _f₁ + √Δt * _r₁

    # Corrector step
    field!(_B, _S₁, sys.interactions)
    @. _S₂ = S + 0.5 * Δt * (_f₁ + f(_S₁, _B, α)) + 0.5 * √Δt * (_r₁ + f(_S₁, _ξ, α))
    @. _S₂ /= norm(_S₂)

    # Swap buffers
    sys.sites, integrator._S₂ = integrator._S₂, sys.sites
    nothing
end