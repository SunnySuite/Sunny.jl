""" Structs and methods for dynamics integrators -- attempted to stay as close as possible
     to the types of David's, but didn't make it parametric over the Float type (T)
"""

""" Standard (non-damped) Landau-Lifshitz dynamics integrates:
      dS/dt = -S × H, where H = B + ∑_neigh S_neigh
"""

abstract type Integrator end

"""
    HeunP(sys::SpinSystem)

Integrates Landau-Lifshitz spin dynamics using the Heun method, with a final
projection step that exactly constrains |S|=1. The method is locally second
order accurate.
"""
mutable struct HeunP{D, L, Db} <: Integrator
    sys :: SpinSystem{D, L, Db}
    _S₁ :: Array{Vec3, Db}
    _S₂ :: Array{Vec3, Db}
    _B  :: Array{Vec3, Db}
    _f₁ :: Array{Vec3, Db}

    function HeunP(sys::SpinSystem{D, L, Db}) where {D, L, Db}
        return new{D, L, Db}(
            sys, zero(sys.sites), zero(sys.sites),
            zero(sys.sites), zero(sys.sites)
        )
    end
end

"""
    LangevinHeunP(sys, kT, α)

Implements Langevin dynamics on `sys` targetting a temperature `kT`, with a
damping coefficient `α`. Provided `α` should not be normalized by the spin
magnitude -- this is done internally.

Uses the 2nd-order Heun + projection scheme."
"""
mutable struct LangevinHeunP{D, L, Db} <: Integrator
    α   :: Float64                # Damping coeff normalized by S
    D   :: Float64                # Stochastic strength normalized by S
    sys :: SpinSystem{D, L, Db}
    _S₁ :: Array{Vec3, Db}        # Intermediate integration variable space
    _S₂ :: Array{Vec3, Db}
    _B  :: Array{Vec3, Db}
    _f₁ :: Array{Vec3, Db}
    _r₁ :: Array{Vec3, Db}
    _ξ  :: Array{Vec3, Db}

    function LangevinHeunP(sys::SpinSystem{D, L, Db}, kT::Float64, α::Float64) where {D, L, Db}
        return new{D, L, Db}(
            α/sys.S, α*kT/((1+α*α)*sys.S), sys,
            zero(sys.sites), zero(sys.sites), zero(sys.sites),
            zero(sys.sites), zero(sys.sites), zero(sys.sites)
        )
    end
end

"""
    SphericalMidpoint(sys::SpinSystem; atol=1e-12)

Integrates Landau-Lifshitz spin dynamics using the spherical-midpoint
integrator, which is symplectic and implicit. Each step is converged to absolute
tolerance `atol`.
"""
mutable struct SphericalMidpoint{D, L, Db} <: Integrator
    sys :: SpinSystem{D, L, Db}
    _S̄  :: Array{Vec3, Db}
    _Ŝ  :: Array{Vec3, Db}
    _S̄′ :: Array{Vec3, Db}
    _B  :: Array{Vec3, Db}
    atol :: Float64

    function SphericalMidpoint(sys::SpinSystem{D, L, Db}; atol=1e-12) where {D, L, Db}
        return new{D, L, Db}(
            sys, zero(sys.sites), zero(sys.sites),
            zero(sys.sites), zero(sys.sites), atol
        )
    end    
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
    field!(_B, S, sys.hamiltonian)
    @. _f₁ = f(S, _B)
    @. _S₁ = S + Δt * _f₁

    # Corrector step
    field!(_B, _S₁, sys.hamiltonian)
    @. _S₂ = normalize(S + 0.5 * Δt * (_f₁ + f(_S₁, _B)))

    # Swap buffers
    sys.sites, integrator._S₂ = integrator._S₂, sys.sites
    nothing
end

function evolve!(integrator::LangevinHeunP, Δt::Float64)
    @unpack α, D, sys, _S₁, _S₂, _B, _f₁, _r₁, _ξ = integrator
    S = sys.sites

    randn!(_ξ)
    _ξ .*= √(2D)

    # Euler step
    field!(_B, S, sys.hamiltonian)
    @. _f₁ = f(S, _B, α)
    @. _r₁ = f(S, _ξ, α)
    @. _S₁ = S + Δt * _f₁ + √Δt * _r₁

    # Corrector step
    field!(_B, _S₁, sys.hamiltonian)
    @. _S₂ = normalize(S + 0.5 * Δt * (_f₁ + f(_S₁, _B, α)) + 0.5 * √Δt * (_r₁ + f(_S₁, _ξ, α)))

    # Swap buffers
    sys.sites, integrator._S₂ = integrator._S₂, sys.sites
    nothing
end

function evolve!(integrator::SphericalMidpoint, Δt::Float64)
    @unpack sys, _S̄, _Ŝ, _S̄′, _B, atol = integrator
    S = sys.sites
    
    # Initial guess for midpoint
    @. _S̄ = S

    max_iters = 100
    for iter in 1:max_iters
        # Integration step for current best guess of midpoint _S̄. Produces
        # improved midpoint estimator _S̄′.
        @. _Ŝ = normalize(_S̄)
        field!(_B, _Ŝ, sys.hamiltonian)
        @. _S̄′ = S + 0.5 * Δt * f(_Ŝ, _B)

        # Convergence is reached if every element of _S̄ and _S̄′ agree to
        # within tolerance.
        converged = true
        for i = 1:length(S)
            converged = converged && isapprox(_S̄[i], _S̄′[i]; atol)
        end

        # If converged, then we can return
        if converged
            # Normalization here should not be necessary in principle, but it
            # could be useful in practice for finite `atol`.
            @. S = normalize(2*_S̄′ - S)
            return
        end

        @. _S̄ = _S̄′
    end

    error("Spherical midpoint method failed to converge to tolerance $atol after $max_iters iterations.")
end


"""
    A sampler which produces new samples using Langevin Landau-Lifshitz dynamics
"""
struct LangevinSampler{D, L, Db} <: AbstractSampler
    integrator :: LangevinHeunP{D, L, Db}
    Δt         :: Float64
    nsteps     :: Int
end

"""
    LangevinSampler(sys, kT, α, Δt, nsteps)

Creates a `LangevinSampler` which samples the spin system's Hamiltonian using Langevin
 dynamics at a temperature `kT`, damping coefficient `α`, and producing a new sample
 by integrating with `nsteps` timesteps of size `Δt`.
"""
function LangevinSampler(sys::SpinSystem{D, L, Db}, kT::Float64, α::Float64, Δt::Float64, nsteps::Int) where {D, L, Db}
    integrator = LangevinHeunP(sys, kT, α)
    LangevinSampler{D, L, Db}(integrator, Δt, nsteps)
end

@inline function set_temp!(sampler::LangevinSampler, kT::Float64)
    α = sampler.integrator.α
    S = sampler.integrator.sys.S
    sampler.integrator.D = α * kT / ((1 + α * α) * S)
    nothing
end

@inline function sample!(sampler::LangevinSampler)
    for _ in 1:sampler.nsteps
        evolve!(sampler.integrator, sampler.Δt)
    end
end

@inline running_energy(sampler::LangevinSampler) = energy(sampler.integrator.sys)
@inline running_mag(sampler::LangevinSampler) = sum(sampler.integrator.sys)
