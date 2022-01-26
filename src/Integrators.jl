""" Structs and methods for dynamics integrators.
    Standard (non-damped) Landau-Lifshitz dynamics integrates:
        dS/dt = -S × B,
    where `B = -∇ H` is the effective field.
"""

abstract type Integrator end

"""
    HeunP(sys::SpinSystem)

Integrates Landau-Lifshitz spin dynamics using the Heun method, with a final
projection step that exactly constrains |S|=1. The method is locally second
order accurate.
"""
mutable struct HeunP <: Integrator
    sys :: SpinSystem
    _S₁ :: Array{Vec3, 4}
    _S₂ :: Array{Vec3, 4}
    _B  :: Array{Vec3, 4}
    _f₁ :: Array{Vec3, 4}

    function HeunP(sys::SpinSystem)
        return new(
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
mutable struct LangevinHeunP <: Integrator
    α   :: Float64        # Damping coeff
    D   :: Float64        # Stochastic strength
    sys :: SpinSystem
    _S₁ :: Array{Vec3, 4} # Intermediate integration variable space
    _S₂ :: Array{Vec3, 4}
    _B  :: Array{Vec3, 4}
    _f₁ :: Array{Vec3, 4}
    _r₁ :: Array{Vec3, 4}
    _ξ  :: Array{Vec3, 4}

    function LangevinHeunP(sys::SpinSystem, kT::Float64, α::Float64)
        return new(
            α, α*kT/(1+α*α), sys,
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
mutable struct SphericalMidpoint <: Integrator
    sys :: SpinSystem
    _S̄  :: Array{Vec3, 4}
    _Ŝ  :: Array{Vec3, 4}
    _S̄′ :: Array{Vec3, 4}
    _B  :: Array{Vec3, 4}
    atol :: Float64

    function SphericalMidpoint(sys::SpinSystem; atol=1e-12)
        return new(
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
    # Normalize fluctuations by 1/√S
    for b in 1:nbasis(sys)
        selectdim(_ξ, 1, b) .*= 1. / sqrt(sys.sites_info[b].S)
    end

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
mutable struct LangevinSampler <: AbstractSampler
    integrator :: LangevinHeunP
    kT         :: Float64
    Δt         :: Float64
    nsteps     :: Int
end

"""
    LangevinSampler(sys, kT, α, Δt, nsteps)

Creates a `LangevinSampler` which samples the spin system's Hamiltonian using Langevin
 dynamics at a temperature `kT`, damping coefficient `α`, and producing a new sample
 by integrating with `nsteps` timesteps of size `Δt`.
"""
function LangevinSampler(sys::SpinSystem, kT::Float64, α::Float64, Δt::Float64, nsteps::Int)
    integrator = LangevinHeunP(sys, kT, α)
    LangevinSampler(integrator, kT, Δt, nsteps)
end

@inline function set_temp!(sampler::LangevinSampler, kT::Float64)
    sampler.kT = kT
    α = sampler.integrator.α
    sampler.integrator.D = α * kT / (1 + α * α)
    nothing
end

get_temp(sampler::LangevinSampler) = sampler.kT
get_system(sampler::LangevinSampler) = sampler.integrator.system

@inline function sample!(sampler::LangevinSampler)
    for _ in 1:sampler.nsteps
        evolve!(sampler.integrator, sampler.Δt)
    end
end
