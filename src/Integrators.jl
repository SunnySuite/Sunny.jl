""" Structs and methods for dynamics integrators.
    Standard (non-damped) Landau-Lifshitz dynamics integrates:
        dS/dt = -S × B,
    where `B = -∇ H` is the effective field.
"""

abstract type Integrator end

# TODO: These integrators should verify that the system is in
#        "dipole-mode". Or, is there some nice way to make them
#        work with both? Using e.g. DipoleView.
"""
    HeunP(sys::SpinSystem)

Integrates Landau-Lifshitz spin dynamics using the Heun method, with a final
projection step that exactly constrains |S|=1. The method is locally second
order accurate.
"""
mutable struct HeunP <: Integrator
    sys :: SpinSystem
    _S₁ :: Array{Vec3, 4}
    _B  :: Array{Vec3, 4}
    _f₁ :: Array{Vec3, 4}

    function HeunP(sys::SpinSystem)
        return new(
            sys, zero(sys._dipoles), zero(sys._dipoles), zero(sys._dipoles)
        )
    end
end

"""
    LangevinHeunP(sys, kT, α)

Implements Langevin dynamics on `sys` targeting a temperature `kT`, with a
damping coefficient `α`. Provided `α` should not be normalized by the spin
magnitude -- this is done internally.

Uses the 2nd-order Heun + projection scheme."
"""
mutable struct LangevinHeunP <: Integrator
    α   :: Float64        # Damping coeff
    D   :: Float64        # Stochastic strength
    sys :: SpinSystem
    _S₁ :: Array{Vec3, 4} # Intermediate integration variable space
    _B  :: Array{Vec3, 4}
    _f₁ :: Array{Vec3, 4}
    _r₁ :: Array{Vec3, 4}
    _ξ  :: Array{Vec3, 4}

    function LangevinHeunP(sys::SpinSystem{0}, kT::Float64, α::Float64)
        return new(
            α, α*kT/(1+α*α), sys,
            zero(sys._dipoles), zero(sys._dipoles), zero(sys._dipoles),
            zero(sys._dipoles), zero(sys._dipoles)
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
            sys, zero(sys._dipoles), zero(sys._dipoles),
            zero(sys._dipoles), zero(sys._dipoles), atol
        )
    end    
end


"""
    LangevinHeunPSUN(sys::SpinSystem{N}, kT::Float64, α::Float64)

Implements Langevin dynamics on `sys` targeting a temperature `kT`, with a
damping coefficient `α`. This integrator is only for Generalized Spin Dynamics.
If sys has type `SpinSystem{0}`, will revert to the `LangevinHeunP` integrator.

Uses the 2nd-order Heun + projection scheme."
"""
mutable struct LangevinHeunPSUN{N} <: Integrator
    kT  :: Float64
    α   :: Float64
    sys :: SpinSystem{N}
    _ℌZ  :: Array{CVec{N}, 4}
    _Z′  :: Array{CVec{N}, 4}
    _ΔZ₁ :: Array{CVec{N}, 4}
    _ΔZ₂ :: Array{CVec{N}, 4}
    _ξ   :: Array{CVec{N}, 4}
    _B   :: Array{Vec3, 4}

    function LangevinHeunPSUN(sys::SpinSystem{N}, kT::Float64, α::Float64) where N
        new{N}(
            kT, α, sys,
            zero(sys._coherents), zero(sys._coherents),
            zero(sys._coherents), zero(sys._coherents),
            zero(sys._coherents), zero(sys._dipoles)
        )
    end
end

# Constructor so SU(N) integrator can be used with identical interface to LL LangevinHeunP integrator
function LangevinHeunP(sys::SpinSystem{N}, kT::Float64, α::Float64) where N
    LangevinHeunPSUN(sys, kT, α)
end


"""
    SchrodingerMidpoint(sys::SpinSystem{N}) where N

Implements Generalized Spin Dynamics using the Schrodinger Midpoint method. No
noise or damping. Only works for SU(N) systems (N != 0). Will revert to SphericalMidpoint
method if sys has a type of `SpinSystem{0}`.
"""
mutable struct SchrodingerMidpoint{N} <: Integrator
    sys :: SpinSystem{N}
    _ΔZ  :: Array{CVec{N}, 4}
    _ℌZ  :: Array{CVec{N}, 4}
    _Zb  :: Array{CVec{N}, 4}
    _Z′  :: Array{CVec{N}, 4}
    _Z″  :: Array{CVec{N}, 4}
    _B   :: Array{Vec3, 4}

    function SchrodingerMidpoint(sys::SpinSystem{N}) where N
        new{N}(
            sys, 
            zero(sys._coherents), zero(sys._coherents),
            zero(sys._coherents), zero(sys._coherents),
            zero(sys._coherents), zero(sys._dipoles)
        )
    end
end

function SchrodingerMidpoint(sys::SpinSystem{0})
    @warn "SchrodingerMidpoint integration is only available for SU(N) systems. Reverting to SphericalMidpoint integrator."
    SphericalMidpoint(sys)
end


@inline f(S, B) = -S × B
@inline f(S, B, α) = -S × (B + α * (S × B))

"""
    evolve!(integrator, Δt)

Performs a single integrator timestep of size Δt.
"""
function evolve!(integrator::HeunP, Δt::Float64)
    (; sys, _S₁, _B, _f₁) = integrator
    S, Z = sys._dipoles, sys._coherents
    
    # Euler step
    field!(_B, S, Z, sys.hamiltonian)
    @. _f₁ = f(S, _B)
    @. _S₁ = S + Δt * _f₁

    # Corrector step
    field!(_B, _S₁, Z, sys.hamiltonian)
    @. S = normalize(S + 0.5 * Δt * (_f₁ + f(_S₁, _B)))

    nothing
end

function evolve!(integrator::LangevinHeunP, Δt::Float64)
    (; α, D, sys, _S₁, _B, _f₁, _r₁, _ξ) = integrator
    S, Z = sys._dipoles, sys._coherents

    randn!(sys.rng, _ξ)
    _ξ .*= √(2D)
    # Normalize fluctuations by 1/√S
    for b in 1:nbasis(sys)
        selectdim(_ξ, 1, b) .*= 1. / sqrt(sys.site_infos[b].κ)
    end

    # Euler step
    field!(_B, S, Z, sys.hamiltonian)
    @. _f₁ = f(S, _B, α)
    @. _r₁ = f(S, _ξ, α)
    @. _S₁ = S + Δt * _f₁ + √Δt * _r₁

    # Corrector step
    field!(_B, _S₁, Z, sys.hamiltonian)
    @. S = normalize(S + 0.5 * Δt * (_f₁ + f(_S₁, _B, α)) + 0.5 * √Δt * (_r₁ + f(_S₁, _ξ, α)))

    nothing
end




function evolve!(integrator::SphericalMidpoint, Δt::Float64)
    @unpack sys, _S̄, _Ŝ, _S̄′, _B, atol = integrator
    S, Z = sys._dipoles, sys._coherents
    
    # Initial guess for midpoint
    @. _S̄ = S

    max_iters = 100
    for iter in 1:max_iters
        # Integration step for current best guess of midpoint _S̄. Produces
        # improved midpoint estimator _S̄′.
        @. _Ŝ = normalize(_S̄)
        field!(_B, _Ŝ, Z, sys.hamiltonian)
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


""" SU(N) integrators and helper functions.
"""

# LinearAlgebra's `normalize!` doesn't normalize to 1 on complex vectors
function normalize!(Z::Array{CVec{N}, 4}) where N
    for i in eachindex(Z)
        Z[i] = Z[i] / √(real(Z[i]' * Z[i]))
    end
end


@inline function _proj(a::T, Z::T)  where T <: CVec
    (a - ((Z' * a) * Z))  
end


function _apply_ℌ!(rhs, sys, B, Z)
    aniso = sys.hamiltonian.sun_aniso
    latsize = size(B)[2:end]
    for (site, Λ) in zip(aniso.sites, aniso.Λs)  # Need to make sure only one aniso per site
        for cell in CartesianIndices(latsize)
             rhs[site, cell] = (-B[site, cell] ⋅ sys.S + Λ) * Z[site, cell] # ∇E = -B
        end
    end
end


function _rhs_langevin!(ΔZ, Z, integrator, Δt)
    (; kT, α, sys, _ℌZ, _ξ, _B) = integrator
    (; _dipoles) = sys

    set_expected_spins!(_dipoles, Z) 
    field!(_B, _dipoles, Z, sys.hamiltonian)
    _apply_ℌ!(_ℌZ, sys, _B, Z)

    for i in eachindex(Z)
        ΔZ′ = -im*√(2*Δt*kT*α)*_ξ[i] - Δt*(im+α)*_ℌZ[i]    
        ΔZ[i] = _proj(ΔZ′, Z[i])
    end 
end


function _rhs_ll!(ΔZ, Z, integrator, Δt)
    (; sys, _ℌZ, _B) = integrator
    (; _dipoles) = sys

    set_expected_spins!(_dipoles, Z) # temporarily de-synchs _dipoles and _coherents
    field!(_B, _dipoles, Z, sys.hamiltonian)
    _apply_ℌ!(_ℌZ, sys, _B, Z)

    for i in eachindex(Z)
        ΔZ[i] = - Δt*im*_ℌZ[i]    
    end 
end


function evolve!(integrator::LangevinHeunPSUN, Δt::Float64)
    (; sys, _Z′, _ΔZ₁, _ΔZ₂, _ξ) = integrator
    Z = sys._coherents

    randn!(sys.rng, _ξ)

    # Prediction
    _rhs_langevin!(_ΔZ₁, Z, integrator, Δt)
    @. _Z′ = Z + _ΔZ₁
    normalize!(_Z′)

    # Correction
    _rhs_langevin!(_ΔZ₂, _Z′, integrator, Δt)
    @. Z += (_ΔZ₁ + _ΔZ₂)/2
    normalize!(Z)

    # Coordinate dipole data
    set_expected_spins!(sys)
end


function evolve!(integrator::SchrodingerMidpoint, Δt::Float64; tol=1e-14, max_iters=100)
    (; sys, _ΔZ, _Zb, _Z′, _Z″) = integrator
    Z = sys._coherents

    @. _Z′ = Z 
    @. _Z″ = Z 

    for i in 1:max_iters
        @. _Zb = (Z + _Z′)/2

        _rhs_ll!(_ΔZ, _Zb, integrator, Δt)

        @. _Z″ = Z + _ΔZ

        if norm(_Z′ - _Z″) < tol
            @. Z = _Z″
            set_expected_spins!(sys)
            return
        end

        _Z′, _Z″ = _Z″, _Z′
    end

    error("SchrodingerMidpoint integrator failed to converge in $max_iters iterations.")
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


