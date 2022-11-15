""" Structs and methods for dynamics integrators.
    Standard (non-damped) Landau-Lifshitz dynamics integrates:
        dS/dt = -S × B,
    where `B = -∇ H` is the effective field.
"""

abstract type Integrator end
abstract type LangevinIntegrator <: Integrator end

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
mutable struct LangevinHeunP <: LangevinIntegrator
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
            α, α*kT, sys,
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

    function SphericalMidpoint(sys::SpinSystem{0}; atol=1e-12)
        return new(
            sys, zero(sys._dipoles), zero(sys._dipoles),
            zero(sys._dipoles), zero(sys._dipoles), atol
        )
    end    
end

function SphericalMidpoint(sys::SpinSystem{N}; atol=1e-12) where N
    error("SphericalMidpoint integrator is not available for SU(N) systems. Use ImplicitMidpoint integrator.")
    nothing
end

@doc raw"""
    ImplicitMidpoint(sys::SpinSystem{N}) where N

Returns an implicit midpoint integrator a `SpinSystem{N}`, where N may be any integer. If
`N=0` (traditional Landua-Lifshitz), it returns a `SphericalMidpoint` integrator,
otherwise a `SchrodingerMidpoint` integrator. Both integrators are sympletic, i.e. they
guarantee that energy error is bounded.

For dissipative dynamics in a themal bath, use the `LangevinHeunP` integrator.
"""
# Not exposing atol option here. If user wants to fiddle with the tolerance, they'LL
# have to use the actual integrator interface (i.e., SchrodingerMidpoint or SphericalMidpoint.)
function ImplicitMidpoint(sys::SpinSystem{N}) where N
    N == 0 ? SphericalMidpoint(sys) : SchrodingerMidpoint(sys)
end


mutable struct LangevinHeunPSUN{N} <: LangevinIntegrator
    kT    :: Float64
    α     :: Float64
    sys   :: SpinSystem{N}
    _ℌZ   :: Array{CVec{N}, 4}
    _Z′   :: Array{CVec{N}, 4}
    _ΔZ₁  :: Array{CVec{N}, 4}
    _ΔZ₂  :: Array{CVec{N}, 4}
    _ξ    :: Array{CVec{N}, 4}
    _B    :: Array{Vec3, 4}
    _ℌ    :: Matrix{ComplexF64}  # Just a buffer
end

"""
    LangevinHeunPSUN(sys::SpinSystem{N}, kT::Float64, α::Float64)

Implements Langevin dynamics on `sys` targeting a temperature `kT`, with a
damping coefficient `α`. This integrator is only for Generalized Spin Dynamics.
If sys has type `SpinSystem{0}`, will revert to the `LangevinHeunP` integrator.

Uses the 2nd-order Heun + projection scheme."
"""
function LangevinHeunPSUN(sys::SpinSystem{N}, kT::Float64, α::Float64) where N
    LangevinHeunPSUN{N}(
        kT, α, sys,
        zero(sys._coherents), zero(sys._coherents),
        zero(sys._coherents), zero(sys._coherents),
        zero(sys._coherents), zero(sys._dipoles),
        zeros(ComplexF64, (N,N))
    )
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
    _ℌ   :: Matrix{ComplexF64}  # Just a buffer
    atol :: Float64
end

function SchrodingerMidpoint(sys::SpinSystem{N}; atol=1e-14) where N
    SchrodingerMidpoint{N}(
        sys, 
        zero(sys._coherents), zero(sys._coherents),
        zero(sys._coherents), zero(sys._coherents),
        zero(sys._coherents), zero(sys._dipoles),
        zeros(ComplexF64, (N,N)), atol
    )
end

function SchrodingerMidpoint(sys::SpinSystem{0}; atol=1e-14)
    error("SchrodingerMidpoint integration is only available for SU(N) systems. Use ImplicitMidpoint integrator.")
    nothing
end


@inline f(S, B) = -S × B
@inline f(S, B, α) = -S × (B + α * (S × B))

# Normalize to κ value given in site_infos. For old LL dynamics only.
function normalize!(S::Array{Vec3, 4}, sys::SpinSystem)
    for site in 1:size(S, 4)
        spin_rescaling = sys.site_infos[site].spin_rescaling
        for cell in CartesianIndices(size(S)[1:3]) 
            S[cell, site] *= spin_rescaling/norm(S[cell, site])
        end
    end
end

normalize_dipoles!(sys::SpinSystem) = normalize!(sys._dipoles, sys)

"""
    evolve!(integrator, Δt)

Performs a single integrator timestep of size Δt.
"""
function evolve!(integrator::HeunP, Δt::Float64)
    (; sys, _S₁, _B, _f₁) = integrator
    S, ℋ = sys._dipoles, sys.hamiltonian
    
    # Euler step
    field!(_B, S, ℋ)
    @. _f₁ = f(S, _B)
    @. _S₁ = S + Δt * _f₁

    # Corrector step
    field!(_B, _S₁, ℋ)
    @. S = S + 0.5 * Δt * (_f₁ + f(_S₁, _B))
    normalize!(S, sys)

    nothing
end

function evolve!(integrator::LangevinHeunP, Δt::Float64)
    (; α, D, sys, _S₁, _B, _f₁, _r₁, _ξ) = integrator
    S, ℋ = sys._dipoles, sys.hamiltonian

    randn!(sys.rng, _ξ)
    _ξ .*= √(2D)

    # Euler step
    field!(_B, S, ℋ)
    @. _f₁ = f(S, _B, α)
    @. _r₁ = f(S, _ξ)   # no absence of α argument -- noise only appears once in rhs.
    @. _S₁ = S + Δt * _f₁ + √Δt * _r₁

    # Corrector step
    field!(_B, _S₁, ℋ)
    @. S = S + 0.5 * Δt * (_f₁ + f(_S₁, _B, α)) + 0.5 * √Δt * (_r₁ + f(_S₁, _ξ))
    normalize!(S, sys)

    nothing
end




function evolve!(integrator::SphericalMidpoint, Δt::Float64)
    (; sys, _S̄, _Ŝ, _S̄′, _B, atol) = integrator
    S, ℋ = sys._dipoles, sys.hamiltonian
    
    # Initial guess for midpoint
    @. _S̄ = S

    max_iters = 100
    for iter in 1:max_iters
        # Integration step for current best guess of midpoint _S̄. Produces
        # improved midpoint estimator _S̄′.
        @. _Ŝ =_S̄
        normalize!(_Ŝ, sys)
        field!(_B, _Ŝ, ℋ)
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
            @. S = 2*_S̄′ - S
            normalize!(S, sys)
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
    @inbounds for i in eachindex(Z)
        Z[i] = Z[i] / √(real(Z[i]' * Z[i]))
    end
end


@inline function _proj(a::T, Z::T)  where T <: CVec
    (a - ((Z' * a) * Z))  
end


function add_dipolar_field!(op::Array{ComplexF64, 2}, B::Sunny.Vec3)
    N = size(op, 1)
    S = (N-1)/2

    # Note indexing by column (hence using conjugate of standard formula)
    @inbounds for j in 1:N
        # Subdiagonal
        if j > 1
            val = 0.5*√(S*(S + 1) - (S - j + 1)*(S - j + 2))
            op[j-1,j] += val*(B[1] - im*B[2])
        end

        # Diagonal
        op[j,j] += B[3]*(S - j + 1)

        # Superdiagonal
        if j < N
            val = 0.5*√(S*(S + 1) - (S - j + 1)*(S - j))
            op[j+1,j] += val*(B[1] + im*B[2])
        end
    end

    nothing
end

@generated function _apply_ℌ!(rhs::Array{CVec{N}, 4}, B::Array{Vec3, 4}, Z::Array{CVec{N}, 4}, integrator)  where {N}

    if integrator <: LangevinIntegrator
        scale_expr = :(site_infos[s].spin_rescaling)
    else
        scale_expr = :(1.0)
    end

    if N < 6 
        S = spin_matrices(N) .|> SMatrix{N, N, ComplexF64, N*N}

        return quote
            (; sys) = integrator
            (; hamiltonian, site_infos, lattice) = sys
            nb = nbasis(lattice)
            Λs′ = hamiltonian.sun_aniso
            Λs = reinterpret(SMatrix{N, N, ComplexF64, N*N},
                             reshape(Λs′, N*N, nb)
            )
            @inbounds for s in 1:nb
                κ = $scale_expr
                for c in eachcellindex(lattice)
                    ℌ = κ * (Λs[s] - ($S)' * B[c,s])
                    rhs[c, s] = ℌ*Z[c, s]
                end
            end
            nothing
        end
    end

    return quote
        (; sys, _ℌ) = integrator
        (; hamiltonian, site_infos, lattice) = sys
        nb = nbasis(lattice)
        Λs = hamiltonian.sun_aniso
        rhs′ = reinterpret(reshape, ComplexF64, rhs) 

        @inbounds for s in 1:nb
            κ = $scale_expr 
            Λ = @view(Λs[:,:,s])
            for c in eachcellindex(lattice)
                @. _ℌ = κ * Λ 
                add_dipolar_field!(_ℌ, -κ*B[c,s]) 
                mul!(@view(rhs′[:, c, s]), _ℌ, Z[c, s])
            end
        end
        nothing
    end
end


function _rhs_langevin!(ΔZ::Array{CVec{N}, 4}, Z::Array{CVec{N}, 4}, integrator::LangevinHeunPSUN, Δt::Float64) where N
    (; kT, α, sys, _ℌZ, _ξ, _B, _ℌ) = integrator
    (; _dipoles) = sys

    set_expected_spins!(_dipoles, Z, sys) 
    field!(_B, _dipoles, sys.hamiltonian)
    _apply_ℌ!(_ℌZ, _B, Z, integrator)

    @inbounds for i in eachindex(Z)
        ΔZ′ = -im*√(2*Δt*kT*α)*_ξ[i] - Δt*(im+α)*_ℌZ[i]    
        ΔZ[i] = _proj(ΔZ′, Z[i])
    end 
    nothing
end


function _rhs_ll!(ΔZ, Z, integrator, Δt)
    (; sys, _ℌZ, _B, _ℌ) = integrator
    (; _dipoles) = sys

    set_expected_spins!(_dipoles, Z, sys) # temporarily de-synchs _dipoles and _coherents
    field!(_B, _dipoles, sys.hamiltonian)
    _apply_ℌ!(_ℌZ, _B, Z, integrator)

    @inbounds for i in eachindex(Z)
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


function evolve!(integrator::SchrodingerMidpoint, Δt::Float64; max_iters=100)
    (; sys, _ΔZ, _Zb, _Z′, _Z″, atol) = integrator
    Z = sys._coherents

    @. _Z′ = Z 
    @. _Z″ = Z 

    for _ in 1:max_iters
        @. _Zb = (Z + _Z′)/2

        _rhs_ll!(_ΔZ, _Zb, integrator, Δt)

        @. _Z″ = Z + _ΔZ

        if norm(_Z′ - _Z″) < atol    # NOTE: This causes allocations.
            @. Z = _Z″
            set_expected_spins!(sys)
            return
        end

        _Z′, _Z″ = _Z″, _Z′
    end

    error("SchrodingerMidpoint integrator failed to converge in $max_iters iterations.")
end



#= Need to figure out a better way to dispatch on LangevinSamplers. Currently
use a Union type to avoid code repetition.

One solution, still kind of ugly: parameterize LangevinSampler with N,
and store N in a dummy variable. Then have the integrator be a Union
type of LangevinHeunP and LangevinHeunPSUN. =#

"""
    A sampler which produces new samples using Langevin Landau-Lifshitz dynamics
"""
mutable struct LangevinSamplerLLD <: AbstractSampler
    integrator :: LangevinHeunP
    kT         :: Float64
    Δt         :: Float64
    nsteps     :: Int
end

mutable struct LangevinSamplerGSD <: AbstractSampler
    integrator :: LangevinHeunPSUN
    kT         :: Float64
    Δt         :: Float64
    nsteps     :: Int
end

LangevinSamplerT = Union{LangevinSamplerLLD, LangevinSamplerGSD}

"""
    LangevinSampler(sys::SpinSystem, kT::Float64, α::Float64, Δt::Float64, nsteps::Int)

Creates a `LangevinSampler` which samples the spin system's Hamiltonian using Langevin
 dynamics at a temperature `kT`, damping coefficient `α`, and producing a new sample
 by integrating with `nsteps` timesteps of size `Δt`.
"""
function LangevinSampler(sys::SpinSystem{0}, kT::Float64, α::Float64, Δt::Float64, nsteps::Int)
    integrator = LangevinHeunP(sys, kT, α)
    LangevinSamplerLLD(integrator, kT, Δt, nsteps)
end
function LangevinSampler(sys::SpinSystem{N}, kT::Float64, α::Float64, Δt::Float64, nsteps::Int) where N
    integrator = LangevinHeunPSUN(sys, kT, α)
    LangevinSamplerGSD(integrator, kT, Δt, nsteps)
end

@inline function set_temp!(sampler::LangevinSamplerLLD, kT::Float64)
    sampler.kT = kT
    α = sampler.integrator.α
    sampler.integrator.D = α * kT
    nothing
end
@inline function set_temp!(sampler::LangevinSamplerGSD, kT::Float64)
    sampler.integrator.kT = sampler.kT = kT
    nothing
end

get_temp(sampler::LangevinSamplerT)  = sampler.kT
get_system(sampler::LangevinSamplerT)  = sampler.integrator.sys
Random.rand!(sampler::LangevinSamplerT) = rand!(sampler.integrator.sys)

@inline function sample!(sampler::LangevinSamplerT) 
    for _ in 1:sampler.nsteps
        evolve!(sampler.integrator, sampler.Δt)
    end
end


