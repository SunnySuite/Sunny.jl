"""
    Langevin(Δt::Float64; λ::Float64, kT::Float64)

Spin dynamics with damping and noise terms that model coupling to an implicit
thermal bath, of strength `λ`. One call to the [`step!`](@ref) function will
advance a [`System`](@ref) by `Δt` units of time. Can be used to sample from the
Boltzmann distribution at temperature `kT`. An alternative approach to sampling
states from thermal equilibrium is [`LocalSampler`](@ref), which proposes local
Monte Carlo moves. For example, use `LocalSampler` to sample Ising-like spins.

Setting `λ = 0` disables coupling to the thermal bath, yielding an
energy-conserving spin dynamics. The `Langevin` integrator uses an explicit
numerical integrator that allows energy drift. Alternatively, the
[`ImplicitMidpoint`](@ref) method can be used, which is more expensive but
prevents energy drift through exact conservation of the symplectic 2-form.

If the [`System`](@ref) has `mode = :dipole`, then the dynamics is the
stochastic Landau-Lifshitz equation,
```math
    d𝐬/dt = -𝐬 × (ξ - 𝐁 + λ 𝐬 × 𝐁),
```
where ``𝐁 = -dE/d𝐬`` is the effective field felt by the expected spin dipole
``𝐬`` and the empirical parameter ``λ`` determines the magnitude of damping.
The components of ``ξ`` are Gaussian white noise, with magnitude ``√(2 k_B T
λ)`` set by a fluctuation-dissipation theorem.

If the `System` has `mode = :SUN`, then this dynamics generalizes [1] to a
stochastic nonlinear Schrödinger equation for SU(_N_) coherent states ``𝐙``,
```math
    d𝐙/dt = -i P [ζ + (1 - i λ̃) ℋ 𝐙].
```
Here, ``P`` projects onto the space orthogonal to ``𝐙``, and ``ζ`` denotes
complex Gaussian white noise with magnitude ``√(2 k_B T λ̃)``. The
local-Hamiltonian ``ℋ`` embeds the energy gradient into the 𝔰𝔲(_N_) Lie
algebra, and generates evolution of spin dipoles, quadrupoles, etc.

When applied to SU(2) coherent states, this generalized dynamics reduces exactly
to the stochastic Landau-Lifshitz equation. The mapping is as follows.
Normalized coherent states ``𝐙`` map to dipole expectation values ``𝐬 = 𝐙^{†}
Ŝ 𝐙``, where spin operators ``Ŝ`` are a spin-``|𝐬|`` representation of
SU(2). The local effective Hamiltonian ``ℋ = -𝐁 ⋅ Ŝ`` generates rotation of
the dipole in analogy to the vector cross product ``S × 𝐁``. The coupling to
the thermal bath maps as ``λ̃ = |𝐬| λ``. Note, however, that the `Langevin`
constructor interprets its `λ` argument as either ``λ`` or ``λ̃``, for modes
`:dipole` or `:SUN`, respectively.

References:

1. [D. Dahlbom et al., Phys. Rev. B 106, 235154 (2022)](https://arxiv.org/abs/2209.01265).
"""
mutable struct Langevin
    Δt  :: Float64
    λ   :: Float64
    kT  :: Float64

    # Averages the square of the dynamical drift term for use in
    # `check_timestep`.
    drift_squared :: OnlineStatistics{Float64}

    function Langevin(Δt; λ, kT)
        Δt <= 0 && error("Select positive Δt")
        return new(Δt, λ, kT, OnlineStatistics{Float64}())
    end
end

"""
    ImplicitMidpoint(Δt::Float64; atol=1e-12) where N

Energy-conserving spin dynamics -- either the Landau-Lifshitz equation, or its
generalization to SU(_N_) coherent states [1]. One call to the [`step!`](@ref)
function will advance a [`System`](@ref) by `Δt` units of time.

Corresponds to the [`Langevin`](@ref) dynamics in the absence of coupling to the
thermal bath (``λ = 0``). Here, however, Sunny uses a more expensive
implicit-midpoint integration scheme that is exactly symplectic [2]. This
approach eliminates energy drift over long simulation trajectories.

References:

1. [H. Zhang and C. D. Batista, Phys. Rev. B 104, 104409 (2021)](https://arxiv.org/abs/2106.14125).
2. [D. Dahlbom et al, Phys. Rev. B 106, 054423 (2022)](https://arxiv.org/abs/2204.07563).
"""
mutable struct ImplicitMidpoint
    Δt   :: Float64
    atol :: Float64

    # Averages the square of the dynamical drift term for use in
    # `check_timestep`.
    drift_squared :: OnlineStatistics{Float64}

    function ImplicitMidpoint(Δt; atol=1e-12)
        Δt <= 0 && error("Select positive Δt")
        return new(Δt, atol, OnlineStatistics{Float64}())
    end    
end


function reset_statistics(integrator::Union{Langevin, ImplicitMidpoint})
    integrator.drift_squared = OnlineStatistics{Float64}()
end

"""
    check_timestep(integrator; tol)

Prints a suggested timestep scale for the relative error tolerance `tol`, based
on statistics collected from previous time integration. A reasonable tolerance
might be of order `tol = 1e-2`, but should be selected on a case-by-case basis.
Note that `tol` is interpreted as the error for the deterministic part of a
single integration timestep, which scales like `Δt²`. Because the stochastic
Heun integration scheme is weakly convergent of order-1, errors in the estimates
of averaged observables may instead scale like ` √tol ~ Δt`. It is the
responsibility of the user to choose `tol` accordingly.

The suggested timestep will be inversely proportional to the largest energy
scale that plays a role in the dynamics. If there is a strong Langvin coupling
``λ`` to the thermal bath, it will effectively rescale the energy scale.

Timestep statistics stored in `integrator` will be reset after calling this
function.
"""
function check_timestep(langevin::Langevin; tol, reset=true)
    (; Δt, λ, kT) = langevin

    drift_magnitude = sqrt(mean(langevin.drift_squared))
    if isnan(drift_magnitude)
        error("Run trajectory first to collect statistics")
    end

    # Heun integration is second-order accurate, so the local error from each
    # deterministic timestep scales as dθ². Angular displacement per timestep dθ
    # scales like dt |drift|, yielding err1 ~ (dt |drift|)^2
    #
    # The "relevant" error introduced by thermal noise should be quantified in
    # terms of the effect on statistical observables. To avoid subtleties in
    # defining this error, we instead naïvely assume it continues to be second
    # order in `dt`. To determine the proportionality constant, consider the
    # high-T limit, where each spin undergoes Brownian motion. Here, the
    # diffusion constant D ~ λ kT sets an inverse time-scale. This implies err2
    # ~ (dt λ kT)².
    #
    # The total error (err1 + err2) should be less than the target tolerance.
    # After some algebra, this implies,
    #
    # dt ≲ sqrt(tol / (c₁² |drift|² + c₂² λ² kT²))
    #
    # for some empirical constants c₁ and c₂.

    c1 = 1.0
    c2 = 1.0
    Δt_bound = sqrt(tol / ((c1*drift_magnitude)^2 + (c2*λ*kT)^2))

    if Δt_bound/2 < Δt < 2Δt_bound
        opinion = "reasonable."
    elseif Δt_bound/5 < Δt < 5Δt_bound
        opinion = "_marginal_."
    elseif Δt <= Δt_bound/5
        opinion = "SMALL!"
    elseif Δt >= 5Δt_bound
        opinion = "LARGE!"
    end

    λkTstr, tolstr, Δtstr, boundstr = number_to_simple_string.((λ*kT, tol, Δt, Δt_bound); digits=4)

    println("Suggest Δt ~ $boundstr for tol = $tolstr.")
    println("Current Δt = $Δtstr is $opinion")
    if c2*λ*kT > c1*drift_magnitude
        println("Thermal noise contributes significantly, with λ*kT = $λkTstr.")
    end

    reset && reset_statistics(langevin)

    return
end


################################################################################
# Dipole integration
################################################################################

@inline rhs_dipole(s, B) = -s × B
@inline rhs_dipole(s, B, λ) = -s × (B + λ * (s × B))

"""
    step!(sys::System, dynamics)

Advance the spin configuration one dynamical time-step. The `dynamics` object
may be a continuous spin dynamics, such as [`Langevin`](@ref) or
[`ImplicitMidpoint`](@ref), or it may be a discrete Monte Carlo sampling scheme
such as [`LocalSampler`](@ref).
"""
function step! end

function step!(sys::System{0}, integrator::Langevin)
    (∇E, s₁, f₁, r₁, ξ) = get_dipole_buffers(sys, 5)
    (; kT, λ, Δt) = integrator
    s = sys.dipoles

    randn!(sys.rng, ξ)
    ξ .*= √(2λ*kT)

    # Euler step
    set_energy_grad_dipoles!(∇E, s, sys)
    @. f₁ = rhs_dipole(s, -∇E, λ)
    @. r₁ = rhs_dipole(s, ξ)   # note absence of λ argument -- noise only appears once in rhs.
    @. s₁ = s + Δt * f₁ + √Δt * r₁

    # Corrector step
    set_energy_grad_dipoles!(∇E, s₁, sys)
    @. s = s + 0.5 * Δt * (f₁ + rhs_dipole(s₁, -∇E, λ)) + 0.5 * √Δt * (r₁ + rhs_dipole(s₁, ξ))
    @. s = normalize_dipole(s, sys.κs)

    # Collect statistics about the magnitude of the drift term squared. For the
    # energy-conserving part, it is important to use |∇E|² rather than |s×∇E|²,
    # because the former gives direct information about the frequency of
    # oscillation, while the latter may artificially vanish (consider, e.g., a
    # spin nearly aligned with an external field).
    for i in eachsite(sys)
        s0 = sys.κs[i] # == norm(s[i])
        ∇E² = ∇E[i]' * ∇E[i]
        accum!(integrator.drift_squared, (1 + (s0*λ)^2) * ∇E²)
    end

    return
end

# The spherical midpoint method, Phys. Rev. E 89, 061301(R) (2014)
# Integrates ds/dt = s × ∂E/∂s one timestep s → s′ via implicit equations
#   s̄ = (s′ + s) / 2
#   ŝ = s̄ / |s̄|
#   (s′ - s)/Δt = 2(s̄ - s)/Δt = - ŝ × B,
# where B = -∂E/∂ŝ.
function step!(sys::System{0}, integrator::ImplicitMidpoint)
    s = sys.dipoles
    (; Δt, atol) = integrator

    (∇E, s̄, ŝ, s̄′) = get_dipole_buffers(sys, 4)
    
    # Initial guess for midpoint
    @. s̄ = s

    max_iters = 100
    for _ in 1:max_iters
        # Integration step for current best guess of midpoint s̄. Produces
        # improved midpoint estimator s̄′.
        @. ŝ = normalize_dipole(s̄, sys.κs)
        set_energy_grad_dipoles!(∇E, ŝ, sys)
        @. s̄′ = s + 0.5 * Δt * rhs_dipole(ŝ, -∇E)

        # If converged, then we can return
        if fast_isapprox(s̄, s̄′,atol=atol* √length(s̄))
            # Normalization here should not be necessary in principle, but it
            # could be useful in practice for finite `atol`.
            @. s = normalize_dipole(2*s̄′ - s, sys.κs)
            return
        end

        @. s̄ = s̄′
    end

    error("Spherical midpoint method failed to converge to tolerance $atol after $max_iters iterations.")
end

function fast_isapprox(x, y; atol)
    acc = 0.
    for i in eachindex(x)
        diff = x[i] - y[i]
        acc += real(dot(diff,diff))
        if acc > atol^2
            return false
        end
    end
    return !isnan(acc)
end


################################################################################
# SU(N) integration
################################################################################

@inline function proj(a::T, Z::T) where T <: CVec
    a - Z * ((Z' * a) / (Z' * Z))
end

function step!(sys::System{N}, integrator::Langevin) where N
    (Z′, ΔZ₁, ΔZ₂, ξ, HZ) = get_coherent_buffers(sys, 5)
    Z = sys.coherents

    randn!(sys.rng, ξ)

    # Prediction
    set_energy_grad_coherents!(HZ, Z, sys)
    rhs_langevin!(ΔZ₁, Z, ξ, HZ, integrator, sys)
    @. Z′ = normalize_ket(Z + ΔZ₁, sys.κs)

    # Correction
    set_energy_grad_coherents!(HZ, Z′, sys)
    rhs_langevin!(ΔZ₂, Z′, ξ, HZ, integrator, sys)
    @. Z = normalize_ket(Z + (ΔZ₁+ΔZ₂)/2, sys.κs)

    # Coordinate dipole data
    @. sys.dipoles = expected_spin(Z)

    # Collect statistics about the magnitude of the drift term squared. Here,
    # |HZ|² plays the role of |∇E|² in the dipole case.
    for i in eachsite(sys)
        HZ² = HZ[i]' * HZ[i]
        accum!(integrator.drift_squared, (1 + integrator.λ^2) * HZ²)
    end

    return
end

function rhs_langevin!(ΔZ::Array{CVec{N}, 4}, Z::Array{CVec{N}, 4}, ξ::Array{CVec{N}, 4},
                       HZ::Array{CVec{N}, 4}, integrator::Langevin, sys::System{N}) where N
    (; kT, λ, Δt) = integrator
    for site in eachsite(sys)
        ΔZ′ = -im*√(2*Δt*kT*λ)*ξ[site] - Δt*(im+λ)*HZ[site]
        ΔZ[site] = proj(ΔZ′, Z[site])
    end
end



# Implicit Midpoint Method applied to the nonlinear Schrödinger dynamics, as
# proposed in Phys. Rev. B 106, 054423 (2022). Integrates dZ/dt = - i H(Z) Z one
# timestep Z → Z′ via the implicit equation
#
#   (Z′-Z)/Δt = - i H(Z̄) Z, where Z̄ = (Z+Z′)/2
#
function step!(sys::System{N}, integrator::ImplicitMidpoint; max_iters=100) where N
    (; atol) = integrator
    (ΔZ, Z̄, Z′, Z″, HZ) = get_coherent_buffers(sys, 5)
    Z = sys.coherents

    @. Z′ = Z 
    @. Z″ = Z 

    for _ in 1:max_iters
        @. Z̄ = (Z + Z′)/2

        set_energy_grad_coherents!(HZ, Z̄, sys)
        rhs_ll!(ΔZ, HZ, integrator, sys)

        @. Z″ = Z + ΔZ

        if fast_isapprox(Z′, Z″, atol=atol*√length(Z′))
            @. Z = normalize_ket(Z″, sys.κs)
            @. sys.dipoles = expected_spin(Z)
            return
        end

        Z′, Z″ = Z″, Z′
    end

    error("Schrödinger midpoint method failed to converge in $max_iters iterations.")
end

function rhs_ll!(ΔZ, HZ, integrator, sys)
    (; Δt) = integrator
    for site in eachsite(sys)
        ΔZ[site] = - Δt*im*HZ[site]
    end
end
