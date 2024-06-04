"""
    Langevin(dt::Float64; damping::Float64, kT::Float64)

An integrator for Langevin spin dynamics using the explicit Heun method. The
`damping` parameter controls the coupling to an implicit thermal bath. One call
to the [`step!`](@ref) function will advance a [`System`](@ref) by `dt` units of
time. Can be used to sample from the Boltzmann distribution at temperature `kT`.
An alternative approach to sampling states from thermal equilibrium is
[`LocalSampler`](@ref), which proposes local Monte Carlo moves. For example, use
`LocalSampler` instead of `Langevin` to sample Ising-like spins.

Setting `damping = 0` disables coupling to the thermal bath, yielding an
energy-conserving spin dynamics. The `Langevin` integrator uses an explicit
numerical integrator which cannot prevent energy drift. Alternatively, the
[`ImplicitMidpoint`](@ref) method can be used, which is more expensive but
prevents energy drift through exact conservation of the symplectic 2-form.

If the [`System`](@ref) has `mode = :dipole`, then the dynamics is the
stochastic Landau-Lifshitz equation,
```math
    d𝐬/dt = -𝐬 × (ξ - 𝐁 + λ 𝐬 × 𝐁),
```
where ``𝐁 = -dE/d𝐬`` is the effective field felt by the expected spin dipole
``𝐬``. The components of ``ξ`` are Gaussian white noise, with magnitude ``√(2
k_B T λ)`` set by a fluctuation-dissipation theorem. The parameter `damping`
sets the phenomenological coupling ``λ`` to the thermal bath. 

If the `System` has `mode = :SUN`, then this dynamics generalizes [1] to a
stochastic nonlinear Schrödinger equation for SU(_N_) coherent states ``𝐙``,
```math
    d𝐙/dt = -i P [ζ + (1 - i λ̃) ℋ 𝐙].
```
Here, ``P`` projects onto the space orthogonal to ``𝐙``, and ``ζ`` denotes
complex Gaussian white noise with magnitude ``√(2 k_B T λ̃)``. The
local-Hamiltonian ``ℋ`` embeds the energy gradient into the 𝔰𝔲(_N_) Lie
algebra, and generates evolution of spin dipoles, quadrupoles, etc. The
parameter `damping` here sets ``λ̃``, which is analogous to ``λ`` above.

When applied to SU(2) coherent states, the generalized spin dynamics reduces
exactly to the stochastic Landau-Lifshitz equation. The mapping is as follows.
Normalized coherent states ``𝐙`` map to dipole expectation values ``𝐬 = 𝐙^{†}
Ŝ 𝐙``, where spin operators ``Ŝ`` are a spin-``|𝐬|`` representation of
SU(2). The local effective Hamiltonian ``ℋ = -𝐁 ⋅ Ŝ`` generates rotation of
the dipole in analogy to the vector cross product ``S × 𝐁``. The coupling to
the thermal bath maps as ``λ̃ = |𝐬| λ``. Note, therefore, that the scaling of
the `damping` parameter varies subtly between `:dipole` and `:SUN` modes.

References:

1. [D. Dahlbom et al., Phys. Rev. B 106, 235154
   (2022)](https://arxiv.org/abs/2209.01265).
"""
mutable struct Langevin
    dt      :: Float64
    damping :: Float64
    kT      :: Float64

    function Langevin(dt; λ=nothing, damping=nothing, kT)
        if !isnothing(λ)
            @warn "`λ` argument is deprecated! Use `damping` instead."
            damping = @something damping λ
        end
        isnothing(damping) && error("`damping` parameter required")
        iszero(damping) && error("Use ImplicitMidpoint instead for energy-conserving dynamics")

        dt <= 0         && error("Select positive dt")
        kT < 0          && error("Select nonnegative kT")
        damping <= 0    && error("Select positive damping")
        return new(dt, damping, kT)
    end
end

function Langevin(; λ=nothing, damping=nothing, kT)
    Langevin(NaN; λ, damping, kT)
end

function Base.copy(dyn::Langevin)
    Langevin(dyn.dt; dyn.damping, dyn.kT)
end

#=
Damping and noise terms may be included through the optional `damping` and `kT`
parameters. In this case, the spin dynamics will coincide with that of
[`Langevin`](@ref), and samples the classical Boltzmann distribution [2].
Relative to the Heun integration method, the implicit midpoint method has a
larger numerical cost, but can achieve much better statistical accuracy,
especially in the limit of small `damping`.

2. [D. Dahlbom et al, Phys. Rev. B 106, 054423
   (2022)](https://arxiv.org/abs/2204.07563).
=#

"""
    ImplicitMidpoint(dt::Float64; atol=1e-12) where N

The implicit midpoint method for integrating the Landau-Lifshitz spin dynamics
or its generalization to SU(_N_) coherent states [1]. One call to the
[`step!`](@ref) function will advance a [`System`](@ref) by `dt` units of time.
This integration scheme is exactly symplectic and eliminates energy drift over
arbitrarily long simulation trajectories.

References:

1. [H. Zhang and C. D. Batista, Phys. Rev. B 104, 104409
   (2021)](https://arxiv.org/abs/2106.14125).
"""
mutable struct ImplicitMidpoint
    dt      :: Float64
    damping :: Float64
    kT      :: Float64
    atol    :: Float64

    function ImplicitMidpoint(dt; damping=0, kT=0, atol=1e-12)
        dt <= 0      && error("Select positive dt")
        kT < 0       && error("Select nonnegative kT")
        damping < 0  && error("Select nonnegative damping")

        # Noise in the implicit midpoint method can be problematic, because rare
        # events can lead to very slow convergence of the fixed point
        # iterations. Perhaps it could be made to work if we clip the sampled
        # noise to a restricted magnitude? For now, simply disable the feature.
        iszero(kT) || error("ImplicitMidpoint with a Langevin thermostat is not currently supported.")

        return new(dt, damping, kT, atol)
    end    
end
ImplicitMidpoint(; atol) = ImplicitMidpoint(NaN; atol)


function check_timestep_available(integrator)
    isnan(integrator.dt) && error("Set integration timestep `dt`.")
end

"""
    suggest_timestep(sys, integrator; tol)

Suggests a timestep for the numerical integration of spin dynamics according to
a given error tolerance `tol`. The `integrator` should be [`Langevin`](@ref) or
[`ImplicitMidpoint`](@ref). The suggested ``dt`` will be inversely proportional
to the magnitude of the effective field ``|dE/d𝐬|`` arising from the current
spin configuration in `sys`. The recommended timestep ``dt`` scales like `√tol`,
which assumes second-order accuracy of the integrator.

The system `sys` should be initialized to an equilibrium spin configuration for
the target temperature. Alternatively, a reasonably timestep estimate can be
obtained from any low-energy spin configuration. For this, one can use
[`randomize_spins!`](@ref) and then [`minimize_energy!`](@ref).

Large `damping` magnitude or target temperature `kT` will tighten the timestep
bound. If `damping` exceeds 1, it will rescale the suggested timestep by an
approximate the factor ``1/damping``. If `kT` is the largest energy scale, then
the suggested timestep will scale like `1/(damping*kT)`. Quantification of
numerical error for stochastic dynamics is subtle. The stochastic Heun
integration scheme is weakly convergent of order-1, such that errors in the
estimates of averaged observables may scale like `dt`. This implies that the
`tol` argument may actually scale like the _square_ of the true numerical error,
and should be selected with this in mind.
"""
function suggest_timestep(sys::System{N}, integrator::Union{Langevin, ImplicitMidpoint}; tol) where N
    (; dt) = integrator
    dt_bound = suggest_timestep_aux(sys, integrator; tol)

    # Print suggestion
    bound_str, tol_str = number_to_simple_string.((dt_bound, tol); digits=4)
    print("Consider dt ≈ $bound_str for this spin configuration at tol = $tol_str.")

    # Compare with existing dt if present
    if !isnan(dt)
        dt_str = number_to_simple_string(dt; digits=4)
        if dt <= dt_bound/2
            println("\nCurrent value dt = $dt_str seems small! Increasing it will make the simulation faster.")
        elseif dt >= 2dt_bound
            println("\nCurrent value dt = $dt_str seems LARGE! Decreasing it will improve accuracy.")
        else
            println(" Current value is dt = $dt_str.")
        end
    else
        println()
    end
end

function suggest_timestep_aux(sys::System{N}, integrator; tol) where N
    (; damping, kT) = integrator
    λ = damping

    # Accumulate statistics regarding Var[∇E]
    acc = 0.0
    if N == 0
        ∇Es, = get_dipole_buffers(sys, 1)
        set_energy_grad_dipoles!(∇Es, sys.dipoles, sys)
        for (κ, ∇E) in zip(sys.κs, ∇Es)
            # In dipole mode, the spin magnitude `κ = |s|` scales the effective
            # damping rate.
            acc += (1 + (κ*λ)^2) * norm(∇E)^2
        end
    else
        ∇Es, = get_coherent_buffers(sys, 1)
        set_energy_grad_coherents!(∇Es, sys.coherents, sys)
        for ∇E in ∇Es
            acc += (1 + λ^2) * norm(∇E)^2
        end
    end

    # `drift_rms` gives the root-mean-squared of the drift term for one
    # integration timestep of the Langevin dynamics. It is associated with the
    # angular velocity dθ/dt where dθ ~ dS/|S| or dZ/|Z| for :dipole or :SUN
    # mode, respectively. In calculating `drift_rms`, it is important to use the
    # energy gradient |∇E| directly, rather than projecting out the component of
    # ∇E aligned with the spin. Without projection, one obtains direct
    # information about the frequency of oscillation. Consider, e.g., a spin
    # approximately aligned with an external field: the precession frequency is
    # given by |∇E| = |B|.
    drift_rms = sqrt(acc/length(eachsite(sys)))
    if iszero(drift_rms)
        error("Cannot suggest a timestep without an energy scale!")
    end

    # In a second-order integrator, the local error from each deterministic
    # timestep scales as dθ². Angular displacement per timestep dθ scales like
    # dt drift_rms, yielding err1 ~ (dt drift_rms)^2
    #
    # Quantifying the "error" introduced by thermal noise is subtle. E.g., for
    # weak convergence, we should consider the effect on statistical
    # observables. We avoid all subtleties by naïvely assuming this error
    # continues to be second order in `dt`. To determine the proportionality
    # constant, consider the high-T limit, where each spin undergoes Brownian
    # motion. Here, the diffusion constant D ~ λ kT sets an inverse time-scale.
    # This implies err2 ~ (dt λ kT)².
    #
    # The total error (err1 + err2) should be less than the target tolerance.
    # After some algebra, this implies,
    #
    # dt ≲ sqrt(tol / (c₁² drift_rms² + c₂² λ² kT²))
    #
    # for some empirical constants c₁ and c₂.
    c1 = 1.0
    c2 = 1.0
    dt_bound = sqrt(tol / ((c1*drift_rms)^2 + (c2*λ*kT)^2))
    return dt_bound
end


function Base.show(io::IO, integrator::Langevin)
    (; dt, damping, kT) = integrator
    dt = isnan(integrator.dt) ? "<missing>" : repr(dt)
    println(io, "Langevin($dt; damping=$damping, kT=$kT)")
end

function Base.show(io::IO, integrator::ImplicitMidpoint)
    (; dt, atol) = integrator
    dt = isnan(integrator.dt) ? "<missing>" : repr(dt)
    println(io, "ImplicitMidpoint($dt; atol=$atol)")
end


################################################################################
# Dipole integration
################################################################################


@inline function rhs_dipole!(Δs, s, ξ, ∇E, integrator)
    (; dt, damping) = integrator
    λ = damping

    if iszero(λ)
        @. Δs = - s × (dt*∇E)
    else
        @. Δs = - s × (ξ + dt*∇E - dt*λ*(s × ∇E))
    end
end

function rhs_sun!(ΔZ, Z, ζ, HZ, integrator)
    (; damping, dt) = integrator
    λ = damping

    if iszero(λ)
        @. ΔZ = - im*dt*HZ
    else
        @. ΔZ = - proj(ζ + dt*(im+λ)*HZ, Z)
    end
end

function fill_noise!(rng, ξ, integrator)
    (; dt, damping, kT) = integrator
    λ = damping

    if iszero(λ) || iszero(kT)
        fill!(ξ, zero(eltype(ξ)))
    else
        randn!(rng, ξ)
        ξ .*= √(2dt*λ*kT)
    end
end


"""
    step!(sys::System, dynamics)

Advance the spin configuration one dynamical time-step. The `dynamics` object
may be a continuous spin dynamics, such as [`Langevin`](@ref) or
[`ImplicitMidpoint`](@ref), or it may be a discrete Monte Carlo sampling scheme
such as [`LocalSampler`](@ref).
"""
function step! end

# Heun integration with normalization

function step!(sys::System{0}, integrator::Langevin)
    check_timestep_available(integrator)

    (s′, Δs₁, Δs₂, ξ, ∇E) = get_dipole_buffers(sys, 5)
    s = sys.dipoles

    fill_noise!(sys.rng, ξ, integrator)

    # Euler prediction step
    set_energy_grad_dipoles!(∇E, s, sys)
    rhs_dipole!(Δs₁, s, ξ, ∇E, integrator)
    @. s′ = normalize_dipole(s + Δs₁, sys.κs)

    # Correction step
    set_energy_grad_dipoles!(∇E, s′, sys)
    rhs_dipole!(Δs₂, s′, ξ, ∇E, integrator)
    @. s = normalize_dipole(s + (Δs₁+Δs₂)/2, sys.κs)

    return
end

function step!(sys::System{N}, integrator::Langevin) where N
    check_timestep_available(integrator)

    (Z′, ΔZ₁, ΔZ₂, ζ, HZ) = get_coherent_buffers(sys, 5)
    Z = sys.coherents

    fill_noise!(sys.rng, ζ, integrator)

    # Euler prediction step
    set_energy_grad_coherents!(HZ, Z, sys)
    rhs_sun!(ΔZ₁, Z, ζ, HZ, integrator)
    @. Z′ = normalize_ket(Z + ΔZ₁, sys.κs)

    # Correction step
    set_energy_grad_coherents!(HZ, Z′, sys)
    rhs_sun!(ΔZ₂, Z′, ζ, HZ, integrator)
    @. Z = normalize_ket(Z + (ΔZ₁+ΔZ₂)/2, sys.κs)

    # Coordinate dipole data
    @. sys.dipoles = expected_spin(Z)

    return
end


# Variants of the implicit midpoint method

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

# The spherical midpoint method, Phys. Rev. E 89, 061301(R) (2014)
# Integrates ds/dt = s × ∂E/∂s one timestep s → s′ via implicit equations
#   s̄ = (s′ + s) / 2
#   ŝ = s̄ / |s̄|
#   (s′ - s)/dt = 2(s̄ - s)/dt = - ŝ × B,
# where B = -∂E/∂ŝ.
function step!(sys::System{0}, integrator::ImplicitMidpoint; max_iters=100)
    check_timestep_available(integrator)

    s = sys.dipoles
    atol = integrator.atol * √length(s)

    (Δs, ŝ, s′, s″, ξ, ∇E) = get_dipole_buffers(sys, 6)

    fill_noise!(sys.rng, ξ, integrator)
    
    @. s′ = s
    @. s″ = s

    for _ in 1:max_iters
        # Current guess for midpoint ŝ
        @. ŝ = normalize_dipole((s + s′)/2, sys.κs)

        set_energy_grad_dipoles!(∇E, ŝ, sys)
        rhs_dipole!(Δs, ŝ, ξ, ∇E, integrator)

        @. s″ = s + Δs

        # If converged, then we can return
        if fast_isapprox(s′, s″; atol)
            # Normalization here should not be necessary in principle, but it
            # could be useful in practice for finite `atol`.
            @. s = normalize_dipole(s″, sys.κs)
            return
        end

        s′, s″ = s″, s′
    end

    error("Spherical midpoint method failed to converge to tolerance $atol after $max_iters iterations.")
end


# Implicit Midpoint Method applied to the nonlinear Schrödinger dynamics, as
# proposed in Phys. Rev. B 106, 054423 (2022). Integrates dZ/dt = - i H(Z) Z one
# timestep Z → Z′ via the implicit equation
#
#   (Z′-Z)/dt = - i H(Z̄) Z, where Z̄ = (Z+Z′)/2
#
function step!(sys::System{N}, integrator::ImplicitMidpoint; max_iters=100) where N
    check_timestep_available(integrator)

    Z = sys.coherents
    atol = integrator.atol * √length(Z)
    
    (ΔZ, Z̄, Z′, Z″, ζ, HZ) = get_coherent_buffers(sys, 6)
    fill_noise!(sys.rng, ζ, integrator)

    @. Z′ = Z 
    @. Z″ = Z 

    for _ in 1:max_iters
        @. Z̄ = (Z + Z′)/2

        set_energy_grad_coherents!(HZ, Z̄, sys)
        rhs_sun!(ΔZ, Z̄, ζ, HZ, integrator)

        @. Z″ = Z + ΔZ

        if fast_isapprox(Z′, Z″; atol)
            @. Z = normalize_ket(Z″, sys.κs)
            @. sys.dipoles = expected_spin(Z)
            return
        end

        Z′, Z″ = Z″, Z′
    end

    error("Schrödinger midpoint method failed to converge in $max_iters iterations.")
end
