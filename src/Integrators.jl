"""
    Langevin(Δt::Float64; λ::Float64, kT::Float64)

Spin dynamics with damping and noise terms that model coupling to an implicit
thermal bath, of strength `λ`. One call to the [`step!`](@ref) function will
advance a [`System`](@ref) by `Δt` units of time. Can be used to sample from the
Boltzmann distribution at temperature `kT`. An alternative approach to sampling
states from thermal equilibrium is [`LocalSampler`](@ref), which proposes local
Monte Carlo moves. For example, use `LocalSampler` instead of `Langevin` to
sample Ising-like spins.

Setting `λ = 0` disables coupling to the thermal bath, yielding an
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

1. [D. Dahlbom et al., Phys. Rev. B 106, 235154
   (2022)](https://arxiv.org/abs/2209.01265).
"""
mutable struct Langevin
    Δt  :: Float64
    λ   :: Float64
    kT  :: Float64

    function Langevin(Δt; λ, kT)
        Δt <= 0   && error("Select positive Δt")
        kT < 0    && error("Select nonnegative kT")
        λ < 0     && error("Select positive damping λ")
        iszero(λ) && error("Use ImplicitMidpoint instead for energy-conserving dynamics")
        λ < 0.1   && @info "Langevin currently uses Heun integration, which loses statistical accuracy at small λ"
        return new(Δt, λ, kT)
    end
end
Langevin(; λ, kT) = Langevin(NaN; λ, kT)


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

1. [H. Zhang and C. D. Batista, Phys. Rev. B 104, 104409
   (2021)](https://arxiv.org/abs/2106.14125).
2. [D. Dahlbom et al, Phys. Rev. B 106, 054423
   (2022)](https://arxiv.org/abs/2204.07563).
"""
mutable struct ImplicitMidpoint
    Δt   :: Float64
    atol :: Float64

    function ImplicitMidpoint(Δt; atol=1e-12)
        Δt <= 0 && error("Select positive Δt")
        return new(Δt, atol)
    end    
end
ImplicitMidpoint(; atol) = ImplicitMidpoint(NaN; atol)

function check_timestep_available(integrator)
    isnan(integrator.Δt) && error("Set integration timestep `Δt`.")
end

"""
    suggest_timestep(sys, integrator; tol)

Suggests a timestep for the numerical integration of spin dynamics according to
a given error tolerance `tol`. The `integrator` should be [`Langevin`](@ref) or
[`ImplicitMidpoint`](@ref). The suggested ``Δt`` will be inversely proportional
to the magnitude of the effective field ``|dE/d𝐬|`` arising from the current
spin configuration in `sys`. The recommended timestep ``Δt`` scales like `√tol`,
which assumes second-order accuracy of the integrator.

The system `sys` should be initialized to an equilibrium spin configuration for
the target temperature. Alternatively, a reasonably timestep estimate can be
obtained from any low-energy spin configuration. For this, one can use
[`randomize_spins!`](@ref) and then [`minimize_energy!`](@ref).

If `integrator` is of type [`Langevin`](@ref), then the damping magnitude ``λ``
and target temperature ``kT`` will tighten the timestep bound. If ``λ`` exceeds
1, it will rescale the suggested timestep by an approximate the factor ``1/λ``.
If ``kT`` is the largest energy scale, then the suggested timestep will scale
like ``1/λkT``. Quantification of numerical error for stochastic dynamics is
subtle. The stochastic Heun integration scheme is weakly convergent of order-1,
such that errors in the estimates of averaged observables may scale like ``Δt``.
This implies that the `tol` argument may actually scale like the _square_ of the
true numerical error, and should be selected with this in mind.
"""
function suggest_timestep(sys::System{N}, integrator::Langevin; tol) where N
    (; Δt, λ, kT) = integrator
    suggest_timestep_aux(sys; tol, Δt, λ, kT)
end
function suggest_timestep(sys::System{N}, integrator::ImplicitMidpoint; tol) where N
    (; Δt) = integrator
    suggest_timestep_aux(sys; tol, Δt, λ=0, kT=0)
end

function suggest_timestep_aux(sys::System{N}; tol, Δt, λ, kT) where N
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
    Δt_bound = sqrt(tol / ((c1*drift_rms)^2 + (c2*λ*kT)^2))

    # Print suggestion
    bound_str, tol_str = number_to_simple_string.((Δt_bound, tol); digits=4)
    print("Consider Δt ≈ $bound_str for this spin configuration at tol = $tol_str.")

    # Compare with existing Δt if present
    if !isnan(Δt)
        Δt_str = number_to_simple_string(Δt; digits=4)
        if Δt <= Δt_bound/2
            println("\nCurrent value Δt = $Δt_str seems small! Increasing it will make the simulation faster.")
        elseif Δt >= 2Δt_bound
            println("\nCurrent value Δt = $Δt_str seems LARGE! Decreasing it will improve accuracy.")
        else
            println(" Current value is Δt = $Δt_str.")
        end
    end
end


function Base.show(io::IO, integrator::Langevin)
    (; Δt, λ, kT) = integrator
    Δt = isnan(integrator.Δt) ? "<missing>" : repr(Δt)
    println(io, "Langevin($Δt; λ=$λ, kT=$kT)")
end

function Base.show(io::IO, integrator::ImplicitMidpoint)
    (; Δt, atol) = integrator
    Δt = isnan(integrator.Δt) ? "<missing>" : repr(Δt)
    println(io, "ImplicitMidpoint($Δt; atol=$atol)")
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
    check_timestep_available(integrator)

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

    return
end

# The spherical midpoint method, Phys. Rev. E 89, 061301(R) (2014)
# Integrates ds/dt = s × ∂E/∂s one timestep s → s′ via implicit equations
#   s̄ = (s′ + s) / 2
#   ŝ = s̄ / |s̄|
#   (s′ - s)/Δt = 2(s̄ - s)/Δt = - ŝ × B,
# where B = -∂E/∂ŝ.
function step!(sys::System{0}, integrator::ImplicitMidpoint)
    check_timestep_available(integrator)

    s = sys.dipoles
    (; Δt, atol) = integrator
    isnan(Δt) && return

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
    check_timestep_available(integrator)

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
    check_timestep_available(integrator)

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
