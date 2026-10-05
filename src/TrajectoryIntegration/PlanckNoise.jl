################################################################################
# Spectrum fitting
################################################################################
function fd_jacobian(model, x, params; relstep=1e-6)
    y0 = model(x, params)
    nparams = length(params)
    J = zeros(eltype(params), length(x), nparams)

    for n in 1:nparams
        dp = relstep * max(abs(params[n]), 1.0)
        params′ = copy(params)
        params′[n] += dp
        J[:,n] .= (model(x, params′) .- y0) ./ dp
    end

    return J
end

function lsq_fit(model, x, y, p0; maxiter=200, tol=1e-10, λ=1e-6)
    params = copy(p0)
    last_loss = Inf

    for _ in 1:maxiter
        residuals = model(x, params) .- y
        loss = sum(abs2, residuals)

        abs(last_loss - loss) < tol * max(1, loss) && break
        last_loss = loss

        J = fd_jacobian(model, x, params)

        # Levenberg-Marquardt-like damped normal equations
        A = J' * J + λ * I
        b = -J' * residuals

        Δp = A \ b
        trial_params = params + Δp

        # Accept only if improvement; otherwise increase damping
        trial_residuals = model(x, trial_params) .- y
        trial_loss = sum(abs2, trial_residuals)

        if trial_loss < loss
            params = trial_params
            λ /= 2
        else
            λ *= 10
        end

        norm(Δp) < tol * max(1, norm(params)) && break
    end

    return params
end

# Target power spectral density for the Planck noise, |ω| n(|ω|), where n is the
# Bose-Einstein distribution. This is the two-sided convention, ⟨ζ(t) ζ(0)⟩ = ∫
# S(ω) e^{-iωt} dω/2π, defined for all ω and even in ω. It approaches kT as ω →
# 0, so that the classical limit recovers the white noise of `Langevin`. Zero
# point fluctuations are omitted.
function planck_spectrum(ω, kT)
    iszero(ω) && return float(kT)
    iszero(kT) && return zero(float(ω))
    return abs(ω) / expm1(abs(ω)/kT)
end

# Analytical expression for the two-sided power spectrum of the sum of two
# noise-driven second-order linear filters with parameters p.
function filter_spectrum(ω, p)
    c1, c2, Ω1, Ω2, Γ1, Γ2 = p
    @. ((2c1^2 * Γ1) / ((Ω1^2 - ω^2)^2 + ω^2 * Γ1^2)) + ((2c2^2 * Γ2) / ((Ω2^2 - ω^2)^2 + ω^2 * Γ2^2))
end


################################################################################
# Planck noise
################################################################################
mutable struct PlanckNoiseGenerator

    # Basic integrator parameters
    dt      :: Float64
    kT      :: Float64
    damping :: Float64

    # Temperature dependent noise-generation parameters
    Ω₁ :: Float64
    Ω₂ :: Float64
    Γ₁ :: Float64
    Γ₂ :: Float64
    c₁ :: Float64
    c₂ :: Float64

    # State 
    ζ    :: Array{Float64, 5}
    ζbuf :: Array{Float64, 5}
    W1   :: Array{Float64, 5}
    W2   :: Array{Float64, 5}
    u1   :: Array{SVector{2, Float64}, 5}
    u2   :: Array{SVector{2, Float64}, 5}

end

# Solves for when the Planck function reaches 1 percent of its maximum value.
# Used to determine range of ω values for fitting the filter responses. Note
# that this can be solved analytically for an arbitrary percentage, but the
# solution requires the Lambert W function, which is not included in
# SpecialFunctions. The constant used here is specifically for acheiving a value
# that is 1% of kT.
ω_cutoff(kT) = 6.4746008706*kT

# Fit the filter parameters to the Planck spectrum at kT = 1. Because the target
# spectrum has the scaling form kT f(ω/kT), the parameters at any other
# temperature follow exactly from these (see `planck_noise_params`). This
# function is not called at runtime; it documents the origin of
# `planck_noise_params_dimensionless`, and may be used to explore improved fits.
function fit_planck_noise_params()
    ωs = range(0.0, ω_cutoff(1.0), 1000)
    ys = planck_spectrum.(ωs, 1.0)
    p0 = [0.3, 1.88, 1.168055, 2.748380, 3.276618, 5.247509]
    params = lsq_fit(filter_spectrum, ωs, ys, p0)
    all(>(0.0), params) || error("Didn't find good parameters")
    c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂ = params
    return (; c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂)
end

# Output of `fit_planck_noise_params()`. For comparison, the parameters of Savin
# et al. [PRB 86, 064305 (2012)], also used by Barker and Bauer, are c = (0.3429,
# 1.8315), Ω = (1.2223, 2.7189), Γ = (3.2974, 5.0142).
const planck_noise_params_dimensionless = (;
    c₁ = 0.30208998812667376, c₂ = 1.8809221271532812,
    Ω₁ = 1.1750451769311079,  Ω₂ = 2.748367409428579,
    Γ₁ = 3.309977828185631,   Γ₂ = 5.247402673567909,
)

# Filter parameters at temperature kT. A process ζ(t) = kT Φ(kT t), where Φ has
# the dimensionless spectrum f(ω), has spectrum kT f(ω/kT). Rescaling time by
# kT scales frequencies Ω, Γ ∝ kT, and the amplitudes then scale as c ∝ kT².
# At kT = 0 all parameters vanish, and the generated noise is identically zero.
function planck_noise_params(kT)
    (; c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂) = planck_noise_params_dimensionless
    return (; c₁ = c₁*kT^2, c₂ = c₂*kT^2, Ω₁ = Ω₁*kT, Ω₂ = Ω₂*kT, Γ₁ = Γ₁*kT, Γ₂ = Γ₂*kT)
end

function PlanckNoiseGenerator(dt; kT, damping, dims)
    c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂ = planck_noise_params(kT) 
    ζ = zeros(3, dims...)
    ζbuf = zeros(3, dims...)
    W1 = zeros(3, dims...)
    W2 = zeros(3, dims...)
    u1 = zeros(SVector{2, Float64}, 3, dims...)
    u2 = zeros(SVector{2, Float64}, 3, dims...)

    PlanckNoiseGenerator(
        dt, kT, damping,
        Ω₁, Ω₂, Γ₁, Γ₂, c₁, c₂,
        ζ, ζbuf, W1, W2, u1, u2,
    )
end

function set_temperature!(cng::PlanckNoiseGenerator, kT)
    cng.kT = kT
    c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂ = planck_noise_params(kT) 

    cng.c₁ = c₁
    cng.c₂ = c₂
    cng.Ω₁ = Ω₁
    cng.Ω₂ = Ω₂
    cng.Γ₁ = Γ₁
    cng.Γ₂ = Γ₂

    # Reset internal state of noise process
    for i in eachindex(cng.u1)
        cng.u1[i] = zero(SVector{2, Float64})
        cng.u2[i] = zero(SVector{2, Float64})
    end

    return nothing
end

# Reallocate the noise state, if needed, so that there are `ncomp` independent
# noise processes per site of a system with `size(sys.dipoles) == dims`. Dipole
# mode uses ncomp = 3 (a noise field coupling to 𝐒). SU(N) mode uses ncomp = N²
# (a Hermitian noise matrix, see `noise_field_times`).
function ensure_noise_dims!(cng::PlanckNoiseGenerator, dims, ncomp=3)
    size(cng.ζ) == (ncomp, dims...) && return
    cng.ζ = zeros(ncomp, dims...)
    cng.ζbuf = zeros(ncomp, dims...)
    cng.W1 = zeros(ncomp, dims...)
    cng.W2 = zeros(ncomp, dims...)
    cng.u1 = zeros(SVector{2, Float64}, ncomp, dims...)
    cng.u2 = zeros(SVector{2, Float64}, ncomp, dims...)
    return
end

function colored_noise_process_rhs(u, W, dt, Ω, Γ)
    SVector{2, Float64}(
          dt*u[2],
         -dt*(Ω^2*u[1] + Γ*u[2]) + √(2Γ*dt)*W
    )
end

function step_pn!(cng::PlanckNoiseGenerator)
    # Sample gaussians for noise processes to drive filters
    (; W1, W2) = cng
    randn!(W1)
    randn!(W2)

    # Advance filter dynamics
    step_pn_aux!(cng)
end

function step_pn!(rng, cng::PlanckNoiseGenerator)
    # Sample gaussians for noise processes to drive filters
    (; W1, W2) = cng
    randn!(rng, W1)
    randn!(rng, W2)

    # Advance filter dynamics
    step_pn_aux!(cng)
end

function step_pn_aux!(cng::PlanckNoiseGenerator)
    (; W1, W2, ζ, dt, u1, u2, c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂) = cng

    for i in eachindex(ζ)
        # Advance first noise-driven resonant filter
        Δ1 = colored_noise_process_rhs(u1[i], W1[i], dt, Ω₁, Γ₁)
        Δ2 = colored_noise_process_rhs(u1[i] + Δ1, W1[i], dt, Ω₁, Γ₁)
        u1[i] += (Δ1 + Δ2)/2

        # Advance second noise-driven resonant filter
        Δ1 = colored_noise_process_rhs(u2[i], W2[i], dt, Ω₂, Γ₂)
        Δ2 = colored_noise_process_rhs(u2[i] + Δ1, W2[i], dt, Ω₂, Γ₂)
        u2[i] += (Δ1 + Δ2)/2

        # Weighted combination of both filter states
        ζ[i] = c₁*u1[i][1] + c₂*u2[i][1]
    end

    return
end


################################################################################
# Langevin integration with Planck noise 
################################################################################
"""
    LangevinPlanck(dt::Float64; damping::Float64, kT::Float64)

An integrator for Langevin spin dynamics, analogous to [`Langevin`](@ref), in
which the white noise is replaced by colored noise whose power spectrum follows
Planck statistics, ``S(ω) = |ω| n(|ω|)``, where ``n`` is the Bose-Einstein
distribution. A weakly damped mode of frequency ``ω`` then acquires the mean
energy ``ω n(ω)`` of a quantum oscillator (without zero-point energy) rather
than the classical ``k_B T``. At high temperatures the noise becomes white and
the dynamics reduces to that of `Langevin`. The colored noise is generated by a
pair of noise-driven second-order filters, following Refs. [1, 2].

In `:SUN` mode, the noise enters as a random Hermitian local Hamiltonian that
couples to all ``N^2-1`` generators of SU(_N_). This is the natural analogue of
the dipole noise field. In the linear regime, each generalized spin-wave mode
then acquires energy ``ω n(ω)``, as in `:dipole` mode. (Replacing the complex
white noise of [`Langevin`](@ref) by colored noise is not equivalent, and would
overheat modes with ``ω ≫ k_B T``.)

The colored noise is combined with memoryless damping. Mode populations are
therefore accurate only in the limit of weak damping. Note that in `:dipole`
mode, the relaxation rate of a mode of frequency ``ω`` is ``|𝐒| λ ω``, where
``λ`` is the `damping` parameter. At finite damping, modes with ``ω ≫ k_B T``
acquire an excess energy that scales linearly with `damping`.

To accurately integrate the noise filters, the timestep should satisfy
`dt ≤ 0.1/kT`. A warning is printed otherwise.

## References

1. [A. V. Savin, Y. A. Kosevich, and A. Cantarero, _Semiquantum molecular
   dynamics simulation of thermal properties and heat transport in
   low-dimensional nanostructures_, Phys. Rev. B **86**, 064305
   (2012)](https://doi.org/10.1103/PhysRevB.86.064305).
2. [J. Barker and G. E. W. Bauer, _Semiquantum thermodynamics of complex
   ferrimagnets_, Phys. Rev. B **100**, 140401(R)
   (2019)](https://doi.org/10.1103/PhysRevB.100.140401).
"""
mutable struct LangevinPlanck <: AbstractIntegrator
    dt              :: Float64
    damping         :: Float64
    kT              :: Float64
    noisesource     :: PlanckNoiseGenerator

    # The noise state is allocated on first use, sized to the system.
    function LangevinPlanck(sys, dt=NaN; λ=nothing, damping=nothing, kT)
        if !isnothing(λ)
            @warn "`λ` argument is deprecated! Use `damping` instead."
            damping = @something damping λ
        end
        isnothing(damping) && error("`damping` parameter required")
        iszero(damping) && error("Use ImplicitMidpoint instead for energy-conserving dynamics")

        dt <= 0         && error("Select positive dt")
        kT < 0          && error("Select nonnegative kT")
        damping <= 0    && error("Select positive damping")

        cng = PlanckNoiseGenerator(dt; kT, damping, dims=(0, 0, 0, 0))
        check_noise_timestep(dt, kT)
        return new(dt, damping, kT, cng)
    end
end

# The copy has fresh noise state, allocated on first use.
function Base.copy(dyn::LangevinPlanck)
    LangevinPlanck(dyn.dt; dyn.damping, dyn.kT)
end

# The noise filters are integrated with the same timestep as the spins. Their
# damping rates reach Γ₂ ≈ 5.2 kT, and `dt ≤ 0.1/kT` reproduces the analytical
# filter spectrum to within statistical error.
function check_noise_timestep(dt, kT)
    if !isnan(dt) && dt > (1 + 1e-12) * 0.1/kT
        @warn "LangevinPlanck timestep dt = $dt exceeds 0.1/kT = $(0.1/kT). The Planck noise filters may be inaccurate." maxlog=1
    end
end

function Base.setproperty!(integrator::LangevinPlanck, sym::Symbol, val)
    cng = integrator.noisesource
    if sym == :dt
        setfield!(integrator, sym, convert(Float64, val))
        cng.dt = val
        check_noise_timestep(integrator.dt, integrator.kT)
    elseif sym == :kT
        val < 0 && error("Select nonnegative kT")
        setfield!(integrator, sym, convert(Float64, val))
        set_temperature!(cng, val)
        check_noise_timestep(integrator.dt, integrator.kT)
    elseif sym == :damping
        setfield!(integrator, sym, convert(Float64, val))
        cng.damping = val
    else
        setfield!(integrator, sym, val)
    end
end

# See `check_noise_timestep`.
noise_timestep_bound(integrator::LangevinPlanck) = iszero(integrator.kT) ? Inf : 0.1/integrator.kT

function Base.show(io::IO, integrator::LangevinPlanck)
    (; dt, damping, kT) = integrator
    dt = isnan(integrator.dt) ? "<missing>" : repr(dt)
    println(io, "LangevinPlanck($dt; damping=$damping, kT=$kT)")
end

@inline function rhs_dipole_pn!(ΔS, S, ξ, ∇E, integrator)
    (; dt, damping) = integrator
    λ = damping

    @. ΔS = - S × (ξ + dt*∇E - dt*λ*(S × ∇E))
end

@inline function advance_and_retrieve_noise!(sys, integrator)
    (; damping, noisesource, dt) = integrator
    cng = noisesource
    ensure_noise_dims!(cng, size(sys.dipoles))
    step_pn!(sys.rng, cng)
    ζ = view(reinterpret(SVector{3, Float64}, cng.ζ), 1, :, :, :, :)
    ζ .*= sqrt(2damping)*dt # Note dt here -- treat as noise field
    return ζ
end

function step!(sys::System{0}, integrator::LangevinPlanck)
    check_timestep_available(integrator)

    (S′, ΔS₁, ΔS₂, ∇E) = get_dipole_buffers(sys, 5)
    S = sys.dipoles

    ζ = advance_and_retrieve_noise!(sys, integrator)

    # Euler prediction step
    set_energy_grad_dipoles!(∇E, S, sys)
    rhs_dipole_pn!(ΔS₁, S, ζ, ∇E, integrator)
    @. S′ = normalize_dipole(S + ΔS₁, sys.κs)

    # Correction step
    set_energy_grad_dipoles!(∇E, S′, sys)
    rhs_dipole_pn!(ΔS₂, S′, ζ, ∇E, integrator)
    @. S = normalize_dipole(S + (ΔS₁+ΔS₂)/2, sys.κs)

    return
end

# For SU(N) coherent states, the noise enters as a random local Hamiltonian,
#
#     dZ/dt = -i P [X(t) + (1 - iλ) ℋ] Z,
#
# with X Hermitian. This "field form" is the direct analogue of the dipole
# noise field, which couples to 𝐒. It generalizes Eq. (42) of Dahlbom et al.,
# PRB 106, 235154 (2022). For white noise it is equivalent to the complex noise
# vector ζ used by `Langevin` (their Eq. 46), but for colored noise it is not.
# The vector form contains multiplicative amplitude noise that samples the
# spectrum at zero frequency, S(0) = kT, which overheats modes with ω ≫ kT. In
# the field form, the corresponding term is pure phase noise.
#
# X is built from N² independent Planck processes per site. The real and
# imaginary parts of each off-diagonal element have spectrum λ S(ω), and each
# diagonal element has spectrum 2λ S(ω). This is a colored Gaussian unitary
# ensemble. It is unitarily invariant and, for S = kT, reproduces the
# normalization ⟨ζ_a* ζ_b⟩ = 2 λ kT δ_ab of `Langevin`. In the linear regime,
# each mode of frequency ω then acquires energy S(ω), with no N-dependent
# factor.
#
# Returns dt⋅X⋅Z for the noise matrix of `site`, where ζ[k, site] holds the k-th
# unit-normalized noise process.
@inline function noise_field_times(ζ, site, Z::CVec{N}, λ, dt) where N
    cdiag = dt*sqrt(2λ)
    coff = dt*sqrt(λ)
    return CVec{N}(ntuple(N) do a
        # Components are ordered: for each a, the diagonal element (a, a),
        # then the real and imaginary parts of (a, b) for b > a.
        acc = zero(ComplexF64)
        for b in 1:N
            acc += noise_matrix_element(ζ, site, a, b, N, cdiag, coff) * Z[b]
        end
        acc
    end)
end

# Linear index into the N² processes for the upper-triangular storage order
# described above.
@inline function noise_offset(a, N)
    # Number of processes used by rows 1..a-1: each row r uses 1 + 2(N-r)
    return (a-1) + (a-1)*(2N - a)
end

@inline function noise_matrix_element(ζ, site, a, b, N, cdiag, coff)
    if a == b
        return complex(cdiag * ζ[noise_offset(a, N) + 1, site])
    elseif a < b
        k = noise_offset(a, N) + 1 + 2(b - a) - 1
        return coff * complex(ζ[k, site], ζ[k+1, site])
    else
        k = noise_offset(b, N) + 1 + 2(a - b) - 1
        return coff * complex(ζ[k, site], -ζ[k+1, site])
    end
end

function step!(sys::System{N}, integrator::LangevinPlanck) where N
    check_timestep_available(integrator)

    (; damping, dt, noisesource) = integrator
    (Z′, ΔZ₁, ΔZ₂, ξ, HZ) = get_coherent_buffers(sys, 5)
    Z = sys.coherents

    ensure_noise_dims!(noisesource, size(Z), N^2)
    step_pn!(sys.rng, noisesource)
    ζ = reshape(noisesource.ζ, N^2, :)

    # Euler prediction step. The noise term of `rhs_sun!` is -P ξ, so pass
    # ξ = i dt X Z.
    for i in eachindex(Z)
        ξ[i] = im * noise_field_times(ζ, i, Z[i], damping, dt)
    end
    set_energy_grad_coherents!(HZ, Z, sys)
    rhs_sun!(ΔZ₁, Z, ξ, HZ, integrator)
    @. Z′ = normalize_ket(Z + ΔZ₁, sys.κs)

    # Correction step. The multiplicative noise is re-evaluated at Z′, with the
    # same noise matrix (Stratonovich-consistent Heun).
    for i in eachindex(Z)
        ξ[i] = im * noise_field_times(ζ, i, Z′[i], damping, dt)
    end
    set_energy_grad_coherents!(HZ, Z′, sys)
    rhs_sun!(ΔZ₂, Z′, ξ, HZ, integrator)
    @. Z = normalize_ket(Z + (ΔZ₁+ΔZ₂)/2, sys.κs)

    # Coordinate dipole data
    sync_dipoles!(sys)

    return
end
