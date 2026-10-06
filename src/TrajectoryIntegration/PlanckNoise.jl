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
    dt :: Float64

    # Temperature, uniform or per site. The filter parameters follow from
    # `planck_noise_params_dimensionless` by scaling with the (local) kT.
    kT :: Union{Float64, Array{Float64, 4}}

    # Noise output and filter state, (ncomp, dims...)
    ζ    :: Array{Float64, 5}
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

# Allocates `ncomp` independent noise processes for each site of a system with
# `size(sys.dipoles) == dims`. Dipole mode uses ncomp = 3 (a noise field
# coupling to 𝐒). SU(N) mode uses ncomp = N² (a Hermitian noise matrix, see
# `noise_field_times`). The temperature `kT` may be a scalar, or an array of
# size `dims`.
function PlanckNoiseGenerator(dt; kT, dims, ncomp=3)
    ζ = zeros(ncomp, dims...)
    W1 = zeros(ncomp, dims...)
    W2 = zeros(ncomp, dims...)
    u1 = zeros(SVector{2, Float64}, ncomp, dims...)
    u2 = zeros(SVector{2, Float64}, ncomp, dims...)
    return PlanckNoiseGenerator(dt, scalar_or_site_array(kT, dims), ζ, W1, W2, u1, u2)
end

# Changes the temperature of every site, preserving the stationarity of the
# noise. Each filter obeys ü + Γ u̇ + Ω² u = √(2Γ) W, with stationary
# distribution ∝ exp[-(u̇² + Ω² u²)/2], so ⟨u²⟩ = 1/Ω² ∝ 1/kT² and ⟨u̇²⟩ = 1.
# Rescaling u by kT_old/kT_new, with u̇ unchanged, maps the stationary state at
# kT_old exactly onto the stationary state at kT_new. Sites with kT_new = 0 are
# reset to zero. Sites heated from kT_old = 0 start from zero, and become
# stationary after a time of order 1/kT_new.
function set_temperature!(cng::PlanckNoiseGenerator, kT)
    (; u1, u2) = cng
    ncomp = size(u1, 1)
    dims = size(u1)[2:end]
    kT_new = scalar_or_site_array(kT, dims)
    for site in 1:prod(dims)
        T₀, T₁ = site_value(cng.kT, site), site_value(kT_new, site)
        T₀ == T₁ && continue
        for k in 1:ncomp
            i = k + (site-1)*ncomp
            if iszero(T₁)
                u1[i] = u2[i] = zero(SVector{2, Float64})
            elseif !iszero(T₀)
                r = T₀ / T₁
                u1[i] = SVector(r*u1[i][1], u1[i][2])
                u2[i] = SVector(r*u2[i][1], u2[i][2])
            end
        end
    end
    cng.kT = kT_new
    return nothing
end

# The noise state is allocated for a specific system. Using the integrator with
# a system of different size or mode is an error.
function check_noise_dims(cng::PlanckNoiseGenerator, dims, ncomp)
    if size(cng.ζ) != (ncomp, dims...)
        error("LangevinPlanck was constructed for a different system. Construct a new integrator for this system.")
    end
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

# Advances each filter by one Heun step. A function barrier selects the uniform
# or site-dependent temperature path.
step_pn_aux!(cng::PlanckNoiseGenerator) = step_pn_aux!(cng, cng.kT)

@inline function advance_filters!((; W1, W2, ζ, dt, u1, u2), i, c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂)
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

function step_pn_aux!(cng::PlanckNoiseGenerator, kT::Float64)
    (; c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂) = planck_noise_params(kT)
    # Load the arrays once; fields of a mutable struct are reloaded otherwise
    arrays = (; cng.W1, cng.W2, cng.ζ, cng.dt, cng.u1, cng.u2)
    for i in eachindex(cng.ζ)
        advance_filters!(arrays, i, c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂)
    end
    return
end

function step_pn_aux!(cng::PlanckNoiseGenerator, kT::Array{Float64, 4})
    ncomp = size(cng.ζ, 1)
    arrays = (; cng.W1, cng.W2, cng.ζ, cng.dt, cng.u1, cng.u2)
    for site in eachindex(kT)
        (; c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂) = planck_noise_params(kT[site])
        for k in 1:ncomp
            advance_filters!(arrays, k + (site-1)*ncomp, c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂)
        end
    end
    return
end

################################################################################
# Langevin integration with Planck noise 
################################################################################
"""
    LangevinPlanck(sys::System, dt::Float64; damping, kT)

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
`dt ≤ 0.1/kT`. A warning is printed otherwise. If `dt` is omitted, it must be
set before integration, e.g., using [`suggest_timestep`](@ref).

The colored noise has internal state for every site of `sys`, which is
allocated by the constructor. The integrator can only be used with `sys`, or
with another system of the same size and mode.

Each of `damping` and `kT` may be a number, or an array of size
`size(sys.dipoles)` to make it site dependent, e.g., for a temperature gradient.
Each site then receives independent noise with the local Planck spectrum. Damping
may vanish on some sites, which then exchange energy with the bath only through
their neighbors. Assigning a new value, `integrator.kT = kT′`, rescales the
internal noise state so that the noise remains stationary at the new
temperatures. An array-valued `integrator.kT` may also be modified in place, in
which case the noise adapts to the new temperatures over a time of order
``1/k_B T``.

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
    damping         :: Union{Float64, Array{Float64, 4}}
    kT              :: Union{Float64, Array{Float64, 4}}
    noisesource     :: PlanckNoiseGenerator

    function LangevinPlanck(dt, damping, kT, noisesource::PlanckNoiseGenerator)
        return new(dt, damping, kT, noisesource)
    end
end

function LangevinPlanck(sys::System{N}, dt=NaN; λ=nothing, damping=nothing, kT) where N
    if !isnothing(λ)
        @warn "`λ` argument is deprecated! Use `damping` instead."
        damping = @something damping λ
    end
    isnothing(damping) && error("`damping` parameter required")
    dt <= 0 && error("Select positive dt")

    dims = size(sys.dipoles)
    damping = validate_damping(damping, dims)
    kT = validate_kT(kT, dims)

    ncomp = N == 0 ? 3 : N^2
    cng = PlanckNoiseGenerator(dt; kT, dims, ncomp)
    check_noise_timestep(dt, kT)
    # An array-valued kT is the noise source's own array, so that in-place
    # modifications take effect.
    return LangevinPlanck(Float64(dt), damping, kT isa Real ? kT : cng.kT, cng)
end

# The copy has fresh noise state, of the same size.
function Base.copy(dyn::LangevinPlanck)
    (; dt, damping, kT) = dyn
    sz = size(dyn.noisesource.ζ)
    cng = PlanckNoiseGenerator(dt; kT, dims=sz[2:end], ncomp=sz[1])
    return LangevinPlanck(dt, copy(damping), kT isa Real ? kT : cng.kT, cng)
end

# The noise filters are integrated with the same timestep as the spins. Their
# damping rates reach Γ₂ ≈ 5.2 kT, and `dt ≤ 0.1/kT` reproduces the analytical
# filter spectrum to within statistical error. With site-dependent kT, the
# bound is set by the hottest site.
function check_noise_timestep(dt, kT)
    kTmax = maximum(kT)
    if !isnan(dt) && dt > (1 + 1e-12) * 0.1/kTmax
        @warn "LangevinPlanck timestep dt = $dt exceeds 0.1/kT = $(0.1/kTmax). The Planck noise filters may be inaccurate." maxlog=1
    end
end

function Base.setproperty!(integrator::LangevinPlanck, sym::Symbol, val)
    cng = integrator.noisesource
    dims = size(cng.ζ)[2:end]
    if sym == :dt
        setfield!(integrator, sym, convert(Float64, val))
        cng.dt = val
        check_noise_timestep(integrator.dt, integrator.kT)
    elseif sym == :kT
        kT = validate_kT(val, dims)
        set_temperature!(cng, kT)
        setfield!(integrator, sym, kT isa Real ? kT : cng.kT)
        check_noise_timestep(integrator.dt, integrator.kT)
    elseif sym == :damping
        setfield!(integrator, sym, validate_damping(val, dims))
    else
        setfield!(integrator, sym, val)
    end
end

# See `check_noise_timestep`.
function noise_timestep_bound(integrator::LangevinPlanck)
    kTmax = maximum(integrator.kT)
    return iszero(kTmax) ? Inf : 0.1/kTmax
end

function Base.show(io::IO, integrator::LangevinPlanck)
    (; dt, damping, kT) = integrator
    dt = isnan(integrator.dt) ? "<missing>" : repr(dt)
    println(io, "LangevinPlanck(sys, $dt; damping=$(param_string(damping)), kT=$(param_string(kT)))")
end

# Scale unit-normalized noise by the amplitude √(2λ) dt
scale_noise!(ζ, damping::Float64, dt) = (ζ .*= sqrt(2damping)*dt)
scale_noise!(ζ, damping::Array{Float64, 4}, dt) = (@. ζ *= sqrt(2damping)*dt)

@inline function rhs_dipole_pn!(ΔS, S, ξ, ∇E, integrator)
    (; dt, damping) = integrator
    λ = damping

    @. ΔS = - S × (ξ + dt*∇E - dt*λ*(S × ∇E))
end

@inline function advance_and_retrieve_noise!(sys, integrator)
    (; damping, noisesource, dt) = integrator
    cng = noisesource
    check_noise_dims(cng, size(sys.dipoles), 3)
    step_pn!(sys.rng, cng)
    ζ = view(reinterpret(SVector{3, Float64}, cng.ζ), 1, :, :, :, :)
    scale_noise!(ζ, damping, dt) # Note dt here -- treat as noise field
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
# Returns dt⋅X⋅Z for the noise matrix of `site`, where ζ[k + (site-1)N²] holds the k-th
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
        return complex(cdiag * ζ[noise_offset(a, N) + 1 + (site-1)*N^2])
    elseif a < b
        k = noise_offset(a, N) + 1 + 2(b - a) - 1
        j = k + (site-1)*N^2
        return coff * complex(ζ[j], ζ[j+1])
    else
        k = noise_offset(b, N) + 1 + 2(a - b) - 1
        j = k + (site-1)*N^2
        return coff * complex(ζ[j], -ζ[j+1])
    end
end

function step!(sys::System{N}, integrator::LangevinPlanck) where N
    check_timestep_available(integrator)
    # Function barrier for the Float64 or per-site `damping` (see `Langevin`)
    (; damping, dt, noisesource) = integrator
    step_planck!(sys, noisesource, dt, damping)
end

function step_planck!(sys::System{N}, noisesource, dt, damping) where N
    integrator = (; dt, damping)
    (Z′, ΔZ₁, ΔZ₂, ξ, HZ) = get_coherent_buffers(sys, 5)
    Z = sys.coherents

    check_noise_dims(noisesource, size(Z), N^2)
    step_pn!(sys.rng, noisesource)
    ζ = noisesource.ζ   # (N², dims...), indexed linearly; `reshape` would allocate

    # Euler prediction step. The noise term of `rhs_sun!` is -P ξ, so pass
    # ξ = i dt X Z.
    for i in eachindex(Z)
        ξ[i] = im * noise_field_times(ζ, i, Z[i], site_value(damping, i), dt)
    end
    set_energy_grad_coherents!(HZ, Z, sys)
    rhs_sun!(ΔZ₁, Z, ξ, HZ, integrator)
    @. Z′ = normalize_ket(Z + ΔZ₁, sys.κs)

    # Correction step. The multiplicative noise is re-evaluated at Z′, with the
    # same noise matrix (Stratonovich-consistent Heun).
    for i in eachindex(Z)
        ξ[i] = im * noise_field_times(ζ, i, Z′[i], site_value(damping, i), dt)
    end
    set_energy_grad_coherents!(HZ, Z′, sys)
    rhs_sun!(ΔZ₂, Z′, ξ, HZ, integrator)
    @. Z = normalize_ket(Z + (ΔZ₁+ΔZ₂)/2, sys.κs)

    # Coordinate dipole data
    sync_dipoles!(sys)

    return
end
