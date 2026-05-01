import NonlinearSolve: NonlinearProblem, solve, NewtonRaphson

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

@inline n(ω, kT) = 1/(exp(ω/kT) - 1)
planck_spectrum(ω, kT) = iszero(ω) ? kT : ω*n(ω, kT) 
planck_spectrum_sym(ω, kT) = ω*coth(ω/(2kT)) 

# Analytical expression for the power spectrum of a second-order linear filter
# with parameters p.
function filter_spectrum(ω, p)
    c1, c2, Ω1, Ω2, Γ1, Γ2 = p
    @. ((2c1^2 * Γ1) / ((Ω1^2 - ω^2)^2 + ω^2 * Γ1^2)) + ((2c2^2 * Γ2) / ((Ω2^2 - ω^2)^2 + ω^2 * Γ2^2))
end


################################################################################
# Colored noise
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
# Used to determine range of ω values for fitting the filter responses.
# Note that this can be solved analytically for an arbitrary percentage,
# but this would rely on an implementation of the Lambert W function, which is
# not included in SpecialFunctions. The constant used here is specifically for
# acheiving a value that is 1% of kT.
function ω_cutoff(kT) 
    return max(6.4746008706*kT, 5.0)
end

@. quad_model(x, p) = p[1] + p[2]*x + p[3]*x*x

function colored_noise_params(kT)
    lim = ω_cutoff(kT)
    ωs = range(0.0, lim, 1000)
    ys = planck_spectrum.(ωs, kT)

    # Seed optimization correctly -- parameters fit as a function of temperature
    c10 = quad_model(kT,     [0.0, 0.0,      0.3] ) 
    c20 = quad_model(kT,     [0.0, 0.0,      1.88] )
    omega10 = quad_model(kT, [0.0, 1.168055, 0.0])
    omega20 = quad_model(kT, [0.0, 2.748380, 0.0])
    gamma10 = quad_model(kT, [0.0, 3.276618, 0.0])
    gamma20 = quad_model(kT, [0.0, 5.247509, 0.0])

    # fit = curve_fit(filter_spectrum, ωs, ys, [c10, c20, omega10, omega20, gamma10, gamma20])
    params = lsq_fit(filter_spectrum, ωs, ys, [c10, c20, omega10, omega20, gamma10, gamma20])
    c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂ = if all(>(0.0), params)
        params
    else
        error("Didn't find good parameters")
    end

    return (; c₁,  c₂, Ω₁, Ω₂, Γ₁, Γ₂)
end

function PlanckNoiseGenerator(dt; kT, damping, dims)
    ## Savin/Barker parameters
    # c₁, c₂ = 1.8315, 0.3429
    # Ω₁, Γ₁ = 2.7189, 5.0142
    # Ω₂, Γ₂ = 1.2223, 3.2974

    c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂ = colored_noise_params(kT) 
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

    # Determine new coefficients for noise process
    ## Savin/Barker parameters
    # c₁, c₂ = 1.8315, 0.3429
    # Ω₁, Γ₁ = 2.7189, 5.0142
    # Ω₂, Γ₂ = 1.2223, 3.2974

    c₁, c₂, Ω₁, Ω₂, Γ₁, Γ₂ = colored_noise_params(kT) 
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
# Langevin integration with colored noise 
################################################################################
mutable struct LangevinPlanck <: AbstractIntegrator
    dt              :: Float64
    damping         :: Float64
    kT              :: Float64
    noisesource     :: PlanckNoiseGenerator

    function LangevinPlanck(dt; λ=nothing, damping=nothing, kT, sysdims)
        if !isnothing(λ)
            @warn "`λ` argument is deprecated! Use `damping` instead."
            damping = @something damping λ
        end
        isnothing(damping) && error("`damping` parameter required")
        iszero(damping) && error("Use ImplicitMidpoint instead for energy-conserving dynamics")

        dt <= 0         && error("Select positive dt")
        kT < 0          && error("Select nonnegative kT")
        damping <= 0    && error("Select positive damping")


        cng = PlanckNoiseGenerator(dt; kT, damping, dims=sysdims)
        return new(dt, damping, kT, cng)
    end
end

function LangevinPlanck(; λ=nothing, damping=nothing, kT, dims)
    LangevinPlanck(NaN; λ, damping, kT, dims)
end

function Base.copy(dyn::LangevinPlanck)
    LangevinPlanck(dyn.dt; dyn.damping, dyn.kT)
end

function Base.setproperty!(integrator::LangevinPlanck, sym::Symbol, val)
    if sym == :dt
        setfield!(integrator, sym, val)
        setfield!(integrator.noisesource, sym, val)
    else
        setfield!(integrator, sym, val)
    end
end

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
    step_pn!(sys.rng, cng)
    ζ = view(reinterpret(SVector{3, Float64}, cng.ζ), 1, :, :, :, :)
    ζ .*= sqrt(2damping)*dt # Note dt here -- treat as noise field
    return ζ
end

function step!(sys::System{0}, integrator::LangevinPlanck)
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

function step!(sys::System{N}, integrator::LangevinPlanck) where N
    (Z′, ΔZ₁, ΔZ₂, ζ, HZ) = get_coherent_buffers(sys, 5)
    Z = sys.coherents

    ζ = advance_and_retrieve_noise!(sys, integrator)

    # Euler prediction step
    set_energy_grad_coherents!(HZ, Z, sys)
    rhs_sun!(ΔZ₁, Z, ζ, HZ, integrator)
    @. Z′ = normalize_ket(Z + ΔZ₁, sys.κs)

    # Correction step
    set_energy_grad_coherents!(HZ, Z′, sys)
    rhs_sun!(ΔZ₂, Z′, ζ, HZ, integrator)
    @. Z = normalize_ket(Z + (ΔZ₁+ΔZ₂)/2, sys.κs)

    # Coordinate dipole data
    @. sys.dipoles = Sunny.expected_spin(Z)

    return
end