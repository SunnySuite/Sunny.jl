function spiral_energy(sys::System{0}; k, axis)
    sys.mode in (:dipole, :dipole_large_S) || error("SU(N) mode not supported")
    sys.latsize == (1, 1, 1) || error("System must have only a single cell")

    check_rotational_symmetry(sys; axis, θ=0.01)

    E, _dEdk = spiral_energy_and_gradient_aux!(nothing, sys; k, axis)
    return E
end

function spiral_energy_per_site(sys::System{0}; k, axis)
    return spiral_energy(sys; k, axis) / natoms(sys.crystal)
end

function spiral_energy_and_gradient_aux!(dEds, sys::System{0}; k, axis)
    E = 0
    accum_grad = !isnothing(dEds)
    if accum_grad
        fill!(dEds, zero(Vec3))
    end
    dEdk = zero(Vec3)

    @assert sys.latsize == (1,1,1)
    Na = natoms(sys.crystal)

    x, y, z = normalize(axis)
    K = Sunny.Mat3([0 -z y; z 0 -x; -y x 0])
    K² = K*K

    for i in 1:Na
        (; onsite, pair) = sys.interactions_union[i]
        Si = sys.dipoles[i]

        # Pair coupling
        for coupling in pair
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break
            (; j, n) = bond
            Sj = sys.dipoles[j]

            # Rotation angle along `axis` for cells displaced by `n`
            θ = 2π * dot(k, n)
            dθdk = 2π*n

            # Rotation as a 3×3 matrix
            s, c = sincos(θ)
            R = I + s*K + (1-c)*K²
            @assert R ≈ axis_angle_to_matrix(axis, θ)
            dRdθ = c*K + s*K²

            # J is invariant under any rotation along `axis`
            J = Mat3(bilin*I)
            @assert R'*J*R ≈ J

            # Accumulate energy and derivatives
            E += Si' * J * (R * Sj)
            @assert Si' * J * (R * Sj) ≈ (R' * Si)' * J * Sj
            if accum_grad
                dEds[i] += J * (R * Sj)
                dEds[j] += J' * (R' * Si)
            end
            dEdθ = Si' * J * (dRdθ * Sj)
            dEdk += dEdθ * dθdk

            @assert iszero(biquad) "Biquadratic interactions not supported"
        end

        # Onsite coupling
        E_aniso, dEds_aniso = energy_and_gradient_for_classical_anisotropy(Si, onsite)
        E += E_aniso

        # Zeeman coupling
        E += sys.extfield[i]' * (sys.units.μB * sys.gs[i] * Si)

        if accum_grad
            dEds[i] += dEds_aniso
            dEds[i] += sys.units.μB * sys.gs[i]' * sys.extfield[i]
        end
    end

    # See "spiral_energy.lyx" for derivation
    if !isnothing(sys.ewald)
        μ = [magnetic_moment(sys, site) for site in eachsite(sys)]

        A0 = sys.ewald.A
        A0 = reshape(A0, Na, Na)

        Ak = Sunny.precompute_dipole_ewald_at_wavevector(sys.crystal, (1,1,1), sys.units.μ0, k)
        Ak = reshape(Ak, Na, Na)

        ϵ = 1e-8
        if norm(k - round.(k)) < ϵ
            for i in 1:Na, j in 1:Na
                E += real(μ[i]' * A0[i, j] * μ[j]) / 2
            end
        elseif norm(2k - round.(2k)) < ϵ
            for i in 1:Na, j in 1:Na
                E += real(μ[i]' * ((I+K²)*A0[i, j]*(I+K²) + K²*Ak[i, j]*K²) * μ[j]) / 2
            end
        else
            for i in 1:Na, j in 1:Na
                E += real(μ[i]' * ((I+K²)*A0[i, j]*(I+K²) + (im*K+K²)*Ak[i, j]*(im*K+K²)/2) * μ[j]) / 2
            end
        end

        if accum_grad
            error("Cannot yet differentiate through Ewald summation")
        end
    end

    return E, dEdk
end

# Sets sys.dipoles and returns k, according to data in params
function unpack_spiral_params!(sys::System{0}, axis, params)
    params = reinterpret(Vec3, params)
    L = length(sys.dipoles)
    for i in 1:L
        u = stereographic_projection(params[i], axis)
        sys.dipoles[i] = sys.κs[i] * u
    end
    return params[end]
end

# Regularizer that blows up (by factor of 10) when |x| → 1. Use x=u⋅axis to
# favor normalized spins `u` orthogonal to `axis`.
reg(x) = 1 / (1 - x^2 + 1/10)
dreg(x) = 2x * reg(x)^2

function spiral_f(sys::System{0}, axis, params, λ)
    k = unpack_spiral_params!(sys, axis, params)
    E, _dEdk = spiral_energy_and_gradient_aux!(nothing, sys; k, axis)
    for s in sys.dipoles
        u = normalize(s)
        E += λ * reg(u⋅axis)
    end
    return E
end

function spiral_g!(G, sys::System{0}, axis, params, λ)
    k = unpack_spiral_params!(sys, axis, params)
    v = reinterpret(Vec3, params)
    G = reinterpret(Vec3, G)

    L = length(sys.dipoles)
    dEds = view(G, 1:L)
    _E, dEdk = spiral_energy_and_gradient_aux!(dEds, sys; k, axis)

    for i in 1:L
        s = sys.dipoles[i]
        u = normalize(s)
        # dE/du' = dE/ds' * ds/du, where s = |s|*u.
        dEdu = dEds[i] * norm(s) + λ * dreg(u⋅axis) * axis
        # dE/dv' = dE/du' * du/dv
        G[i] = vjp_stereographic_projection(dEdu, v[i], axis)
    end
    G[end] = dEdk
end

function minimize_energy_spiral!(sys, axis; maxiters=10_000, k_guess=randn(sys.rng, 3))
    axis = normalize(axis)

    sys.mode in (:dipole, :dipole_large_S) || error("SU(N) mode not supported")
    sys.latsize == (1, 1, 1) || error("System must have only a single cell")
    norm([s × axis for s in sys.dipoles]) > 1e-12 || error("Spins cannot be exactly aligned with polarization axis")

    # Note: if k were fixed, we could check θ = 2πkᵅ for each component α, which
    # is a weaker constraint.
    check_rotational_symmetry(sys; axis, θ=0.01)

    L = natoms(sys.crystal)
    
    params = fill(zero(Vec3), L+1)
    for i in 1:L
        params[i] = inverse_stereographic_projection(normalize(sys.dipoles[i]), axis)
    end
    params[end] = k_guess

    local λ::Float64
    f(params) = spiral_f(sys, axis, params, λ)
    g!(G, params) = spiral_g!(G, sys, axis, params, λ)

    # Minimize f, the energy of a spiral
    options = Optim.Options(; iterations=maxiters)

    # LBFGS does not converge to high precision, but ConjugateGradient can fail
    # to converge: https://github.com/JuliaNLSolvers/LineSearches.jl/issues/175.
    # TODO: Call only ConjugateGradient when issue is fixed.
    method = Optim.LBFGS(; linesearch=Optim.LineSearches.BackTracking(order=2))
    λ = 1 * abs(spiral_energy_per_site(sys; k=k_guess, axis)) # regularize at some energy scale
    res0 = Optim.optimize(f, g!, collect(reinterpret(Float64, params)), method, options)
    λ = 0 # disable regularization
    res = Optim.optimize(f, g!, Optim.minimizer(res0), Optim.ConjugateGradient(), options)

    k = unpack_spiral_params!(sys, axis, Optim.minimizer(res))

    if Optim.converged(res)
        # For aesthetics, wrap k components to [1-ϵ, -ϵ)
        return wrap_to_unit_cell(k; symprec=1e-6)
    else
        println(res)
        error("Optimization failed to converge within $maxiters iterations.")
    end
end
