function spiral_energy(sys::System{0}, k, axis; exchange_only=false, check_symmetry=true)
    @assert sys.mode in (:dipole, :dipole_large_S) "SU(N) mode not supported"
    @assert sys.latsize == (1, 1, 1) "System must have only a single cell"

    # Optionally disable symmetry check for speed
    check_symmetry && check_rotational_symmetry(sys; axis, θ=0.01)

    E, _dEdk = spiral_energy_and_gradient!(nothing, sys, k, axis; exchange_only)
    return E
end


function spiral_energy_and_gradient!(dEds, sys::System{0}, k, axis; exchange_only)
    E = 0
    accum_grad = !isnothing(dEds)
    if accum_grad
        fill!(dEds, zero(Vec3))
    end
    dEdk = zero(Vec3)

    for i in 1:natoms(sys.crystal)
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
            x, y, z = normalize(axis)
            s, c = sincos(θ)
            K = [0 -z y; z 0 -x; -y x 0]
            K² = K*K
            R = I + s*K + (1-c)*K²
            @assert R ≈ axis_angle_to_matrix(axis, θ)
            dRdθ = c*K + s*K²

            # J is invariant under any rotation along `axis`
            J = Mat3(bilin*I)
            @assert R'*J*R ≈ J

            # Accumulate energy and derivatives
            E += Si' * (J * R) * Sj # ≈ (R * Si)' * (J * R) * (R * Sj)
            if accum_grad
                dEds[i] += (J * R) * Sj
                dEds[j] += (J * R)' * Si
            end
            dEdθ = Si' * (J * dRdθ) * Sj
            dEdk += dEdθ * dθdk

            @assert iszero(biquad) "Biquadratic interactions not supported"
        end

        if !exchange_only
            # Onsite coupling
            E_aniso, dEds_aniso = energy_and_gradient_for_classical_anisotropy(Si, onsite)
            E += E_aniso

            # Zeeman coupling
            E -= sys.extfield[i]' * (sys.units.μB * sys.gs[i] * Si)

            if accum_grad
                dEds[i] += dEds_aniso
                dEds[i] -= sys.units.μB * sys.gs[i]' * sys.extfield[i]
            end
        end
    end

    return E, dEdk
end

function optim_set_spins_spiral!(sys::System{0}, axis, params)
    L = length(sys.dipoles)
    for i in 1:L
        u = stereographic_projection(params[i], axis)
        sys.dipoles[i] = sys.κs[i] * u
    end
end


function optim_set_gradient_spiral!(G, sys::System{0}, axis, params)
    params = reinterpret(Vec3, params)
    G = reinterpret(Vec3, G)
    optim_set_spins_spiral!(sys, axis, params)
    k = params[end]

    L = length(sys.dipoles)
    dEds = view(G, 1:L)
    _E, dEdk = spiral_energy_and_gradient!(dEds, sys, k, axis; exchange_only=false)
    G[end] = dEdk

    for i in 1:L
        # dE/du' = dE/ds' * ds/du, where s = |s|*u.
        dEdu = dEds[i] * sys.κs[i]
        # dE/dα' = dE/du' * du/dα
        G[i] = vjp_stereographic_projection(dEdu, params[i], axis)
    end

    println("Evaluating gradient at k = $k with G = $(G[end])")
end

function optim_energy_spiral(sys::System{0}, axis, params)
    params = reinterpret(Vec3, params)
    optim_set_spins_spiral!(sys, axis, params)
    k = params[end]
    E, _dEdk = spiral_energy_and_gradient!(nothing, sys, k, axis; exchange_only=false)
    println("Evaluating energy at k = $k with E = $E")
    return E
end


function minimize_energy_spiral_aux!(sys, axis; maxiters, g_tol, k_guess, allow_canting)
    L = natoms(sys.crystal)
    
    params = fill(zero(Vec3), L+1)
    for i in 1:L
        params[i] = inverse_stereographic_projection(normalize(sys.dipoles[i]), axis)
    end
    params[end] = k_guess

    # println("Using params ", params)

    f(params) = optim_energy_spiral(sys, axis, params)
    g!(G, params) = optim_set_gradient_spiral!(G, sys, axis, params)

    function g(params)
        G = copy(params)
        optim_set_gradient_spiral!(G, sys, axis, params)
        return G
    end

    # αs = range(0, 0.02, length=5)
    # G = g(params)
    # xs = [params - α*G for α in αs]
    # println("Test eval f = $(f(params)) g = $(g(params))")
    # println("Perturb = $(f.(xs))")

    # Optimize to the full system energy
    options = Optim.Options(; iterations=maxiters, show_trace=true)
    res = Optim.optimize(f, g!, collect(reinterpret(Float64, params)), Optim.ConjugateGradient(), options)

    params = reinterpret(Vec3, Optim.minimizer(res))
    optim_set_spins_spiral!(sys, axis, params)

    if Optim.converged(res)
        # println(res)

        # Return optimal k. For aesthetics, wrap components to [1-ϵ, -ϵ)
        return wrap_to_unit_cell(params[end]; symprec=1e-7)
    else
        println(res)
        error("Optimization failed to converge within $maxiters iterations.")
    end
end


function minimize_energy_spiral!(sys, axis; maxiters=10_000, g_tol=1e-10, k_guess=randn(sys.rng, 3))
    @assert sys.mode in (:dipole, :dipole_large_S) "SU(N) mode not supported"
    @assert sys.latsize == (1, 1, 1) "System must have only a single cell"

    # TODO: If k were fixed, we could check θ = 2πkᵅ for each component α, which
    # may be a weaker constraint.
    check_rotational_symmetry(sys; axis, θ=0.01)

    # First perform an optimization step where canting is dissallowed. This
    # avoids the singularity for polar angles θ = {0, π}, and seems to make
    # convergence more robust. Alternatively, could use a stereographic project
    # as in minimize_energy!.
    k_guess = minimize_energy_spiral_aux!(sys, axis; maxiters, g_tol, k_guess, allow_canting=false)
    return minimize_energy_spiral_aux!(sys, axis; maxiters, g_tol, k_guess, allow_canting=true)
end

