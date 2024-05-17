function spiral_energy(sys::System, k, axis; exchange_only=false, check_symmetry=true)
    @assert sys.mode in (:dipole, :dipole_large_S) "SU(N) mode not supported"
    @assert sys.latsize == (1, 1, 1) "System must have only a single cell"
    
    # Optionally disable symmetry check for speed
    check_symmetry && check_rotational_symmetry(sys; axis, θ=0.01)

    E = 0
    for i in 1:natoms(sys.crystal)
        (; onsite, pair) = sys.interactions_union[i]
        Si = sys.dipoles[i]

        # Pair coupling
        for coupling in pair
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break
            (; j, n) = bond
            # R along `axis` for cells displaced by `n`
            θ = 2π * dot(k, n)
            R = axis_angle_to_matrix(axis, θ)

            # J is invariant under any `axis` rotations
            J = Mat3(bilin*I)
            @assert R'*J*R ≈ J

            Sj = sys.dipoles[j]
            E += Si' * (J * R) * Sj # == (R * Si)' * (J * R) * (R * Sj)

            @assert iszero(biquad) "Biquadratic interactions not supported"
        end

        if !exchange_only
            # Onsite coupling
            E += energy_and_gradient_for_classical_anisotropy(Si, onsite)[1]

            # Zeeman coupling
            E -= sys.extfield[1, 1, 1, i]' * magnetic_moment(sys, (1, 1, 1, i))
        end
    end

    return E
end

function minimize_energy_spiral_aux!(sys, axis; maxiters, g_tol, k_guess, allow_canting)
    nspin = natoms(sys.crystal)
    
    # Any rotation that maps spins into "local" frame: R * axis = [0, 0, 1]
	R = rotation_between_vectors(axis, [0, 0, 1])

    # Polar and azimuthal angle in local frame
    ϕs = Float64[]
    θs = Float64[]
    for s in sys.dipoles
        s′ = R * s
        push!(ϕs, atan(s′[2], s′[1]))
        push!(θs, angle_between_vectors(s′, Vec3(0, 0, 1)))
    end

    # Full set of optimization parameters
    params = [ϕs; θs; k_guess]

    function set_dipoles_from_params!(params)
        ϕ = view(params, 1:nspin)
        θ = view(params, nspin+1:2nspin)
        for i in 1:nspin
            sinϕ, cosϕ = sincos(ϕ[i])
            sinθ, cosθ = allow_canting ? sincos(θ[i]) : (1.0, 0.0)
            S = R' * Vec3(sinθ*cosϕ, sinθ*sinϕ, cosθ)
            set_dipole!(sys, S, (1,1,1,i))
        end
    end

    # Optimize to the full system energy
    options = Optim.Options(; g_tol, iterations=maxiters)
    res = Optim.optimize(params, Optim.ConjugateGradient(), options) do params
        set_dipoles_from_params!(params)
        k = Vec3(params[end-2],params[end-1], params[end])
        spiral_energy(sys, k, axis; check_symmetry=false)
    end

    params = Optim.minimizer(res)
    set_dipoles_from_params!(params)

    if Optim.converged(res)
        # println(res)

        # Return optimal k. For aesthetics, wrap components to [1-ϵ, -ϵ)
        k = Vec3(params[end-2], params[end-1], params[end])
        return wrap_to_unit_cell(k; symprec=1e-7)
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

