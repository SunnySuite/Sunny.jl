@testitem "Spin Scaling" begin
    include("shared.jl")

    function make_exchange_interactions()
        J  = 1.0   # Anti-ferro nearest neighbor
        J′ = -1.0  # Ferro next-nearest neighbor
        K  = 1.0   # Scale of Kitaev term
        Γ  = 0.0   # Off-diagonal exchange, not used
        J_exch = [J     Γ   0.0;
                Γ     J   0.0;
                0.0  0.0  J+K]
        return [exchange(J_exch, Bond(1, 2, [0,0,0])),
                heisenberg(J′, Bond(1, 1, [1,0,0]))]
    end


    function make_test_system_lld(; spin_rescaling=1.0)
        cryst = Sunny.fcc_crystal()

        # Exchange interactions
        exchange_interactions = make_exchange_interactions()

        # Quartic anisotropy
        D = 1.0 
        quartic_interactions = [anisotropy(D*(𝒮[1]^4+𝒮[2]^4+𝒮[3]^4), 1, "quartic")]

        interactions_all = [exchange_interactions..., quartic_interactions...]
        dims = (3,3,3)

        return SpinSystem(cryst,
                        interactions_all,
                        dims,
                        [SiteInfo(1; spin_rescaling)]
        )
    end


    function make_test_system_gsd(; spin_rescaling=1.0, N=2)
        cryst = Sunny.fcc_crystal()

        # Exchange interactions
        exchange_interactions = make_exchange_interactions()

        # Quartic anisotropy
        quartic_sun = anisotropy(-𝒮[3]^4, 1, "quartic")

        dims = (3,3,3)
        interactions_all = [exchange_interactions..., quartic_sun]

        return SpinSystem(cryst,
                        interactions_all,
                        dims,
                        [SiteInfo(1; N, spin_rescaling)]
        )
    end

    const spin_rescalings = [0.2, 2.1]

    function spin_magnitude_stability_tester(sys_maker, integrators)
        for integrator in integrators
            for spin_rescaling in spin_rescalings
                sys = sys_maker(; spin_rescaling)
                rand!(sys)
                mags = norm.(sys.dipoles)
                for _ ∈ 1:100
                    step!(sys, integrator)
                end
                @test mags ≈ norm.(sys.dipoles)
            end
        end
    end

    function test_spin_magnitude_stability()
        kT = 0.1
        λ  = 0.1
        Δt = 0.01

        integrators = [LangevinHeunP(kT, λ, Δt), ImplicitMidpoint(Δt)]

        spin_magnitude_stability_tester(make_test_system_lld, integrators)
        spin_magnitude_stability_tester(make_test_system_gsd, integrators)
    end

    test_spin_magnitude_stability()


    function test_energy_scaling()
        Ns = [0, 5, 6]

        cryst = Sunny.fcc_crystal()
        dims = (2,2,2)

        Λ = 𝒪[4,0]+5𝒪[4,4]

        interactions = [heisenberg(1.0, Bond(1,2,[0,0,0])),
                        biquadratic(1.1, Bond(1,2,[0,0,0])),
                        anisotropy(Λ, 1)]

        for N in Ns
            # In SU(N) mode, rescaling is applied to expectation values, ⟨A⟩ ->
            # spin_rescaling ⟨A⟩, for arbitrary operator A.
            powers = [2, 4, N==0 ? 4 : 1]

            for (interaction, power) in zip(interactions, powers)
                for spin_rescaling in spin_rescalings
                    sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1; N)]; seed=0)
                    rand!(sys)
                    E₀ = energy(sys)

                    sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1; N, spin_rescaling)]; seed=0)
                    rand!(sys)
                    E₁ = energy(sys)

                    @test E₁ ≈ spin_rescaling^power * E₀
                end
            end
        end
    end

    test_energy_scaling()

    """Generates a trajectory for a single spin in the presence of an 
    external magnetic field. Rescales resulting spin magnitude so trajectories
    with different scalings can be directly compared.
    """
    function generate_scaled_zeeman_trajectory(spin_rescaling, θ, Δt; N=0, dur=10.0)
        cryst = Sunny.cubic_crystal()
        dims = (1,1,1)
        interactions = [external_field([0.0, 0.0, 10.0])]

        sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N, spin_rescaling)])
        spin_mag = spin_rescaling * (N == 0 ? 1 : (N-1)/2)
        Sunny.set_dipole!(sys, CartesianIndex((1,1,1,1)), spin_mag * [0, sin(θ), cos(θ)])

        integrator = ImplicitMidpoint(Δt)

        numsteps = round(Int, dur/Δt) 
        ts = (0:numsteps) .* Δt
        S = zeros(Sunny.Vec3, numsteps+1)
        S[1] = sys.dipoles[1]

        for i in 1:numsteps
            step!(sys, integrator)
            S[i+1] = sys.dipoles[1]
        end

        return (;
            xs = [S[1]/spin_rescaling for S ∈ S],
            ys = [S[1]/spin_rescaling for S ∈ S],
            zs = [S[1]/spin_rescaling for S ∈ S],
            ts
        ) 
    end

    """Tests invariance of spin dynamics under spin rescaling 
    in the presence of a Zeeman term. Tests both LLD and GSD. 
    """
    function test_scaling_zeeman()
        Δt = 0.001
        θ = (π/4 - π/32)*rand() + π/32  # amount to tilt spin in zy-plane
        spin_rescaling = 2.1
        Ns = [0, 2]

        for N ∈ Ns
            (; xs) = generate_scaled_zeeman_trajectory(1.0, θ, Δt; N)
            xs_1 = xs
            (; xs) = generate_scaled_zeeman_trajectory(spin_rescaling, θ, Δt; N)
            xs_2 = xs

            rms = √sum( (xs_2 .- xs_1) .^2 )

            @test rms < 1e-10 
        end
    end

    test_scaling_zeeman()

    """Generate a trajectory for a system with only quadratic interactions. Results are rescaled 
    so results with different spin magnitudes can be compared directly.
    """
    function generate_scaled_quadratic_trajectory(spin_rescaling, Δt; N=0, dur=10.0)
        cryst = Sunny.cubic_crystal()
        dims = (4,4,3)
        interactions = Sunny.AbstractInteraction[
            heisenberg(1.0, Bond(1,1,[1,0,0])),
        ]
        if N == 0   # "Quadratic anisotropy" only scales quadratically for old dynamics
            push!(interactions, quadratic_anisotropy(1.0*I(3), 1))
        end

        sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N, spin_rescaling)]; seed=111)
        enable_dipole_dipole!(sys)

        rand!(sys)

        integrator = ImplicitMidpoint(Δt/spin_rescaling)

        numsteps = round(Int, dur/Δt) 
        ts = (0:numsteps) .* Δt
        S = zeros(Sunny.Vec3, numsteps+1)
        S[1] = sys.dipoles[1]

        for i in 1:numsteps
            step!(sys, integrator)
            S[i+1] = sys.dipoles[1]
        end

        return (;
            xs = [S[1]/spin_rescaling for S ∈ S],
            ys = [S[1]/spin_rescaling for S ∈ S],
            zs = [S[1]/spin_rescaling for S ∈ S],
            ts 
        ) 
    end

    """Test invariance of dynamics (with Hamiltonian that is quadratic in spins) under 
    the rescaling of spin magnitudes.
    """
    function test_scaling_quadratic()
        Δt = 0.01
        spin_rescaling = 2.1
        Ns = [0, 2]

        for N ∈ Ns
            (; xs) = generate_scaled_quadratic_trajectory(1.0, Δt; N)
            xs_1 = xs
            (; xs) = generate_scaled_quadratic_trajectory(spin_rescaling, Δt; N)
            xs_2 = xs

            rms = √sum( (xs_2 .- xs_1) .^2 )

            @test rms < 1e-8
        end
    end

    test_scaling_quadratic()

end
