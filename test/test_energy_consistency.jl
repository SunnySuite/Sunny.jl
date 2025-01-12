@testitem "Energy consistency" begin
    include("shared.jl")


    function make_system(; mode, inhomog)
        cryst = Sunny.diamond_crystal()
        sys = System(cryst, [1 => Moment(s=2, g=2)], mode; dims=(2, 2, 2), seed=0)
        add_linear_interactions!(sys, mode)
        add_quadratic_interactions!(sys, mode)
        add_quartic_interactions!(sys, mode)
        enable_dipole_dipole!(sys, 0.5)

        # Random field
        for site in eachsite(sys)
            set_field_at!(sys, 0.1*randn(sys.rng, 3), site)
        end
        # Random spin rescaling
        rand!(sys.rng, sys.κs)
        # Random spins
        randomize_spins!(sys)

        if !inhomog
            return sys
        else
            # Add some inhomogeneous interactions
            sys2 = to_inhomogeneous(sys)
            @test energy(sys2) ≈ energy(sys)
            set_vacancy_at!(sys2, (1,1,1,2))
            set_exchange_at!(sys2, 0.5, (1,1,1,1), (2,1,1,2); offset=(1, 0, 0))
            set_pair_coupling_at!(sys2, (Si, Sj) -> 0.7*(Si'*Sj)^2, (2,2,1,2), (2,1,1,3); offset=(0,-1,0))

            set_onsite_coupling_at!(sys2, S -> 0.4*(S[1]^4+S[2]^4+S[3]^4), (2,2,1,4))
            return sys2
        end
    end


    # Tests that SphericalMidpoint conserves energy to a certain tolerance.
    function test_conservation(sys)
        NITERS = 1500
        dt     = 0.002
        energies = Float64[]

        integrator = ImplicitMidpoint(dt)
        for _ in 1:NITERS
            step!(sys, integrator)
            push!(energies, energy(sys))
        end

        # Fluctuations should scale like square-root of system size,
        # independent of trajectory length
        sqrt_size = sqrt(length(sys.dipoles))
        ΔE = (maximum(energies) - minimum(energies)) / sqrt_size
        @test ΔE < 1e-3
    end

    test_conservation(make_system(; mode=:SUN, inhomog=false))
    test_conservation(make_system(; mode=:dipole, inhomog=true))


    # Tests that energy deltas are consistent with total energy
    function test_delta(sys)
        for _ in 1:5
            # Pick a random site, try to set it to a random spin
            site = rand(sys.rng, eachsite(sys))
            spin = Sunny.randspin(sys, site)

            ΔE = Sunny.local_energy_change(sys, site, spin)

            E0 = energy(sys)
            sys.dipoles[site]   = spin.S
            sys.coherents[site] = spin.Z
            E1 = energy(sys)
            ΔE_ref = E1 - E0

            @test isapprox(ΔE, ΔE_ref; atol=1e-12)
        end
    end

    test_delta(make_system(; mode=:SUN, inhomog=true))
    test_delta(make_system(; mode=:dipole, inhomog=false))
end
