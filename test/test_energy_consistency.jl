@testitem "Energy consistency" begin
    include("shared.jl")


    function make_system(; mode, inhomog)
        cryst = Sunny.diamond_crystal()
        sys = System(cryst, (3, 3, 3), [SpinInfo(1, S=2)], mode; seed=0)
        add_linear_interactions!(sys, mode)
        add_quadratic_interactions!(sys, mode)
        add_quartic_interactions!(sys, mode)
        enable_dipole_dipole!(sys)

        # Random field
        for site in all_sites(sys)
            set_external_field_at!(sys, randn(sys.rng, 3), site)
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
            set_vacancy_at!(sys2, (1,1,1,1))
            set_exchange_at!(sys2, 0.5, (1,1,1,1), (2,1,1,2))
            set_biquadratic_at!(sys2, 0.7, (3,2,1,2), (3,1,1,3))
            set_anisotropy_at!(sys2, 0.4*(𝒮[1]^4+𝒮[2]^4+𝒮[3]^4), (2,2,2,4))
            return sys2
        end
    end


    # Tests that SphericalMidpoint conserves energy to a certain tolerance.
    function test_conservation(sys)
        NITERS = 5_000
        Δt     = 0.002
        energies = Float64[]

        integrator = ImplicitMidpoint(Δt)
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
            site = rand(sys.rng, all_sites(sys))
            spin = Sunny.randspin(sys, site)
            
            ΔE = Sunny.local_energy_change(sys, site, spin)

            E0 = energy(sys)
            sys.dipoles[site]   = spin.s
            sys.coherents[site] = spin.Z
            E1 = energy(sys)
            ΔE_ref = E1 - E0

            @test isapprox(ΔE, ΔE_ref; atol=1e-12)
        end
    end

    test_delta(make_system(; mode=:SUN, inhomog=true))
    test_delta(make_system(; mode=:dipole, inhomog=false))
end
