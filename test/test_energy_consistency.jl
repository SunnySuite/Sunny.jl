@testitem "Energy consistency" begin
    include("shared.jl")


    function make_system(; mode)
        cryst = Sunny.diamond_crystal()
        sys = SpinSystem(cryst, (3, 3, 3); mode, seed=0)
        add_linear_interactions!(sys, mode)
        add_quadratic_interactions!(sys, mode)
        add_quartic_interactions!(sys, mode)
        enable_dipole_dipole!(sys)

        rand!(sys.rng, sys.κs)
        randomize_spins!(sys)
        return sys
    end


    # Tests that SphericalMidpoint conserves energy to a certain tolerance.
    function test_conservation(sys)
        NITERS = 5_000
        Δt     = 0.005
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
        @test ΔE < 1e-2
    end

    test_conservation(make_system(; mode=:SUN))
    test_conservation(make_system(; mode=:dipole))


    # Tests that energy deltas are consistent with total energy
    function test_delta(sys)
        for _ in 1:5
            # Pick a random site, try to set it to a random spin
            idx = rand(sys.rng, CartesianIndices(sys.dipoles))
            state = Sunny.random_state(sys, idx)
            
            ΔE = Sunny.local_energy_change(sys, idx, state)

            E0 = energy(sys)
            sys.dipoles[idx]   = state.s
            sys.coherents[idx] = state.Z
            E1 = energy(sys)
            ΔE_ref = E1 - E0

            @test isapprox(ΔE, ΔE_ref; atol=1e-12)
        end
    end

    test_delta(make_system(; mode=:SUN))
    test_delta(make_system(; mode=:dipole))
end
