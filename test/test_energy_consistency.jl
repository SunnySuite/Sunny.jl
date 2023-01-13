@testitem "Energy consistency" begin
    include("shared.jl")


    function make_system(; SUN)
        cryst = Sunny.diamond_crystal()
        ints = Sunny.AbstractInteraction[]
        add_linear_interactions!(ints, SUN)
        add_quadratic_interactions!(ints, SUN)
        add_quartic_interactions!(ints, SUN)
        sys = SpinSystem(cryst, ints, (3, 3, 3); seed=0)
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

    test_conservation(make_system(; SUN=true))
    test_conservation(make_system(; SUN=false))


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

    test_delta(make_system(; SUN=true))
    test_delta(make_system(; SUN=false))
end
