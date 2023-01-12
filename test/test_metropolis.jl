@testitem "Delta energy consistency" begin
    include("shared.jl")

    function test_local_energy_change(; SUN)
        cryst = Sunny.diamond_crystal()
        ints = Sunny.AbstractInteraction[]
        add_linear_interactions!(ints, SUN)
        add_quadratic_interactions!(ints, SUN)
        add_quartic_interactions!(ints, SUN)
        sys = SpinSystem(cryst, ints, (5, 5, 5); seed=0)
        enable_dipole_dipole!(sys)

        rand!(sys.rng, sys.κs)
        randomize_spins!(sys)

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

    test_local_energy_change(; SUN=false)
    test_local_energy_change(; SUN=true)

end
