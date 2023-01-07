@testitem "Metropolis Sampling" begin
    include("shared.jl")

    function test_local_energy_change()
        SUN = false
        cryst = Sunny.diamond_crystal()
        ints = Sunny.AbstractInteraction[]
        add_linear_interactions!(ints, SUN)
        add_quadratic_interactions!(ints, SUN)
        sys = SpinSystem(cryst, ints, (5, 5, 5); seed=1111)

        for _ in 1:3
            randomize_spins!(sys)
            for _ in 1:50
                # Pick a random site, try to set it to a random spin
                idx = rand(CartesianIndices(sys.dipoles))
                newspin = Sunny._random_spin(sys, idx)

                func_diff = Sunny.local_energy_change(sys, idx, newspin)

                orig_energy = energy(sys)
                sys.dipoles[idx] = newspin
                new_energy = energy(sys)

                actual_diff = new_energy - orig_energy

                if !(func_diff ≈ actual_diff)
                    return false
                end
            end
        end

        return true
    end

    @test test_local_energy_change()

    "Tests that set_temp!/get_temp behave as expected"
    function test_set_get_temp_metropolis()
        SUN = false
        cryst = Sunny.diamond_crystal()
        ints = Sunny.AbstractInteraction[]
        add_linear_interactions!(ints, SUN)
        add_quadratic_interactions!(ints, SUN)
        sys = SpinSystem(cryst, ints, (5, 5, 5); seed=1111)

        for sampler_type in [MetropolisSampler, IsingSampler]
            sampler = sampler_type(sys, 1.0, 1)
            for kT in [1.1, 10.3]
                set_temp!(sampler, kT)
                # approximate because `sampler` stores 1/kT
                if !(get_temp(sampler) ≈ kT)
                    return false
                end
            end
        end

        return true
    end

    @test test_set_get_temp_metropolis()

end