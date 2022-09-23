@testitem "Metropolis Sampling" begin
include("test_shared.jl")

function test_local_energy_change()
    system = produce_example_system()

    for _ in 1:3
        rand!(system)
        for _ in 1:50
            # Pick a random site, try to set it to a random spin
            randsite = CartesianIndex(Tuple(mod1.(rand(Int, 4), size(system))))
            N = system.site_infos[1].N 
            newspin = Sunny._random_spin(system.rng, Val(N))

            func_diff = Sunny.local_energy_change(system, randsite, newspin)

            orig_energy = energy(system)
            system._dipoles[randsite] = newspin
            new_energy = energy(system)

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
    system = produce_example_system()
    
    Random.seed!(1111)
    rand_kTs = rand(100)

    for sampler_type in [MetropolisSampler, IsingSampler]
        sampler = sampler_type(system, 1.0, 1)
        for kT in rand_kTs
            set_temp!(sampler, kT)
            # The approx here is due to numerical issues: the temp returned by
            #  get_temp is actually 1 / (1 / kT) from set_temp!.
            if !(get_temp(sampler) ≈ kT)
                return false
            end
        end
    end

    return true
end

@test test_set_get_temp_metropolis()

end