@testset "Metropolis Sampling" begin

function test_local_energy_change()
    cryst = Sunny.diamond_crystal()
    latsize = [5, 5, 5]

    interactions = [
        diamond_test_exchanges()...,
        external_field([0, 0, 1])
    ]

    system = SpinSystem(cryst, interactions, latsize)

    for _ in 1:3
        rand!(system)
        for _ in 1:50
            # Pick a random site, try to set it to a random spin
            randsite = CartesianIndex(Tuple(mod1.(rand(Int, 4), size(system))))
            newspin = Sunny._random_spin()

            func_diff = Sunny.local_energy_change(system, randsite, newspin)

            orig_energy = energy(system)
            system[randsite] = newspin
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

end