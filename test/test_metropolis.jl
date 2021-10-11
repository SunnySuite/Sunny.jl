function test_local_energy_change()
    lat_vecs = lattice_vectors(1.0, 1.0, 2.0, 90., 90., 120.)
    basis_vecs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    latsize = [5, 5, 5]
    cryst = Crystal(lat_vecs, basis_vecs)

    # Make one of each kind of interaction with some arbitrary numbers
    B = ExternalField([-2.4, -2.4, 0.])
    J1 = heisenberg(1.5, Bond(1, 2, [0, 0, 0]))
    J2 = exchange(diagm([1.2, 2.6, -1.4]), Bond(1, 1, [1, 0, 0]))
    J3mat = [0.5 -1.2 0.0; -1.2 0.7 0.0; 0.0 0.0 2.5]
    J3 = exchange(J3mat, Bond(1, 1, [0, 1, 0]))
    D = onsite_anisotropy([0.0, 0.0, -2.3], 1)

    interactions = [B, J1, J2, J3, D]

    system = SpinSystem(cryst, interactions, latsize)

    for _ in 1:3
        rand!(system)
        for _ in 1:50
            # Pick a random site, try to set it to a random spin
            randsite = CartesianIndex(Tuple(mod1.(rand(SVector{4, Int}), size(system))))
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