using StaticArrays

function test_local_energy_change()
    lat_vecs = lattice_vectors(1.0, 1.0, 2.0, 90., 90., 120.)
    basis_vecs = [SA[0.0, 0.0, 0.0], SA[0.5, 0.5, 0.5]]
    latsize = SA[5, 5, 5]
    lattice = Lattice(lat_vecs, basis_vecs, latsize)
    cryst = Crystal(lattice)

    # Make one of each kind of interaction with some arbitrary numbers
    B = ExternalField(SA[-2.4, -2.4, 0.])
    J1 = Heisenberg(1.5, cryst, Bond{3}(1, 2, [0, 0, 0]))
    J2 = DiagonalCoupling(SA[1.2, 1.2, -1.4], cryst, Bond{3}(1, 1, [1, 0, 0]))
    J3mat = SA[0.5 0.0 0.0; 0.0 0.7 1.3; 0.0 -1.3 2.5]
    J3 = GeneralCoupling(J3mat, cryst, Bond{3}(1, 1, [0, 0, 1]))
    D = OnSite(SA[0.0, 0.0, -2.3])
    ℋ = Hamiltonian{3}([B, J1, J2, J3, D])

    system = SpinSystem(lattice, ℋ)

    for _ in 1:3
        rand!(system)
        for _ in 1:50
            # Pick a random site, try to set it to a random spin
            randsite = CartesianIndex(Tuple(mod1.(rand(SVector{4, Int}), size(system))))
            newspin = FastDipole._random_spin()

            func_diff = FastDipole.local_energy_change(system, randsite, newspin)

            orig_energy = energy(system)
            system[randsite] = newspin
            new_energy = energy(system)

            actual_diff = new_energy - orig_energy

            @assert func_diff ≈ actual_diff
        end
    end
end