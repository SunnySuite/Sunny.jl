@testitem "Dimensional Scaling" begin
    include("shared.jl")

    # These tests checks that varying μ0, μB results in the correct scale
    # transformations in various energies / fields. Note:
    # does not check that the actual absolute values are
    # correct.

    crystal = Sunny.diamond_crystal()
    latsize = (4, 4, 4)

    units = [Sunny.Units.meV, Sunny.Units.theory]

    function collect_energy_and_forces(test, units)
        sys = SpinSystem(crystal, latsize; units, seed=0)
        if test == :zeeman
            set_external_field!(sys, randn(sys.rng, Sunny.Vec3))
        elseif test == :exchange
            add_exchange_interactions!(sys, false)
        else test == :dipdip
            enable_dipole_dipole!(sys)
        end
        randomize_spins!(sys)
        return energy(sys), forces(sys)
    end

    # All exchange interactions should be invariant under changes to physical
    # constants
    function validate_exchanges_scaling()
        E1, B1 = collect_energy_and_forces(:exchange, units[1])
        E2, B2 = collect_energy_and_forces(:exchange, units[2])
        @test E1 ≈ E2
        @test B1 ≈ B2
    end

    validate_exchanges_scaling()


    # Zeeman energies/forces should scale linearly with μB, invariant to μ0
    function validate_zeeman_scaling()
        E1, B1 = collect_energy_and_forces(:zeeman, units[1])
        E2, B2 = collect_energy_and_forces(:zeeman, units[2])

        @test E1 / units[1].μB ≈ E2 / units[2].μB
        @test B1 / units[1].μB ≈ B2 / units[2].μB
    end

    validate_zeeman_scaling()


    # Dipole-dipole interactions should scale linearly with μ0, and
    # quadratically with μB
    function validate_dipole_scaling()
        E1, B1 = collect_energy_and_forces(:dipdip, units[1])
        E2, B2 = collect_energy_and_forces(:dipdip, units[2])

        (; μB, μ0) = units[1]
        c1 = μ0*μB^2
        (; μB, μ0) = units[2]
        c2 = μ0*μB^2
        @test E1 / c1 ≈ E2 / c2
        @test B1 / c1 ≈ B2 / c2
    end

    validate_dipole_scaling()

end
