@testitem "Unit Scaling" begin
include("test_shared.jl")

# These tests checks that varying μ0, μB results in the correct scale
# transformations in various energies / fields. Note:
# does not check that the actual absolute values are
# correct.

crystal = Sunny.diamond_crystal()
latsize = (4, 4, 4)


consts = [Sunny.CONSTS_meV, Sunny.CONSTS_ONES]

function collect_energy_and_field(ints, dipole_dipole, consts)
    sys = SpinSystem(crystal, ints, latsize; seed=1111, consts)
    dipole_dipole && enable_dipole_dipole!(sys)
    rand!(sys)
    return energy(sys), field(sys)
end

# All exchange interactions should be invariant under changes to physical
# constants
function validate_exchanges_scaling()
    ints = diamond_test_exchanges()
    E1, B1 = collect_energy_and_field(ints, false, consts[1])
    E2, B2 = collect_energy_and_field(ints, false, consts[2])
    @test E1 ≈ E2
    @test B1 ≈ B2
end

validate_exchanges_scaling()


# Zeeman energies/fields should scale linearly with μB, invariant to μ0
function validate_field_scaling()
    ext_field = external_field(rand(3))

    E1, B1 = collect_energy_and_field([ext_field], false, consts[1])
    E2, B2 = collect_energy_and_field([ext_field], false, consts[2])

    @test E1 / consts[1].μB ≈ E2 / consts[2].μB
    @test B1 / consts[1].μB ≈ B2 / consts[2].μB
end

validate_field_scaling()


# Dipole-dipole interactions should scale linearly with μ0, and
#  quadratically with μB
function validate_dipole_scaling()
    E1, B1 = collect_energy_and_field(Sunny.AbstractInteraction[], true, consts[1])
    E2, B2 = collect_energy_and_field(Sunny.AbstractInteraction[], true, consts[2])

    (; μB, μ0) = consts[1]
    c1 = μ0*μB^2
    (; μB, μ0) = consts[2]
    c2 = μ0*μB^2
    @test E1 / c1 ≈ E2 / c2
    @test B1 / c1 ≈ B2 / c2
end

validate_dipole_scaling()

end
