# These tests checks that varying μ0, μB results in the correct scale
# transformations in various energies / fields. Note:
# does not check that the actual absolute values are
# correct.

seed = 1111
μBs = [1, 2, 1, 4]
μ0s = [1, 1, 2, 3]
crystal = Sunny.diamond_crystal()
latsize = (4, 4, 4)


function collect_energy_and_field(int, μB, μ0)
    sys = SpinSystem(crystal, [int], latsize; μB=μB, μ0=μ0)
    Random.seed!(seed)
    rand!(sys)
    return energy(sys), field(sys)
end

# All exchange interactions should be invariant under μ0, μB changes
function test_exchanges_scaling(μBs, μ0s)
    exchange_ints = diamond_test_exchanges()
    
    # Check that for all values of μB, μ0, all energies
    #  and fields are unchanged compared to first set of values.
    ref_μB, ref_μ0 = μBs[1], μ0s[1]
    @testset "Exchange Scaling" begin
        for int in exchange_ints
            (ref_E, ref_B) = collect_energy_and_field(int, ref_μB, ref_μ0)
            for (μB, μ0) in zip(μBs[2:end], μ0s[2:end])
                (E, B) = collect_energy_and_field(int, μB, μ0)
                @test E ≈ ref_E
                @test B ≈ ref_B
            end
        end
    end
end

# Zeeman energies/fields should scale linearly with μB, invariant to μ0
function test_field_scaling(μBs, μ0s)
    ref_μB, ref_μ0 = μBs[1], μ0s[1]

    ext_field = external_field(rand(3))

    (ref_E, ref_B) = collect_energy_and_field(ext_field, ref_μB, ref_μ0)
    @testset "Zeeman Scaling" begin
        for (μB, μ0) in zip(μBs[2:end], μ0s[2:end])
            (E, B) = collect_energy_and_field(ext_field, μB, μ0)
            @test E ≈ (μB / ref_μB) * ref_E
            @test B ≈ (μB / ref_μB) .* ref_B
        end
    end
end

# Dipole-dipole interactions should scale linearly with μ0, and
#  quadratically with μB
function test_dipole_scaling(μBs, μ0s)
    ref_μB, ref_μ0 = μBs[1], μ0s[1]

    dip_dip = dipole_dipole()
    
    (ref_E, ref_B) = collect_energy_and_field(dip_dip, ref_μB, ref_μ0)
    @testset "Dipole-Dipole Scaling" begin
        for (μB, μ0) in zip(μBs[2:end], μ0s[2:end])
            (E, B) = collect_energy_and_field(dip_dip, μB, μ0)
            @test E ≈ (μ0 / ref_μ0) * (μB / ref_μB)^2 * ref_E
            @test B ≈ (μ0 / ref_μ0) * (μB / ref_μB)^2 .* ref_B
        end
    end
end


@testset verbose=true "Unit Scaling" begin
    test_exchanges_scaling(μBs, μ0s)
    test_field_scaling(μBs, μ0s)
    test_dipole_scaling(μBs, μ0s)
end