@testitem "Dynamics" begin
include("test_shared.jl")

"Tests that SphericalMidpoint conserves energy for simple forces to a certain tolerance."
function test_spherical_midpoint()
    crystal = Sunny.diamond_crystal()

    interactions = [
        diamond_test_exchanges()...,
        external_field([0, 0, 1])
    ]

    sys = SpinSystem(crystal, interactions, (5, 5, 5); seed=0)
    rand!(sys)

    NITERS = 5_000
    Δt     = 0.01
    energies = Float64[]

    integrator = ImplicitMidpoint(Δt)
    for _ in 1:NITERS
        step!(sys, integrator)
        push!(energies, energy(sys))
    end

    # Check that the energy hasn't fluctuated much. Expected fluctuations should
    # scale like square-root of system size.
    sqrt_size = sqrt(length(sys.dipoles))
    ΔE = (maximum(energies) - minimum(energies)) / sqrt_size
    @test ΔE < 1e-2
end

test_spherical_midpoint()

"Tests that SphericalMidpoint conserves energy for dipole forces (w/ Fourier) to a certain tolerance."
function test_dipole_ft()
    crystal = Sunny.diamond_crystal()
    interactions = Sunny.AbstractInteraction[]
    sys = SpinSystem(crystal, interactions, (5, 5, 5); seed=0)
    rand!(sys)
    enable_dipole_dipole!(sys)

    NITERS = 5_000
    Δt     = 0.01
    energies = Float64[]

    # TODO: Benchmark and compare with previous versions
    integrator = ImplicitMidpoint(Δt)
    for _ in 1:NITERS
        step!(sys, integrator)
        push!(energies, energy(sys))
    end

    # Check that the energy hasn't fluctuated much
    sqrt_size = sqrt(length(sys.dipoles))
    ΔE = (maximum(energies) - minimum(energies)) / sqrt_size
    @test ΔE < 2e-2
end

test_dipole_ft()

"Tests that set_temp!/get_temp behave as expected"
function test_set_get_temp_langevin()
    integrator = LangevinHeunP(1.0, 1.0, 1.0)
    sampler = LangevinSampler(integrator, 1)
    for kT in [1.1, 10.3]
        set_temp!(sampler, kT)
        if get_temp(sampler) != kT
            return false
        end
    end

    return true
end

@test test_set_get_temp_langevin()

end
