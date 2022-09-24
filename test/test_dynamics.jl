@testitem "Dynamics" begin
include("test_shared.jl")

"Tests that SphericalMidpoint conserves energy for simple forces to a certain tolerance."
function test_spherical_midpoint()
    crystal = Sunny.diamond_crystal()

    interactions = [
        diamond_test_exchanges()...,
        external_field([0, 0, 1])
    ]

    sys = SpinSystem(crystal, interactions, (5, 5, 5))
    Random.seed!(0)
    rand!(sys)

    NITERS = 5_000
    Δt     = 0.01
    energies = Float64[]

    integrator = Sunny.SphericalMidpoint(sys)
    for it in 1:NITERS
        evolve!(integrator, Δt)
        # Compute the energy
        push!(energies, energy(sys))
    end

    # Check that the energy hasn't fluctuated much. Expected fluctuations should
    # scale like square-root of system size.
    sqrt_size = sqrt(length(sys._dipoles))
    ΔE = (maximum(energies) - minimum(energies)) / sqrt_size
    @test ΔE < 1e-2
end

test_spherical_midpoint()

"Tests that SphericalMidpoint conserves energy for dipole forces (w/ Fourier) to a certain tolerance."
function test_dipole_ft()
    crystal = Sunny.diamond_crystal()
    interactions = [dipole_dipole()]
    sys = SpinSystem(crystal, interactions, (5, 5, 5))
    Random.seed!(0)
    rand!(sys)

    NITERS = 5_000
    Δt     = 0.01
    energies = Float64[]

    integrator = Sunny.SphericalMidpoint(sys)
    for it in 1:NITERS
        evolve!(integrator, Δt)
        # Compute the energy
        push!(energies, energy(sys))
    end

    # Check that the energy hasn't fluctuated much
    sqrt_size = sqrt(length(sys._dipoles))
    ΔE = (maximum(energies) - minimum(energies)) / sqrt_size
    @test ΔE < 2e-2
end

test_dipole_ft()

"Tests that set_temp!/get_temp behave as expected"
function test_set_get_temp_langevin()
    system = produce_example_system()
    
    Random.seed!(1111)
    rand_kTs = rand(100)

    sampler = LangevinSampler(system, 1.0, 1.0, 1.0, 1)
    for kT in rand_kTs
        set_temp!(sampler, kT)
        if get_temp(sampler) != kT
            return false
        end
    end

    return true
end

@test test_set_get_temp_langevin()

end

# == These should not be @test-ed, but are for manual inspection == #

"Measure timings for speed disparity between real and fourier-space dipole interactions"
function time_real_fourier_dipole()
    lat_vecs = [1 0 0; 0 1 0; 0 0 1]
    positions = [[0, 0, 0], [1, 1, 1]/2]
    crystal = Crystal(lat_vecs, positions)

    Ls = [2, 6, 10, 14, 18, 22, 26, 30]
    real_times = Float64[]
    ft_times = Float64[]
    for L in Ls
        println("Measuring L = $L")

        # We'll manually set interactions, since DipoleRealCPU
        #  isn't actually a create-able force in the current
        #  SpinSystem/HamiltonianCPU constructor.
        interactions = Sunny.Interaction[]
        latsize = (L, L, L)
        sys = SpinSystem(crystal, interactions, latsize)
        Random.seed!(0)
        rand!(sys)

        # Since we just want to time integration speed,
        #  it's alright for the Ewald summation to be
        #  very inaccurate (truncated at a single cell)
        dip_dip = dipole_dipole(;extent=1)
        empty_h = Sunny.HeisenbergCPU{3}[]
        empty_d = Sunny.DiagonalCouplingCPU{3}[]
        empty_g = Sunny.GeneralCouplingCPU{3}[]

        if L <= 14
            println("\tSetting up real-space...")
            dipr = Sunny.DipoleRealCPU(dip_dip, crystal, latsize)
            ℋ = Sunny.HamiltonianCPU(
                nothing, empty_h, empty_d, empty_g, dipr
            )
            sys.hamiltonian = ℋ
            integrator = HeunP(sys)

            println("\tStarting real-space...")
            evolve!(integrator, 0.001)
            _, t, _, _ = @timed for _ in 1:100
                evolve!(integrator, 0.001)
            end
            push!(real_times, t)
        end

        println("\tSetting up Fourier-space...")
        dipf = Sunny.DipoleFourierCPU(dip_dip, crystal, latsize)
        ℋ = Sunny.HamiltonianCPU(
            nothing, empty_h, empty_d, empty_g, dipf
        )
        sys.hamiltonian = ℋ
        integrator = HeunP(sys)

        println("\tStarting Fourier-space...")
        evolve!(integrator, 0.001)
        _, t, _, _ = @timed for _ in 1:100
            evolve!(integrator, 0.001)
        end
        push!(ft_times, t)
    end
    return real_times, ft_times
end