@testitem "Dynamics" begin
    include("shared.jl")

    "Tests that SphericalMidpoint conserves energy for simple forces to a certain tolerance."
    function test_spherical_midpoint(; SUN)
        cryst = Sunny.diamond_crystal()
        ints = Sunny.AbstractInteraction[]
        add_linear_interactions!(ints, SUN)
        add_quadratic_interactions!(ints, SUN)
        add_quartic_interactions!(ints, SUN)
        sys = SpinSystem(cryst, ints, (3, 3, 3); SUN, seed=0)
        enable_dipole_dipole!(sys)

        randomize_spins!(sys)

        NITERS = 5_000
        Δt     = 0.005
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

    test_spherical_midpoint(; SUN=false)
    test_spherical_midpoint(; SUN=true)


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
