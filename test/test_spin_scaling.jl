@testitem "Spin Scaling" begin
    include("shared.jl")

    # Check that magnitude of coherent (SUN=true) or dipole (SUN=false) is
    # invariant under the dynamics
    function test_spin_magnitude_stability()
        cryst = Sunny.fcc_crystal()
        kT = 0.1
        λ  = 0.1
        Δt = 0.01
        integrators = (LangevinHeunP(kT, λ, Δt), ImplicitMidpoint(Δt))

        for integrator in integrators
            for SUN in (true, false)
                ints = Sunny.AbstractInteraction[]
                add_linear_interactions!(ints, SUN)
                add_quadratic_interactions!(ints, SUN)
                add_quartic_interactions!(ints, SUN)
                sys = SpinSystem(cryst, ints, (3,3,3), [SiteInfo(1; S=5/2)]; SUN, seed=0)
                randomize_spins!(sys)
                mags1 = norm.(SUN ? sys.coherents : sys.dipoles)
                for _ in 1:100
                    step!(sys, integrator)
                end
                mags2 = norm.(SUN ? sys.coherents : sys.dipoles)
                @test mags1 ≈ mags2
            end
        end
    end

    test_spin_magnitude_stability()


    # Check that each energy term rescales properly with κ
    function test_energy_scaling()
        function gen_energy(κ, int_adder, SUN)
            cryst = Sunny.diamond_crystal()
            ints = Sunny.AbstractInteraction[]
            int_adder(ints, SUN)
            sys = SpinSystem(cryst, ints, (2,2,2), [SiteInfo(1; S=5/2)]; SUN, seed=0)
            sys.κs .= κ
            randomize_spins!(sys)
            return energy(sys)
        end

        κ = 2.0
        for SUN in (true, false)
            E1 = gen_energy(1, add_linear_interactions!, SUN)
            E2 = gen_energy(κ, add_linear_interactions!, SUN)
            @test E1 ≈ E2 / κ

            E1 = gen_energy(1, add_quadratic_interactions!, SUN)
            E2 = gen_energy(κ, add_quadratic_interactions!, SUN)
            @test E1 ≈ E2 / κ^2

            E1 = gen_energy(1, add_quartic_interactions!, SUN)
            E2 = gen_energy(κ, add_quartic_interactions!, SUN)
            @test E1 ≈ E2 / κ^4
        end
    end

    test_energy_scaling()


    # Check that a scaling of κ corresponds to an appropriate rescaling of dynamical time
    # TODO: Figure out scaling for Langevin dynamics?
    function test_dynamics_scaling()
        function gen_trajectory(κ, Δt, int_adder, SUN)
            cryst = Sunny.diamond_crystal()
            ints = Sunny.AbstractInteraction[]
            int_adder(ints, SUN)
            sys = SpinSystem(cryst, ints, (4,3,2), [SiteInfo(1; S=5/2)]; SUN, seed=0)
            sys.κs .= κ
            randomize_spins!(sys)
            integrator = ImplicitMidpoint(Δt)
            for _ in 1:100
                step!(sys, integrator)
            end
            return first(sys.dipoles)
        end
    
        κ = 2.0
        Δt = 0.005
        for SUN in (true, false)
            s1 = gen_trajectory(1, Δt, add_linear_interactions!, SUN)
            s2 = gen_trajectory(κ, Δt, add_linear_interactions!, SUN)
            @test s1 ≈ s2/κ

            s1 = gen_trajectory(1, Δt, add_quadratic_interactions!, SUN)
            s2 = gen_trajectory(κ, Δt/κ, add_quadratic_interactions!, SUN)
            @test s1 ≈ s2/κ

            s1 = gen_trajectory(1, Δt, add_quartic_interactions!, SUN)
            s2 = gen_trajectory(κ, Δt/κ^3, add_quartic_interactions!, SUN)
            @test s1 ≈ s2/κ
        end
    end

    test_dynamics_scaling()

end
