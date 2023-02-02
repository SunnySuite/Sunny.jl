@testitem "Spin Scaling" begin
    include("shared.jl")

    # Check that magnitude of coherent (SUN=true) or dipole (SUN=false) is
    # invariant under the dynamics
    function test_spin_magnitude_stability()
        cryst = Sunny.diamond_crystal()
        kT = 0.1
        λ  = 0.1
        Δt = 0.01
        integrators = (LangevinHeunP(kT, λ, Δt), ImplicitMidpoint(Δt))

        for integrator in integrators
            for mode in (:SUN, :dipole)
                sys = System(cryst, (3,3,3), [SpinInfo(1, S=5/2)], mode; seed=0)
                add_linear_interactions!(sys, mode)
                add_quadratic_interactions!(sys, mode)
                add_quartic_interactions!(sys, mode)
                randomize_spins!(sys)
                mags1 = norm.(mode==:SUN ? sys.coherents : sys.dipoles)
                for _ in 1:100
                    step!(sys, integrator)
                end
                mags2 = norm.(mode==:SUN ? sys.coherents : sys.dipoles)
                @test mags1 ≈ mags2
            end
        end
    end

    test_spin_magnitude_stability()


    # Check that each energy term rescales properly with κ
    function test_energy_scaling()
        function gen_energy(κ, adder, mode)
            cryst = Sunny.diamond_crystal()
            sys = System(cryst, (2,2,2), [SpinInfo(1, S=5/2)], mode; seed=0)
            adder(sys, mode)
            sys.κs .= κ
            randomize_spins!(sys)
            return energy(sys)
        end

        κ = 2.0
        for mode in (:SUN, :dipole)
            E1 = gen_energy(1, add_linear_interactions!, mode)
            E2 = gen_energy(κ, add_linear_interactions!, mode)
            @test E1 ≈ E2 / κ

            E1 = gen_energy(1, add_quadratic_interactions!, mode)
            E2 = gen_energy(κ, add_quadratic_interactions!, mode)
            @test E1 ≈ E2 / κ^2

            E1 = gen_energy(1, add_quartic_interactions!, mode)
            E2 = gen_energy(κ, add_quartic_interactions!, mode)
            @test E1 ≈ E2 / κ^4
        end
    end

    test_energy_scaling()


    # Check that a scaling of κ corresponds to an appropriate rescaling of dynamical time
    # TODO: Figure out scaling for Langevin dynamics?
    function test_dynamics_scaling()
        function gen_trajectory(κ, Δt, adder, mode)
            cryst = Sunny.diamond_crystal()
            sys = System(cryst, (4,3,2), [SpinInfo(1, S=5/2)], mode; seed=0)
            adder(sys, mode)
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
        for mode in (:SUN, :dipole)
            s1 = gen_trajectory(1, Δt, add_linear_interactions!, mode)
            s2 = gen_trajectory(κ, Δt, add_linear_interactions!, mode)
            @test s1 ≈ s2/κ

            s1 = gen_trajectory(1, Δt, add_quadratic_interactions!, mode)
            s2 = gen_trajectory(κ, Δt/κ, add_quadratic_interactions!, mode)
            @test s1 ≈ s2/κ

            s1 = gen_trajectory(1, Δt, add_quartic_interactions!, mode)
            s2 = gen_trajectory(κ, Δt/κ^3, add_quartic_interactions!, mode)
            @test s1 ≈ s2/κ
        end
    end

    test_dynamics_scaling()

end
