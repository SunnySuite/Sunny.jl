@testitem "Kappa rescaling" begin
    include("shared.jl")

    # Check that magnitude of coherent (SUN=true) or dipole (SUN=false) is
    # invariant under the dynamics
    let
        cryst = Sunny.diamond_crystal()
        kT = 0.1
        λ = 0.1
        Δt = 0.01
        integrators = (Langevin(Δt; kT, λ), ImplicitMidpoint(Δt))

        for integrator in integrators
            for mode in (:SUN, :dipole)
                sys = System(cryst, (3, 3, 3), [SpinInfo(1; S=5 / 2, g=2)], mode; seed=0)
                randn!(sys.κs)
                add_linear_interactions!(sys, mode)
                add_quadratic_interactions!(sys, mode)
                add_quartic_interactions!(sys, mode)
                randomize_spins!(sys)
                mags1 = norm.(mode == :SUN ? sys.coherents : sys.dipoles)
                for _ in 1:100
                    step!(sys, integrator)
                end
                mags2 = norm.(mode == :SUN ? sys.coherents : sys.dipoles)
                @test mags1 ≈ mags2
            end
        end
    end

    # Check that each energy term rescales properly with κ
    let
        function gen_energy(κ, adder, mode)
            cryst = Sunny.diamond_crystal()
            sys = System(cryst, (2, 2, 2), [SpinInfo(1; S=5 / 2, g=2)], mode; seed=0)
            # κ must be set before anisotropy operators are added, otherwise we
            # will lose information about how to rescale the Stevens expansion.
            sys.κs .= κ
            adder(sys, mode)
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

    # Check that a scaling of κ corresponds to an appropriate rescaling of dynamical time
    let
        function gen_trajectory(κ, Δt, adder, mode)
            cryst = Sunny.diamond_crystal()
            sys = System(cryst, (4, 3, 2), [SpinInfo(1; S=5 / 2, g=2)], mode; seed=0)
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
            @test s1 ≈ s2 / κ

            s1 = gen_trajectory(1, Δt, add_quadratic_interactions!, mode)
            s2 = gen_trajectory(κ, Δt / κ, add_quadratic_interactions!, mode)
            @test s1 ≈ s2 / κ

            s1 = gen_trajectory(1, Δt, add_quartic_interactions!, mode)
            s2 = gen_trajectory(κ, Δt / κ^3, add_quartic_interactions!, mode)
            @test s1 ≈ s2 / κ
        end
    end
end

@testitem "Anisotropy rescaling" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]], "P1")
    S = 3
    λ = Sunny.anisotropy_renormalization(S)

    for k in (2, 4, 6)
        c = randn(2k + 1)
        E1, E2 = map([:dipole, :dipole_large_S]) do mode
            sys = System(cryst, (1, 1, 1), [SpinInfo(1; S, g=2)], mode)
            O = stevens_matrices(spin_label(sys, 1))
            set_onsite_coupling!(sys, sum(c[k-q+1] * O[k, q] for q in -k:k), 1)
            return energy(sys)
        end
        @test E1 ≈ λ[k] * E2
    end
end
