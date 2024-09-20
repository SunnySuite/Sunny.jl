@testitem "Kappa rescaling" begin
    include("shared.jl")

    cryst = Sunny.diamond_crystal()
    damping = 0.1
    dt = 0.005
    
    # Check that magnitude of coherent (SUN=true) or dipole (SUN=false) is
    # invariant under the dynamics
    let
        integrators = (Langevin(dt; damping, kT=0.1), ImplicitMidpoint(dt; damping))

        for integrator in integrators
            for mode in (:SUN, :dipole)
                sys = System(cryst, [1 => Moment(s=3/2, g=2)], mode; dims=(2, 2, 2), seed=0)
                randn!(sys.κs)
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


    # Check that each energy term rescales properly with κ
    let
        function gen_energy(κ, adder, mode)
            sys = System(cryst, [1 => Moment(s=3/2, g=2)], mode; dims=(2, 2, 2), seed=0)
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
        function gen_trajectory(κ, dt, adder, mode)
            sys = System(cryst, [1 => Moment(s=3/2, g=2)], mode; dims=(2, 2, 2), seed=0)
            adder(sys, mode)
            sys.κs .= κ
            randomize_spins!(sys)
            integrator = ImplicitMidpoint(dt)
            for _ in 1:100
                step!(sys, integrator)
            end
            return first(sys.dipoles)
        end
    
        κ = 2.0
        for mode in (:SUN, :dipole)
            s1 = gen_trajectory(1, dt, add_linear_interactions!, mode)
            s2 = gen_trajectory(κ, dt, add_linear_interactions!, mode)
            @test s1 ≈ s2/κ

            s1 = gen_trajectory(1, dt, add_quadratic_interactions!, mode)
            s2 = gen_trajectory(κ, dt/κ, add_quadratic_interactions!, mode)
            @test s1 ≈ s2/κ

            s1 = gen_trajectory(1, dt, add_quartic_interactions!, mode)
            s2 = gen_trajectory(κ, dt/κ^3, add_quartic_interactions!, mode)
            @test s1 ≈ s2/κ
        end
    end
end


@testitem "Anisotropy rescaling" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]], "P1")
    s = 3
    λ = Sunny.rcs_factors(s)
    
    for k in (2, 4, 6)
        c = randn(2k+1)
        E1, E2 = map([:dipole, :dipole_uncorrected]) do mode
            sys = System(cryst, [1 => Moment(; s, g=2)], mode)
            O = stevens_matrices(spin_label(sys, 1))
            set_onsite_coupling!(sys, sum(c[k-q+1]*O[k, q] for q in -k:k), 1)
            return energy(sys)
        end
        @test E1 ≈ λ[k] * E2
    end
end


@testitem "Biquadratic renormalization" begin
    cryst = Sunny.square_crystal()

    s = 3/2
    sys1 = System(cryst, [1 => Moment(; s, g=2)], :dipole_uncorrected, seed=0)
    sys2 = System(cryst, [1 => Moment(; s, g=2)], :dipole, seed=0)

    # Reference scalar biquadratic without renormalization
    set_exchange!(sys1, 0, Bond(1, 1, [1,0,0]); biquad=1)
    @test sys1.interactions_union[1].pair[1].bilin ≈ 0
    @test sys1.interactions_union[1].pair[1].biquad ≈ 1

    # Same thing, but with renormalization
    rcs = Sunny.rcs_factors(s)[2]^2
    set_exchange!(sys2, 0, Bond(1, 1, [1,0,0]); biquad=1)
    @test rcs ≈ (1-1/2s)^2
    @test sys2.interactions_union[1].pair[1].bilin ≈ -1/2
    @test sys2.interactions_union[1].pair[1].biquad ≈ 1 * rcs

    # Same thing, but with explicit operators. Internally, Sunny decomposes
    # (S1'*S2)^2 into two parts, and only the first gets rescaled by the `rcs`
    # factor:
    #   1. (S1'*S2)^2 + S1'*S2/2 (a pure coupling of Stevens quadrupoles)
    #   2. -S1'*S2/2             (a Heisenberg shift)
    set_pair_coupling!(sys2, (S1, S2) -> (S1'*S2)^2, Bond(1, 1, [1,0,0]))
    @test sys2.interactions_union[1].pair[1].bilin ≈ -1/2
    @test sys2.interactions_union[1].pair[1].biquad ≈ 1 * rcs
end

@testitem "Biquadratic renormalization 2" begin
    # Simple dimer model
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0], [0.3, 0, 0]]; types=["A", "B"])
    s1 = 3/2
    s2 = 2
    v1 = randn(3)
    v2 = randn(3)
    biquad = 1.2
    bond = Bond(1, 2, [0, 0, 0])

    # Biquadratic energy in dipole mode with RCS
    sys = System(cryst, [1 => Moment(s=s1, g=-1), 2 => Moment(s=s1, g=-1)], :dipole)
    set_dipole!(sys, v1, (1, 1, 1, 1))
    set_dipole!(sys, v2, (1, 1, 1, 2))
    set_exchange!(sys, 0.0, bond; biquad)
    E_dipole = energy(sys)

    # Biquadratic energy in SU(N) mode is same
    sys = System(cryst, [1 => Moment(s=s1, g=-1), 2 => Moment(s=s1, g=-1)], :SUN)
    set_dipole!(sys, v1, (1, 1, 1, 1))
    set_dipole!(sys, v2, (1, 1, 1, 2))
    set_exchange!(sys, 0.0, bond; biquad)
    E_SUN_1 = energy(sys)
    set_pair_coupling!(sys, (Si, Sj) -> biquad * (Si'*Sj)^2, bond)
    E_SUN_2 = energy(sys)
    @test E_dipole ≈ E_SUN_1 ≈ E_SUN_2

    # Biquadratic energy in large-s limit is classical formula
    sys = System(cryst, [1 => Moment(s=s1, g=-1), 2 => Moment(s=s2, g=-1)], :dipole_uncorrected)
    set_dipole!(sys, v1, (1, 1, 1, 1))
    set_dipole!(sys, v2, (1, 1, 1, 2))
    set_exchange!(sys, 0.0, bond; biquad)
    E_large_s = energy(sys)
    @test E_large_s ≈ biquad * (sys.dipoles[1]' * sys.dipoles[2])^2
end
