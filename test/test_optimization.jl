@testitem "FM optimization" begin
    cryst = Crystal(lattice_vectors(1, 1, 2, 90, 90, 90), [[0, 0, 0]])

    for mode in (:dipole, :SUN)
        # H = -∑Sᵢ⋅Sⱼ - ∑(Sᵢᶻ)² on 2D square lattice (z-polarized ground state)
        sys = System(cryst, [1 => Moment(; s=3/2, g=2)], mode)
        reshape_supercell(sys, [[1, 0, 0] [1, 1, 0] [0, 0, 1]])
        set_exchange!(sys, -1, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> -S[3]^2, 1)

        # From initially x-polarized state, optimization without perturbation
        # can't break symmetry.
        polarize_spins!(sys, [1, 0, 0])
        @test minimize_energy!(sys; jitter=0).converged
        @test all(sys.dipoles) do S
            x = (mode == :dipole) ? 3/2 : 1.487610951765511
            S ≈ [x, 0, 0]
        end

        # From random initial condition, optimization reliably reaches ±ẑ
        # ground state
        randomize_spins!(sys)
        @test minimize_energy!(sys).converged
        @test all(sys.dipoles) do S
            abs.(S) ≈ [0, 0, 3/2]
        end
        @test energy_per_site(sys) ≈ -6.75
    end
end


@testitem "Canted AFM optimization" begin
    cryst = Crystal(lattice_vectors(1, 1, 2, 90, 90, 90), [[0, 0, 0]])

    for mode in (:dipole, :SUN)
        sys = System(cryst, [1 => Moment(; s=3/2, g=2)], mode; dims=(4, 4, 1))
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        set_field!(sys, [1.0, 0, 0])
        enable_dipole_dipole!(sys, 1.0)

        # Unstable canted FM phase. No (π, π) symmetry breaking. In SU(N) mode,
        # some weight shifts from dipole to quadrupole sector.
        polarize_spins!(sys, [0, 0, 1])
        @test minimize_energy!(sys; jitter=0).converged
        ref = (mode == :dipole) ? 0.6327349993784377 : -0.619368052000758
        @test energy_per_site(sys) ≈ ref

        # Stable canted AFM phase. Small perturbation δ breaks (π, π) symmetry.
        polarize_spins!(sys, [0, 0, 1])
        @test minimize_energy!(sys; jitter=1e-8).converged
        @test energy_per_site(sys) ≈ -5.7060439873282
    end
end


@testitem "Vacancies" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(3, 3, 1), seed=1)
    set_exchange!(sys, 1, Bond(1, 1, [1, 0, 0]))
    sys2 = to_inhomogeneous(sys)
    sys2.κs[1, 1, 1, 1] = 0.01
    randomize_spins!(sys2)
    minimize_energy!(sys2)
    energy_per_site(sys2)
    @test energy_per_site(sys2) ≈ -0.8889

    set_vacancy_at!(sys2, (1, 1, 1, 1))
    randomize_spins!(sys2)
    minimize_energy!(sys2)
    @test energy_per_site(sys2) ≈ -8/9
end
