# TODO: Investigate TestItemRunner slowdown. Runtime is 1.1s on Sunny 0.7,
# mainly due to type inference. But the same code, compiled in a function from
# the terminal, is a small fraction of a second.
@testitem "Magnetization Observables" begin
    using LinearAlgebra

    positions = [[0.,0,0], [0.1, 0.2, 0.3]]
    cryst = Crystal(I(3), positions, 1)
    g1 = [0.8 3.2 5.1; -0.3 0.6 -0.1; 1.0 0. 0.]
    g2 = [1.8 2.2 -3.1; -0.3 2.6 0.1; 1.0 0. 0.3]

    for mode = [:SUN, :dipole]
        moments = [1 => Moment(s=3/2, g=g1), 2 => Moment(s=(mode == :SUN ? 3/2 : 1/2), g=g2)]
        sys = System(cryst, moments, mode)

        set_dipole!(sys, [1,3,2], (1,1,1,1))
        set_dipole!(sys, [3,4,5], (1,1,1,2))

        # Dipole magnetization observables (classical)
        sc = SampledCorrelationsStatic(sys; measure=ssf_custom((q, ssf) -> ssf, sys; apply_g=true))
        add_sample!(sc, sys)
        res = intensities_static(sc, [[0,0,0]])
        corr_mat = res.data[1]

        # Compute magnetization correlation "by hand", averaging over sites
        mag_corr = sum([sys.gs[i] * sys.dipoles[i] * (sys.gs[j] * sys.dipoles[j])' for i = 1:2, j = 1:2]) / Sunny.natoms(cryst)
        mag_corr_time_natoms = mag_corr * Sunny.natoms(cryst)
        @test isapprox(corr_mat, mag_corr_time_natoms)

        # For spin wave theory, check that `apply_g=true` is equivalent to
        # setting `apply_g` and manually contracting spin indices with the
        # g-tensor. This only works when g is equal among sites. TODO: Test with
        # anisotropic g-tensors.
        moments_homog = [1 => Moment(s=3/2, g=g1), 2 => Moment(s=(mode == :SUN ? 3/2 : 1/2), g=g1)]
        sys_homog = System(cryst, moments_homog, mode)

        measure = ssf_custom((q, ssf) -> ssf, sys_homog; apply_g=false)
        swt = SpinWaveTheory(sys_homog; measure)
        res1 = intensities_bands(swt, [[0,0,0]])

        measure = ssf_custom((q, ssf) -> ssf, sys_homog; apply_g=true)
        swt = SpinWaveTheory(sys_homog; measure)
        res2 = intensities_bands(swt, [[0,0,0]])

        @test isapprox(g1 * res1.data[1] * g1', res2.data[1])
        @test isapprox(g1 * res1.data[2] * g1', res2.data[2])
    end
end

@testitem "Available Energies Dirac Identity" begin
     # Create a dummy SampledCorrelations object
    cryst = Sunny.cubic_crystal()
    sys = System(cryst, [1 => Moment(s=1/2, g=2)], :SUN; seed=0)
    dt = 0.08
    sc = SampledCorrelations(sys; dt, energies=range(0.0, 10.0, 100), measure=ssf_perp(sys))

    ωs = Sunny.available_energies(sc; negative_energies=true)
    dts = 0:(sc.dt * sc.measperiod):3
    vals = sum(exp.(im .* ωs .* dts'), dims=1)[:]

    # Verify it made a delta function
    @test vals[1] ≈ length(ωs)
    @test all(abs.(vals[2:end]) .< 1e-12)
end

@testitem "Sum rule with reshaping" begin
    s = 1/2
    g = 2.3
    cryst = Sunny.diamond_crystal()
    sys = System(cryst, [1 => Moment(; s, g)], :SUN)
    randomize_spins!(sys)
    sc = SampledCorrelationsStatic(sys; measure=ssf_trace(sys; apply_g=true))
    add_sample!(sc, sys)

    # For the diamond cubic crystal, reciprocal space is periodic over a
    # distance of 4 BZs. Verify that the average intensity matches the expected
    # sum rule.
    counts = (4, 4, 4)
    qs = Sunny.available_wave_vectors(sc.parent; counts)
    res = intensities_static(sc, qs[:])
    @test sum(res.data) / length(qs) ≈ Sunny.natoms(cryst) * s^2 * g^2

    # Repeat the same calculation for a primitive cell.
    shape = [0 1 1; 1 0 1; 1 1 0] / 2
    sys_prim = reshape_supercell(sys, shape)
    sys_prim = repeat_periodically(sys_prim, (4, 4, 4))
    sc_prim = SampledCorrelationsStatic(sys_prim; measure=ssf_trace(sys_prim; apply_g=true))
    add_sample!(sc_prim, sys_prim)

    qs = Sunny.available_wave_vectors(sc_prim.parent; counts)
    res_prim = intensities_static(sc_prim, qs[:])
    @test sum(res_prim.data) / length(qs) ≈ Sunny.natoms(cryst) * s^2 * g^2
end

@testitem "Sublattice sum rule" begin
    latvecs = lattice_vectors(1, 1, 2, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0], [0, 0, 0.5]]; types=["A", "B"])
    s = 1/2
    g1 = 2
    g2 = 2.3
    sys = System(cryst, [1 => Moment(; s, g=g1), 2 => Moment(; s, g=g2)], :dipole)
    randomize_spins!(sys)

    sc = SampledCorrelationsStatic(sys; measure=ssf_trace(sys))
    qs = Sunny.available_wave_vectors(sc.parent; counts=(1, 1, 2))
    add_sample!(sc, sys)
    res = intensities_static(sc, qs[:])
    @test sum(res.data) / length(qs) ≈ s^2 * (g1^2 + g2^2)

    # Just atoms A
    formfactors = [1 => one(FormFactor), 2 => zero(FormFactor)]
    sc = SampledCorrelationsStatic(sys; measure=ssf_trace(sys; formfactors))
    add_sample!(sc, sys)
    res = intensities_static(sc, qs[:])
    @test sum(res.data) / length(qs) ≈ s^2 * g1^2

    # Just atoms B
    formfactors = [1 => zero(FormFactor), 2 => one(FormFactor)]
    sc = SampledCorrelationsStatic(sys; measure=ssf_trace(sys; formfactors))
    add_sample!(sc, sys)
    res = intensities_static(sc, qs[:])
    @test sum(res.data) / length(qs) ≈ s^2 * g2^2
end

@testitem "Dynamic sum rule" begin
    s = 3/2
    g = 2.3
    cryst = Sunny.cubic_crystal()
    sys = System(cryst, [1 => Moment(; s, g)], :dipole)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys)

    sc = SampledCorrelations(sys; dt=0.1, energies=range(0.0, 1.0, 3), measure=ssf_trace(sys))
    add_sample!(sc, sys)

    qs = Sunny.available_wave_vectors(sc)
    res = intensities(sc, qs[:]; energies=:available_with_negative, kT=nothing)

    # Integrate over energies to get static intensity, then average over sampled
    # q values to get sum rule.
    @test sum(res.data) * sc.Δω / length(qs) ≈ (s*g)^2
end
