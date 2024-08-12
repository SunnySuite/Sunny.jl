@testitem "Magnetization Observables" begin
    using LinearAlgebra

    positions = [[0.,0,0], [0.1, 0.2, 0.3]]
    cryst = Crystal(I(3), positions, 1)
    g1 = [0.8 3.2 5.1; -0.3 0.6 -0.1; 1.0 0. 0.]
    g2 = [1.8 2.2 -3.1; -0.3 2.6 0.1; 1.0 0. 0.3]    

    for mode = [:SUN, :dipole]
        infos = [SpinInfo(1, S=3/2, g=g1), SpinInfo(2, S=(mode == :SUN ? 3/2 : 1/2), g=g2)]
        sys = System(cryst, (1,1,1), infos, mode)

        set_dipole!(sys, [1,3,2], (1,1,1,1))
        set_dipole!(sys, [3,4,5], (1,1,1,2))

        # Dipole magnetization observables (classical)
        sc = SampledCorrelations(sys; measure=ssf_custom((q, ssf) -> ssf, sys; apply_g=true), energies=nothing)
        add_sample!(sc, sys)
        is = intensities_instant(sc, Sunny.QPoints([[0,0,0]]); kT=nothing)
        corr_mat = is.data[1]

        # Compute magnetization correlation "by hand", averaging over sites
        mag_corr = sum([sys.gs[i] * sys.dipoles[i] * (sys.gs[j] * sys.dipoles[j])' for i = 1:2, j = 1:2]) / Sunny.natoms(cryst)
        mag_corr_time_natoms = mag_corr * Sunny.natoms(cryst)
        @test isapprox(corr_mat, mag_corr_time_natoms)

        # For spin wave theory, check that `apply_g=true` is equivalent to
        # setting `apply_g` and manually contracting spin indices with the
        # g-tensor. This only works when g is homogeneous. TODO: Test with
        # inhomogeneous g-tensors.
        infos_homog = [SpinInfo(1, S=3/2, g=g1), SpinInfo(2, S=(mode == :SUN ? 3/2 : 1/2), g=g1)]
        sys_homog = System(cryst, (1,1,1), infos_homog, mode)

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
    latsize = (1,1,1)
    cryst = Sunny.cubic_crystal()
    sys = System(cryst, latsize, [SpinInfo(1; S = 1/2, g=2)], :SUN; seed = 0)
    dt = 0.1
    sc = SampledCorrelations(sys; dt, energies=range(0.0, 10.0, 100), measure=ssf_perp(sys))

    ωs = available_energies(sc;negative_energies=true)
    dts = 0:(sc.dt * sc.measperiod):3
    vals = sum(exp.(im .* ωs .* dts'),dims = 1)[:]

    # Verify it made a delta function
    @test vals[1] ≈ length(ωs)
    @test all(isapprox.(0,vals[2:end];atol = 1e-12))
end

@testitem "Polyatomic sum rule" begin
    sys = System(Sunny.diamond_crystal(),(4,1,1),[SpinInfo(1,S=1/2,g=2)],:SUN,seed=1)
    randomize_spins!(sys)
    sc = SampledCorrelations(sys; dt=1.0, energies=range(0.0, 1.0, 3), measure=ssf_trace(sys; apply_g=true))
    add_sample!(sc, sys)

    sum_rule_ixs = [1, 4, 6]  # indices for zz, yy, xx
    sub_lat_sum_rules = sum(sc.data[sum_rule_ixs,:,:,:,:,:,:], dims=[1,4,5,6,7])[1,:,:,1,1,1,1]

    Δq³ = 1/prod(sys.latsize) # Fraction of a BZ
    n_all_ω = size(sc.data, 7)
    # Intensities in sc.data are a density in q, but already integrated over dω
    # bins, and then scaled by n_all_ω. Therefore, we need the factor below to
    # convert the previous sum to an integral. (See same logic in
    # intensities_interpolated function.)
    sub_lat_sum_rules .*= Δq³ / n_all_ω

    # SU(N) sum rule for S = 1/2:
    # ⟨∑ᵢSᵢ²⟩ = 3/4 on every site, but because we're classical, we
    # instead compute ∑ᵢ⟨Sᵢ⟩² = (1/2)^2 = 1/4 since the ⟨Sᵢ⟩ form a vector with
    # length (1/2). Since the actual observables are the magnetization M = gS, we
    # need to include the g factor. This is the equal-space-and-time correlation value:
    gS_squared = (2 * 1/2)^2

    expected_sum = gS_squared
    # This sum rule should hold for each sublattice, independently, and only
    # need to be taken over a single BZ (which is what sc.data contains) to hold:
    [sub_lat_sum_rules[i,i] for i in 1:Sunny.natoms(sc.crystal)] ≈ expected_sum * ones(ComplexF64,Sunny.natoms(sc.crystal))

    # The polyatomic sum rule demands going out 4 BZ's for the diamond crystal
    # since there is an atom at relative position [1/4, 1/4, 1/4]. It also
    # requires integrating over the full sampling frequency range, in this
    # case by going over both positive and negative energies.
    nbzs = (4, 4, 4)
    qs = available_wave_vectors(sc; bzsize=nbzs)
    is = intensities(sc, Sunny.QPoints(qs[:]); energies=:available_with_negative, kT=nothing)
    calculated_sum = sum(is.data) * Δq³ * sc.Δω

    # This tests that `negative_energies = true` spans exactly one sampling frequency
    expected_multi_BZ_sum = gS_squared * prod(nbzs) # ⟨S⋅S⟩
    expected_multi_BZ_sum_times_natoms = expected_multi_BZ_sum * Sunny.natoms(sc.crystal) # Nₐ×⟨S⋅S⟩
    @test calculated_sum ≈ expected_multi_BZ_sum_times_natoms
end
