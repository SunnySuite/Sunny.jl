@testitem "Contractors" begin
    using LinearAlgebra
    cryst = Crystal(I(3), [[0.,0,0]], 1)
    sys = System(cryst, (1,1,1), [SpinInfo(1,S=1,g=2)], :SUN)
    sc = instant_correlations(sys)
    @test_nowarn intensity_formula(sc,(:Sx,:Sz))
    @test_nowarn intensity_formula(sc,:trace)
    @test_nowarn intensity_formula(sc,:perp)
    @test_nowarn intensity_formula(sc,:full)
    swt = SpinWaveTheory(sys)
    @test_nowarn intensity_formula(swt,(:Sx,:Sz); kernel = delta_function_kernel)
    @test_nowarn intensity_formula(swt,:trace; kernel = delta_function_kernel)
    @test_nowarn intensity_formula(swt,:perp; kernel = delta_function_kernel)
    @test_nowarn intensity_formula(swt,:full; kernel = delta_function_kernel)

end

@testitem "Dipole Factor Ordering" begin
    using LinearAlgebra
    obs = Sunny.parse_observables(3; observables=nothing, correlations=nothing, g=nothing)
    dipoleinfo = Sunny.DipoleFactor(obs)
    fake_intensities = [1. 3 5; 0 7 11; 0 0 13]
    # This is the order expected by contract(...)
    fake_dipole_elements = fake_intensities[[1,4,5,7,8,9]]
    k = Sunny.Vec3(1,2,3)
    mat = Sunny.polarization_matrix(k)
    out = Sunny.contract(fake_dipole_elements, k, dipoleinfo)
    out_true = dot(Hermitian(fake_intensities), mat)
    @test isapprox(out, out_true)
end


@testitem "Magnetization Observables" begin
    using LinearAlgebra

    g_factor = [0.8 3.2 5.1; -0.3 0.6 -0.1; 1.0 0. 0.]
    g_factor_2 = [1.8 2.2 -3.1; -0.3 2.6 0.1; 1.0 0. 0.3]
    cryst = Crystal(I(3), [[0.,0,0],[0.1,0.2,0.3]], 1)

    for mode = [:SUN, :dipole]
        sys = System(cryst, (1,1,1), [SpinInfo(1,S=3/2, g=g_factor), SpinInfo(2, S=mode == :SUN ? 3/2 : 1/2, g=g_factor_2)], mode)

        set_dipole!(sys,[1,3,2], (1,1,1,1))
        set_dipole!(sys,[3,4,5], (1,1,1,2))

        # Dipole magnetization observables (classical)
        sc = instant_correlations(sys)
        add_sample!(sc, sys)

        formula = intensity_formula(sc,:full)
        corr_mat = instant_intensities_interpolated(sc, [[0,0,0]], formula,interpolation=:round)[1]

        # Compute magnetization correlation "by hand", averaging over sites
        mag_corr = sum([sys.gs[i] * sys.dipoles[i] * (sys.gs[j] * sys.dipoles[j])' for i = 1:2, j = 1:2]) / Sunny.natoms(cryst)
        mag_corr_time_natoms = mag_corr * Sunny.natoms(cryst)
        @test isapprox(corr_mat, mag_corr_time_natoms)

        # Spin wave theory only gives the "transverse part" which is difficult to calculate.
        # So we compare spin correlations vs magnetization correlations externally.
        sys_homog_g = System(cryst, (1,1,1), [SpinInfo(1,S=3/2,g=g_factor), SpinInfo(2, S=mode == :SUN ? 3/2 : 1/2, g=g_factor)], mode)
        swt = SpinWaveTheory(sys_homog_g; apply_g = false)
        formula = intensity_formula(swt,:full, kernel=delta_function_kernel)
        disp, is_spin_spin = intensities_bands(swt,[[0,0,0]], formula)

        swt = SpinWaveTheory(sys_homog_g; apply_g = true)
        formula = intensity_formula(swt, :full, kernel=delta_function_kernel)
        disp, is_mag_mag = intensities_bands(swt, [[0,0,0]], formula)

        # TODO: Can the ground truth be calculated explicitly from the sys.dipoles, and inhomogeneous g_factor be used?
        @test isapprox(g_factor * is_spin_spin[1] * g_factor', is_mag_mag[1])
        @test isapprox(g_factor * is_spin_spin[2] * g_factor', is_mag_mag[2])
    end

end
@testitem "Available Energies Dirac Identity" begin
     # Create a dummy SampledCorrelations object
    latsize = (1,1,1)
    cryst = Sunny.cubic_crystal()
    sys = System(cryst, latsize, [SpinInfo(1; S = 1/2, g=2)], :SUN; seed = 0)
    dt = 0.1
    sc = dynamical_correlations(sys; dt, ωmax=10, nω=100)

    ωs = available_energies(sc;negative_energies=true)
    dts = 0:(dt * sc.measperiod):3
    vals = sum(exp.(im .* ωs .* dts'),dims = 1)[:]

    # Verify it made a delta function
    @test vals[1] ≈ length(ωs)
    @test all(isapprox.(0,vals[2:end];atol = 1e-12))
end

@testitem "Polyatomic sum rule" begin
    sys = System(Sunny.diamond_crystal(),(4,1,1),[SpinInfo(1,S=1/2,g=2)],:SUN,seed=1)
    randomize_spins!(sys)
    sc = dynamical_correlations(sys; dt=1, nω=3, ωmax=1)
    add_sample!(sc, sys)

    sum_rule_ixs = Sunny.Trace(sc.observables).indices
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
    @test [sub_lat_sum_rules[i,i] for i in 1:Sunny.natoms(sc.crystal)] ≈ expected_sum * ones(ComplexF64,Sunny.natoms(sc.crystal))

    formula = intensity_formula(sc,:trace)
    # The polyatomic sum rule demands going out 4 BZ's for the diamond crystal
    # since there is an atom at relative position [1/4, 1/4, 1/4]. It also
    # requires integrating over the full sampling frequency range, in this
    # case by going over both positive and negative energies.
    params_pasr = unit_resolution_binning_parameters(sc;negative_energies = true)
    params_pasr.binstart[1:3] .= -params_pasr.binwidth[1:3] ./ 2
    params_pasr.binend[1:3] .+= 3

    # This should result in spanning exactly 4x4x4 BZ's
    nbzs = (params_pasr.binwidth .* params_pasr.numbins)[1:3]
    @test nbzs ≈ [4.0,4.0,4.0]

    # This tests that `negative_energies = true` spans exactly one sampling frequency
    nfs = params_pasr.binwidth[4] * params_pasr.numbins[4] / (sc.Δω * size(sc.data,7))
    @test nfs ≈ 1
    is, counts = intensities_binned(sc,params_pasr,formula)
    expected_multi_BZ_sum = gS_squared * prod(nbzs) * nfs # ⟨S⋅S⟩
    expected_multi_BZ_sum_times_natoms = expected_multi_BZ_sum * Sunny.natoms(sc.crystal) # Nₐ×⟨S⋅S⟩
    @test sum(is) ≈ expected_multi_BZ_sum_times_natoms
end
