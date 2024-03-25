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
        @test_broken isapprox(corr_mat, mag_corr) # Sunny bug: classical intensity off by natoms factor
        @test isapprox(corr_mat / Sunny.natoms(cryst), mag_corr)

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
