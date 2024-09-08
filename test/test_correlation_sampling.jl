@testitem "Correlation sampling" begin
    using LinearAlgebra

    function simple_model_fcc(; mode, seed=111)
        dims = (4, 4, 4)
        J = 1.0

        # FCC with nonstandard, primitive lattice vectors
        latvecs = [[1, 1, 0] [0, 1, 1] [1, 0, 1]] / 2
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)

        s = mode==:SUN ? 1/2 : 1
        κ = mode==:SUN ? 2 : 1
        sys = System(cryst, [1 => Moment(; s, g=2)], mode; dims, seed)
        sys.κs .= κ
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        return sys
    end

    function thermalize_simple_model!(sys; kT)
        dt = 0.05
        nsteps = 1000
        langevin = Langevin(dt; damping=0.1, kT)
        for _ in 1:nsteps
            step!(sys, langevin)
        end
    end

    # Test classical sum rule in SU(N) mode
    sys = simple_model_fcc(; mode=:SUN)
    thermalize_simple_model!(sys; kT=0.1)
    energies=range(0, 10.0, 100)
    sc = SampledCorrelations(sys; dt=0.08, energies, measure=ssf_trace(sys; apply_g=false))
    Δω = sc.Δω
    add_sample!(sc, sys)
    qgrid = Sunny.available_wave_vectors(sc)[:]
    Δq³ = 1/prod(sys.dims) # Fraction of a BZ
    Sqw = intensities(sc, qgrid; energies=:available_with_negative, kT=nothing)
    expected_sum_rule = Sunny.norm2(sys.dipoles[1]) # S^2 classical sum rule
    @test isapprox(sum(Sqw.data) * Δq³ * Δω, expected_sum_rule; atol=1e-12)

    # Test classical sum rule in dipole mode
    sys = simple_model_fcc(; mode=:dipole)
    thermalize_simple_model!(sys; kT=0.1)
    sc = SampledCorrelations(sys; dt=0.08, energies, measure=ssf_trace(sys; apply_g=false))
    add_sample!(sc, sys)
    Sqw = intensities(sc, qgrid; energies=:available_with_negative, kT=nothing)
    total_intensity_trace = sum(Sqw.data)
    expected_sum_rule = Sunny.norm2(sys.dipoles[1]) # S^2 classical sum rule
    @test isapprox(total_intensity_trace * Δq³ * Δω, expected_sum_rule; atol=1e-12)

    # Test perp reduces intensity
    sc.measure = ssf_perp(sys; apply_g=false)
    Sqw_perp = intensities(sc, qgrid; energies=:available_with_negative, kT=nothing)
    total_intensity_unpolarized = sum(Sqw_perp.data)
    @test total_intensity_unpolarized < total_intensity_trace

    # Test diagonal elements are approximately real (at one wave vector)
    function ssf_diag_imag(sys::System{N}; apply_g=false) where N
        return ssf_custom(sys; apply_g) do _, ssf
            tr(imag(ssf))
        end
    end
    sc.measure = ssf_diag_imag(sys; apply_g=false)
    res = intensities(sc, [[0.25, 0.5, 0]]; energies=:available, kT=nothing)
    res.data
    @test sum(res.data) < 1e-14

    # Test form factor correction works and is doing something.
    formfactors = [1 => FormFactor("Fe2")]
    sc.measure = ssf_trace(sys; apply_g=false, formfactors)
    res = intensities(sc, qgrid; energies=:available_with_negative, kT=nothing)
    total_intensity_ff = sum(res.data)
    @test total_intensity_ff != total_intensity_trace

    # Test static from dynamic intensities working
    sc.measure = ssf_trace(sys; apply_g=false)
    res_static = intensities_static(sc, qgrid; kT=nothing)
    total_intensity_static = sum(res_static.data)
    @test isapprox(total_intensity_static, total_intensity_trace * sc.Δω; atol=1e-9)  # Order of summation can lead to very small discrepancies

    # Test quantum-to-classical increases intensity
    res_static_c2q = intensities_static(sc, qgrid; kT=0.1)
    total_intensity_static_c2q = sum(res_static_c2q.data)
    @test total_intensity_static_c2q > total_intensity_static

    # Test static intensities working
    sys = simple_model_fcc(; mode=:dipole)
    thermalize_simple_model!(sys; kT=0.1)
    res = SampledCorrelationsStatic(sys; measure=ssf_trace(sys; apply_g=false))
    add_sample!(res, sys)
    true_static_vals = intensities_static(res, qgrid)
    true_static_total = sum(true_static_vals.data)
    @test isapprox(true_static_total / prod(sys.dims), 1.0; atol=1e-12)

    # Test whether two distinct wave vectors referencing the same underlying
    # data point in sc.data are in fact treated differently.
    formfactors = [1 => FormFactor("Fe2")]
    sc.measure = ssf_trace(sys; apply_g=false, formfactors)
    res = intensities_static(sc, [[0, 0, 1/2], [1, 0, 1/2]]; kT=nothing)
    @test res.data[1] != res.data[2]
end

@testitem "Merge correlations" begin
    # Set up a system.
    sys = System(Sunny.diamond_crystal(), [1 => Moment(s=3/2, g=2)], :dipole; dims=(2, 2, 2), seed=101)
    set_exchange!(sys, 0.6498, Bond(1, 3, [0, 0, 0]))
    randomize_spins!(sys)

    # Set up Langevin sampler.
    dt_langevin = 0.07
    langevin = Langevin(dt_langevin; damping=0.1, kT=0.1723)

    measure = ssf_trace(sys)
    sc0 = SampledCorrelations(sys; measure, energies=range(0, 5.5, 25), dt=0.12, calculate_errors=true)
    sc1 = SampledCorrelations(sys; measure, energies=range(0, 5.5, 25), dt=0.12, calculate_errors=true)
    sc2 = SampledCorrelations(sys; measure, energies=range(0, 5.5, 25), dt=0.12, calculate_errors=true)

    for _ in 1:4_000
        step!(sys, langevin)
    end
    add_sample!(sc0, Sunny.clone_system(sys))
    add_sample!(sc1, Sunny.clone_system(sys))

    for _ in 1:2
        for _ in 1:4_000
            step!(sys, langevin)
        end
        add_sample!(sc0, Sunny.clone_system(sys))
        add_sample!(sc2, Sunny.clone_system(sys))
    end

    # Merge correlations and check if result equal to running calculation.
    sc_merged = merge_correlations([sc1, sc2])
    @test sc0.data ≈ sc_merged.data
    @test sc0.M ≈ sc_merged.M
end

@testitem "Sampled correlations reference" begin
    sys = System(Sunny.diamond_crystal(), [1 => Moment(s=3/2, g=2)], :dipole; dims=(2, 3, 4), seed=101)
    set_exchange!(sys, 0.6498, Bond(1, 3, [0,0,0]))
    randomize_spins!(sys)

    sc = SampledCorrelations(sys; energies=range(0.0, 5.5, 10), dt=0.12, measure=ssf_perp(sys))
    add_sample!(sc, sys)
    qs = [[0.0, 0.0, 0.0], [-0.2, 0.4, -0.1]]
    res = intensities(sc, qs; energies=:available, kT=nothing)

    # Compare with reference
    data_golden = [33.52245944537883 31.523781055757002; 16.76122972268928 16.188214427928443; 1.6337100022968968e-14 5.3112747876921045; 1.590108516238891e-13 1.8852219123773621; -1.5194916032483394e-14 0.06400012935688847; -6.217248937900877e-15 -0.01803103943766904; 7.863406895322359e-14 -0.04445301974061088; 7.014086281484114e-14 0.0025512102338097653; 3.195591939212742e-14 -0.02515685630480813; 3.6201681383269686e-14 0.023924996100518413]
    @test res.data ≈ data_golden
end
