@testitem "Correlation sampling" begin
    using LinearAlgebra

    function simple_model_fcc(; mode, seed=111)
        latsize = (4, 4, 4)
        J = 1.0

        # FCC with nonstandard, primitive lattice vectors
        latvecs = [[1, 1, 0] [0, 1, 1] [1, 0, 1]] / 2
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)

        S = mode==:SUN ? 1/2 : 1
        κ = mode==:SUN ? 2 : 1
        sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], mode; seed)
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
    

    # Test sum rule with custom observables 
    sys = simple_model_fcc(; mode=:SUN)
    thermalize_simple_model!(sys; kT=0.1)
    S = spin_matrices(1/2)
    observables = Dict(:Sx => S[1], :Sy => S[2], :Sz => S[3])
    sc = dynamical_correlations(sys; nω=100, ωmax=10.0, dt=0.1, apply_g=false, observables)
    Δω = sc.Δω
    add_sample!(sc, sys)
    qgrid = available_wave_vectors(sc)
    Δq³ = 1/prod(sys.latsize) # Fraction of a BZ
    formula = intensity_formula(sc,:trace)
    Sqw = intensities_interpolated(sc, qgrid, formula; negative_energies=true)
    expected_sum_rule = Sunny.norm2(sys.dipoles[1]) # S^2 classical sum rule
    @test isapprox(sum(Sqw) * Δq³ * Δω, expected_sum_rule; atol=1e-12)
    params = unit_resolution_binning_parameters(sc;negative_energies=true)
    Sqw_integrated, counts = intensities_binned(sc, params, formula)
    @test isapprox(sum(Sqw_integrated), expected_sum_rule; atol=1e-12)

    # Test sum rule with default observables in dipole mode 
    sys = simple_model_fcc(; mode=:dipole)
    thermalize_simple_model!(sys; kT=0.1)
    sc = dynamical_correlations(sys; dt=0.1, nω=100, ωmax=10.0, apply_g=false)
    add_sample!(sc, sys)
    trace_formula = intensity_formula(sc,:trace)
    Sqw = intensities_interpolated(sc, qgrid, trace_formula; negative_energies=true)
    total_intensity_trace = sum(Sqw)
    expected_sum_rule = Sunny.norm2(sys.dipoles[1]) # S^2 classical sum rule
    @test isapprox(sum(Sqw) * Δq³ * Δω, expected_sum_rule; atol=1e-12)
    # Test binned version doesn't require ΔqΔω measure
    params = unit_resolution_binning_parameters(sc;negative_energies=true)
    Sqw_integrated, counts = intensities_binned(sc, params, trace_formula)
    @test isapprox(sum(Sqw_integrated), expected_sum_rule; atol=1e-12)

    # Test perp reduces intensity
    perp_formula = intensity_formula(sc,:perp)
    vals = intensities_interpolated(sc, qgrid, perp_formula; negative_energies=true) 
    total_intensity_unpolarized = sum(vals)
    @test total_intensity_unpolarized < total_intensity_trace


    # Test diagonal elements are approximately real (at one wave vector)
    diag_elems = [(α,α) for α in keys(sc.observables.observable_ixs)]
    formula_imaginary_parts = intensity_formula(sc,diag_elems) do k,ω,corr
        sum(abs.(imag.(corr)))
    end
    intensities_symmetric = intensities_interpolated(sc, [(0.25, 0.5, 0)], formula_imaginary_parts)
    @test sum(imag(intensities_symmetric)) < 1e-15


    # Test form factor correction works and is doing something. ddtodo: add example with sublattice
    formfactors = [FormFactor("Fe2")]
    vals = intensities_interpolated(sc, qgrid, intensity_formula(sc,:trace; formfactors); negative_energies=true)
    total_intensity_ff = sum(vals)
    @test total_intensity_ff != total_intensity_trace


    # Test path function and interpolation working (no correctness implied here)
    qs, _ = reciprocal_space_path(sc.crystal, [(0, 0, 0), (0, 1, 0), (1, 1, 0)], 20)
    vals = intensities_interpolated(sc, qs, perp_formula; interpolation=:linear)
    @test length(size(vals)) == 2 


    # Test static from dynamic intensities working
    static_vals = instant_intensities_interpolated(sc, qgrid, trace_formula; negative_energies=true)
    total_intensity_static = sum(static_vals)
    @test isapprox(total_intensity_static, total_intensity_trace; atol=1e-9)  # Order of summation can lead to very small discrepancies

    # Test instant intensities working
    sys = simple_model_fcc(; mode=:dipole)
    thermalize_simple_model!(sys; kT=0.1)
    ic = instant_correlations(sys; apply_g=false)
    add_sample!(ic, sys)
    true_static_vals = instant_intensities_interpolated(ic, qgrid, intensity_formula(ic,:trace))
    true_static_total = sum(true_static_vals)
    @test isapprox(true_static_total / prod(sys.latsize), 1.0; atol=1e-12)
end

@testitem "Merge correlations" begin
    # Set up a system.
    sys = System(Sunny.diamond_crystal(), (2,2,2), [SpinInfo(1; S=3/2, g=2)], :dipole, seed=101)
    set_exchange!(sys, 0.6498, Bond(1, 3, [0,0,0]))
    randomize_spins!(sys)

    # Set up Langevin sampler.
    dt_langevin = 0.07 
    langevin = Langevin(dt_langevin; damping=0.1, kT=0.1723)

    sc0 = dynamical_correlations(sys; nω=25, ωmax=5.5, dt=2dt_langevin, calculate_errors=true)
    sc1 = dynamical_correlations(sys; nω=25, ωmax=5.5, dt=2dt_langevin, calculate_errors=true)
    sc2 = dynamical_correlations(sys; nω=25, ωmax=5.5, dt=2dt_langevin, calculate_errors=true)

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
    sys = System(Sunny.diamond_crystal(), (2,3,4), [SpinInfo(1; S=3/2, g=2)], :dipole, seed=101)
    set_exchange!(sys, 0.6498, Bond(1, 3, [0,0,0]))
    randomize_spins!(sys)
    
    sc = dynamical_correlations(sys; nω=10, ωmax=5.5, dt=0.14)
    add_sample!(sc, sys)
    qs = [[0.0, 0.0, 0.0], [-0.2, 0.4, -0.1]]
    data = intensities_interpolated(sc, qs, intensity_formula(sc, :trace; kT=0.1723); interpolation=:linear)
    # println(round.(data; digits=10))

    # Compare with reference
    data_golden = [52.0366191891 92.1657422033 0.0 -0.0 -0.0 0.0 -0.0 -0.0 -0.0 -0.0; 32.9297460003 86.5770399929 73.5353464564 33.8358366952 20.3323754023 13.7426991056 5.4654945884 1.8418530116 1.1247602598 0.5785036261]
    @test data ≈ data_golden
end
