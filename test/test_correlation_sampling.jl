@testitem "Asymmetric Correlations" begin

    latsize = (1,1,1)
    cryst = Sunny.fcc_primitive_crystal()
    sys = System(cryst, latsize, [SpinInfo(1; S = 1/2, g=2)], :SUN; seed = 0)
    sc = dynamical_correlations(sys; Δt = 0.1, ωmax = 10.0, nω=100, observables = [:A => I(2), :B => I(2)])
    ts = range(0,1,length = size(sc.samplebuf,6))
    #As = cos.(10π .* ts)
    #Bs = 2 .* sin.(10π .* ts)
    As = exp.(-(ts .- 0.15).^2 ./ (2 * 0.05^2))
    Bs = exp.(-(ts .- 0.35).^2 ./ (2 * 0.1^2))
    sc.samplebuf[1,1,1,1,1,:] .= As
    sc.samplebuf[2,1,1,1,1,:] .= Bs
    Sunny.accum_sample!(sc)
    real_data = real(ifft(sc.data,7))
    real_data .*= 199*199

    # Reference calculation
    q11 = zeros(199)
    q12 = zeros(199)
    q21 = zeros(199)
    q22 = zeros(199)
    for t = 0:198
        for tau = 0:198
            q11[1+t] += As[1+(t+tau)] * As[1+(tau)]
            q12[1+t] += As[1+(t+tau)] * Bs[1+(tau)]
            q21[1+t] += Bs[1+(t+tau)] * As[1+(tau)]
            q22[1+t] += Bs[1+(t+tau)] * Bs[1+(tau)]
        end
    end

    @test isapprox(q11,real_data[sc.observables.correlations[CartesianIndex(1,1)],1,1,1,1,1,:];atol = 1e-8)
    @test isapprox(q12,real_data[sc.observables.correlations[CartesianIndex(1,2)],1,1,1,1,1,:];atol = 1e-8)
    @test isapprox(q21,real_data[sc.observables.correlations[CartesianIndex(2,1)],1,1,1,1,1,:];atol = 1e-8)
    @test isapprox(q22,real_data[sc.observables.correlations[CartesianIndex(2,2)],1,1,1,1,1,:];atol = 1e-8)
end

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
        Δt = 0.05  # Time step for thermalization
        λ  = 0.1
        nsteps = 1000  # Number of steps between MC samples
        langevin = Langevin(Δt; kT, λ)
        for _ in 1:nsteps
            step!(sys, langevin)
        end
    end
    

    # Test sum rule with custom observables 
    sys = simple_model_fcc(; mode=:SUN)
    thermalize_simple_model!(sys; kT=0.1)
    S = spin_matrices(1/2)
    observables = Dict(:Sx => S[1], :Sy => S[2], :Sz => S[3])
    sc = dynamical_correlations(sys; nω=100, ωmax=10.0, Δt=0.1, apply_g=false, observables)
    add_sample!(sc, sys)
    qgrid = available_wave_vectors(sc)
    formula = intensity_formula(sc,:trace)
    vals = intensities_interpolated(sc, qgrid, formula; negative_energies=true)
    total_intensity_trace = sum(vals)
    @test isapprox(total_intensity_trace / prod(sys.latsize), 1.0; atol=1e-12)

    # Test sum rule with default observables in dipole mode 
    sys = simple_model_fcc(; mode=:dipole)
    thermalize_simple_model!(sys; kT=0.1)
    sc = dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0, apply_g=false)
    add_sample!(sc, sys)
    trace_formula = intensity_formula(sc,:trace)
    vals = intensities_interpolated(sc, qgrid, trace_formula; negative_energies=true)
    total_intensity_trace = sum(vals)
    @test isapprox(total_intensity_trace / prod(sys.latsize), 1.0; atol=1e-12)


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
    @test isapprox(total_intensity_static, total_intensity_trace; atol=1e-12)  # Order of summation can lead to very small discrepancies

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
    Δt_langevin = 0.07 
    langevin = Langevin(Δt_langevin; kT=0.1723, λ=0.1)

    sc0 = dynamical_correlations(sys; nω=25, ωmax=5.5, Δt=2Δt_langevin, calculate_errors=true)
    sc1 = dynamical_correlations(sys; nω=25, ωmax=5.5, Δt=2Δt_langevin, calculate_errors=true)
    sc2 = dynamical_correlations(sys; nω=25, ωmax=5.5, Δt=2Δt_langevin, calculate_errors=true)

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
    sys = System(Sunny.diamond_crystal(), (3,3,3), [SpinInfo(1; S=3/2, g=2)], :dipole, seed=101)
    set_exchange!(sys, 0.6498, Bond(1, 3, [0,0,0]))
    randomize_spins!(sys)

    Δt_langevin = 0.07 
    kT = 0.1723
    λ  = 0.1
    langevin = Langevin(Δt_langevin; kT, λ)

    # Thermalize
    for _ in 1:4000
        step!(sys, langevin)
    end

    sc = dynamical_correlations(sys; nω=25, ωmax=5.5, Δt=2Δt_langevin)
    add_sample!(sc, sys)
    qs = [[0.0, 0.0, 0.0], [0.2, 0.2, 0], [0.2, -0.4, 0.1]]
    data = intensities_interpolated(sc, qs, intensity_formula(sc,:trace; kT); interpolation=:linear)
    
    refdata = [4.254269966153992 -1.2843038669061228e-14 -5.394619322304545e-15 -7.216011361021513e-15 2.2240831664270433e-15 3.4163924368862625e-15 -2.213918174402409e-15 2.5176941933533583e-15 -1.6967936802985614e-15 -5.809570898878131e-16 9.221497776781295e-17 -1.1157998428324525e-15 -2.379142682899133e-15 2.8471327104209816e-15 1.7428601452501232e-15 -1.2189646774710288e-15 4.2787677309219294e-15 -1.0973563787964497e-15 7.054433862891246e-16 4.380204032382625e-17 1.7981890238059382e-15 8.714300653810433e-16 2.028726078134631e-15 3.1814113498015657e-16 -1.2725645399205768e-15; 0.691817031974635 0.022003904265315233 0.015402770605279703 0.01416546690741492 0.010268650053888165 0.018517150867665515 0.03052522221486546 0.16122718426522498 1.758858093426611 1.477056212473994 0.17592893333848564 0.9125908527099499 0.13542564389015566 0.02791055330977605 0.020926562473942117 0.020478478271579087 0.019894295068977528 0.01461089841500453 0.017161028051175586 0.02342242221257691 0.015976223751103427 0.016509502868504597 0.014721445623676686 0.020611082181934776 0.027610813222342722; 0.008860945396244786 0.027072547009685838 0.02407959453553224 0.03716552552216021 0.026340120173756747 0.03450469651089368 0.029861166686208994 0.18115847147999378 0.9207275891278851 0.9146627283758494 0.43083559486712847 8.309560489648792 2.096790402760553 2.008750574046653 3.2097755895043982 1.929761239453324 0.3246121762953159 0.09985482227045933 0.07389371860037938 0.05252516529923479 0.041887672434665145 0.04800612203272573 0.03783926774429406 0.04171088150600197 0.04288203460674578]   

    # Compare with reference 
    @test isapprox(data, refdata; atol=1e-12)
end
