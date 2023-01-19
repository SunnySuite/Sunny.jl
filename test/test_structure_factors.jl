@testitem "Structure factors" begin

    function simple_model_sf(; mode, seed=111)
        latsize = (4,4,4)
        J = 1.0
        cryst = Sunny.fcc_primitive_crystal()
        S = mode==:SUN ? 1/2 : 1
        κ = mode==:SUN ? 2 : 1
        sys = System(cryst, latsize, [SiteInfo(1; S)]; mode, seed)
        sys.κs .= κ
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        return sys
    end

    function thermalize_simple_model!(sys; kT)
        Δt = 0.05  # Time step for thermalization
        λ  = 0.1
        nsteps = 1000  # Number of steps between MC samples
        integrator = LangevinHeunP(kT, λ, Δt)
        for _ ∈ 1:nsteps
            step!(sys, integrator)
        end
    end
    

    # Test sum rule with custom observables 
    sys = simple_model_sf(; mode=:SUN)
    S = Sunny.spin_matrices(2)
    ops = cat(S...; dims=3)
    sf = StructureFactor(sys; ωmax=10.0, gfactor=false, ops)
    thermalize_simple_model!(sys; kT=0.1)
    add_trajectory!(sf, sys)
    intensities = intensity_grid(sf; contraction=:trace, negative_energies=true)
    total_intensity_trace = sum(intensities)
    @test isapprox(total_intensity_trace / prod(sys.latsize), 1.0; atol=1e-12)


    # Test sum rule with default observables in dipole mode 
    sys = simple_model_sf(; mode=:dipole)
    sf = StructureFactor(sys; ωmax=10.0, gfactor=false)
    thermalize_simple_model!(sys; kT=0.1)
    add_trajectory!(sf, sys)
    intensities = intensity_grid(sf; contraction=:trace, negative_energies=true)
    total_intensity_trace = sum(intensities)
    @test isapprox(total_intensity_trace / prod(sys.latsize), 1.0; atol=1e-12)


    # Test depolarize reduces intensity
    intensities = intensity_grid(sf; contraction=:depolarize, negative_energies=true)
    total_intensity_depolarized = sum(intensities)
    @test total_intensity_depolarized < total_intensity_trace


    # Test diagonal elements are approximately real (at one wave vector)
    for α ∈ 1:3
        intensities = get_intensities(sf, (0.25, 0.5, 0); contraction=(α,α))
        @test sum(imag(intensities)) < 1e-15
    end


    # Test form factor correction works and is doing something. ddtodo: add example with sublattice
    formfactors = [FormFactor(1, "Fe2")]
    intensities = intensity_grid(sf; contraction=:trace, formfactors, negative_energies=true)
    total_intensity_ff = sum(intensities)
    @test total_intensity_ff != total_intensity_trace


    # Test path function and interpolation working
    points = [(0, 0, 0), (0, 1, 0), (1, 1, 0)]
    intensities = path(sf, points; density=20, interpolation=:linear)
    @test length(size(intensities)) == 2 


    # Test static intensities working
    qs = Sunny.qgrid(sf) 
    static_intensities = get_static_intensities(sf, qs; negative_energies=true)
    total_intensity_static = sum(static_intensities)
    @test total_intensity_static ≈ total_intensity_trace  # Order of summation can lead to very small discrepancies
end


@testitem "Structure factor reference" begin
    using DelimitedFiles

    function diamond_model(; J, dims = (3,3,3), kwargs...)
        crystal = Sunny.diamond_crystal()
        S = 3/2
        sys = System(crystal, dims, [SiteInfo(1; S)]; mode=:dipole, kwargs...)
        sys.κs .= 3/2  # This can be removed when potential bug is fixed
        set_exchange!(sys, J, Bond(1, 3, [0,0,0]))
        randomize_spins!(sys)
        return sys
    end

    seed = 101
    J = Sunny.meV_per_K * 7.5413 
    sys = diamond_model(; J, seed)

    Δt_therm = 0.07 
    kT = Sunny.meV_per_K * 2. # Units of meV
    λ  = 0.1
    integrator = LangevinHeunP(kT, λ, Δt_therm)

    # Thermalize
    for _ ∈ 1:3000
        step!(sys, integrator)
    end

    # Calculate a path
    sampler = LangevinSampler(integrator, 1000)
    sf = calculate_structure_factor(sys, sampler; numω=25, ωmax=5.5, numsamps=1)
    points = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.0, 0.0]]
    intensities = path(sf, points; interpolation = :linear, contraction = :trace, kT)

    # Compare with reference 
    refdata = readdlm(joinpath(@__DIR__, "..", "src", "StructureFactors", "data", "sf_ref.dat"))
    @test isapprox(intensities, refdata; atol=1e-12)
end