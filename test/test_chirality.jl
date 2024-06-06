@testitem "Spin precession handedness" begin
    using LinearAlgebra

    for units in (Units.meV, Units.theory)
        crystal = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0, 0, 0]])
        sys_dip = System(crystal, (1, 1, 1), [SpinInfo(1; S=1, g=2)], :dipole; units)
        sys_sun = System(crystal, (1, 1, 1), [SpinInfo(1; S=1, g=2)], :SUN; units)
    
        B = [0, 0, 1] / abs(units.μB)
        set_external_field!(sys_dip, B)
        set_external_field!(sys_sun, B)
    
        ic = [1/√2, 0, 1/√2]
        set_dipole!(sys_dip, ic, (1, 1, 1, 1))
        set_dipole!(sys_sun, ic, (1, 1, 1, 1))
    
        integrator = ImplicitMidpoint(0.05)
        for _ in 1:5
            step!(sys_dip, integrator)
            step!(sys_sun, integrator)
        end
    
        dip_is_lefthanded = B ⋅ (ic × magnetic_moment(sys_dip, (1,1,1,1))) < 0
        sun_is_lefthanded = B ⋅ (ic × magnetic_moment(sys_sun, (1,1,1,1))) < 0
    
        @test dip_is_lefthanded == sun_is_lefthanded == true
    end
end


@testitem "DM chain" begin
    latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]], "P1")
    sys = System(cryst, (1,1,1), [SpinInfo(1,S=1,g=1)], :dipole; units=Units.theory)
    D = 1
    B = 10.0
    set_exchange!(sys, dmvec([0, 0, D]), Bond(1, 1, [0, 0, 1]))
    set_external_field!(sys, [0, 0, B])

    # Above the saturation field, the ground state is fully polarized, with no
    # energy contribution from the DM term.

    randomize_spins!(sys)
    minimize_energy!(sys)
    @test energy_per_site(sys) ≈ -B
    qpoints = [[0, 0, -1/2], [0, 0, 1/2]]
    path, xticks = reciprocal_space_path(cryst, qpoints, 50)
    swt = SpinWaveTheory(sys)
    formula = intensity_formula(swt, :trace; kernel=delta_function_kernel)
    disp, intens = intensities_bands(swt, path, formula)
    disp_ref = [B + 2D*sin(2π*q[3]) for q in path]
    intens_ref = [1.0 for _ in path]
    @test disp[:,1] ≈ disp_ref
    @test intens[:,1] ≈ intens_ref

    # Below the saturation field, the ground state is a canted spiral

    sys2 = resize_supercell(sys, (1, 1, 4))
    B = 1
    set_external_field!(sys2, [0, 0, B])
    randomize_spins!(sys2)
    minimize_energy!(sys2)
    @test energy_per_site(sys2) ≈ -5/4
    swt = SpinWaveTheory(sys2)
    formula = intensity_formula(swt, :trace; kernel=delta_function_kernel)
    qs = [[0,0,-1/3], [0,0,1/3]]
    disp2, intens2 = intensities_bands(swt, qs, formula)
    disp2_ref = [3.0133249314 2.5980762316 1.3228756763 0.6479760935
                 3.0133249314 2.5980762316 1.3228756763 0.6479760935]
    intens2_ref = [0.0292617379 0.4330127014 0.0 0.8804147011
                   0.5292617379 0.4330127014 0.0 0.3804147011]
    @test disp2 ≈ disp2_ref
    @test intens2 ≈ intens2_ref

    # Perform the same calculation with Single-Q functions

    sys3 = resize_supercell(sys2, (1, 1, 1))
    axis = [0, 0, 1]
    randomize_spins!(sys3)
    k = Sunny.minimize_energy_spiral!(sys3, axis; k_guess=randn(3))
    @test k[3] ≈ 3/4
    @test Sunny.spiral_energy_per_site(sys3; k, axis) ≈ -5/4
    swt = SpinWaveTheory(sys3)
    formula = Sunny.intensity_formula_spiral(swt, :perp; k, axis, kernel=delta_function_kernel)
    disp3, intens3 = intensities_bands(swt, qs, formula)
    # TODO: Make the dispersions and intensities match! Some variation of this
    # should work.
    #=
    @assert disp2 ≈ disp3
    @assert intens2 ≈ intens3
    =#
end
