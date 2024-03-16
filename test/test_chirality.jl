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
