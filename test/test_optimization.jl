@testitem "Optimization" begin

    # H = -∑Sᵢ⋅Sⱼ - ∑(Sᵢᶻ)² on 2D square lattice (z-polarized ground state)
    function simple_sys(; mode=:dipole, dims=(4,4,1), seed=nothing)
        cryst = Crystal(lattice_vectors(1,1,2,90,90,90), [[0,0,0]])
        sys = System(cryst, dims, [SpinInfo(1,S=1)], mode; seed) 
        set_exchange!(sys, -1, Bond(1,1,[1,0,0]))
        set_onsite_coupling!(sys, -1*(Sunny.spin_operators(sys, 1)[3])^2, 1)
        sys
    end

    function is_z_polarized(sys; S=1.0)
        for site in all_sites(sys)
            !isapprox(abs(sys.dipoles[site][3]), S) && return false
        end
        true
    end

    sys_dip = simple_sys(; mode=:dipole, seed=101)
    sys_sun = simple_sys(; mode=:SUN, seed=102)

    # Thermalize near ground state a little bit and apply `minimize_energy_touchup!`
    Δt = 0.05
    integrator = Langevin(Δt; λ=0.1, kT = 0.1)

    for _ in 1:1000
        step!(sys_dip, integrator)
        step!(sys_sun, integrator)
    end

    Sunny.minimize_energy_touchup!(sys_dip)
    Sunny.minimize_energy_touchup!(sys_sun)

    @test is_z_polarized(sys_dip)
    @test is_z_polarized(sys_sun)

    # Randomize initial condition and use `minimize_energy!` 
    randomize_spins!(sys_dip)
    randomize_spins!(sys_sun)

    Sunny.minimize_energy!(sys_dip)
    Sunny.minimize_energy!(sys_sun)

    @test is_z_polarized(sys_dip)
    @test is_z_polarized(sys_sun)
end