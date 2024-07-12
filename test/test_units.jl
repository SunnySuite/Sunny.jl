@testitem "Units" begin
    function f(units)
        crystal = Sunny.diamond_crystal(; a=units.angstrom)
        sys = System(crystal, (4, 4, 4), [SpinInfo(1, S=1, g=2)], :dipole; seed=0)
        randomize_spins!(sys)
        set_exchange!(sys, 2 * units.K, Bond(1, 2, [0,0,0]))
        set_field!(sys, [0, 0, 1] * units.T)
        enable_dipole_dipole!(sys, units.vacuum_permeability)
        E = energy(sys) / units.meV
        ∇E = Sunny.energy_grad_dipoles(sys) / units.THz
        return (E, ∇E)
    end

    (E1, ∇E1) = f(Units(:meV, :angstrom))
    (E2, ∇E2) = f(Units(:THz, :nm))

    @test E1 ≈ E2
    @test ∇E1 ≈ ∇E2
end
