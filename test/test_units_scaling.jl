@testitem "Dimensional Scaling" begin
    crystal = Sunny.diamond_crystal()
    sys = System(crystal, (4, 4, 4), [SpinInfo(1, S=1, g=2)], :dipole; seed=0)
    randomize_spins!(sys)

    enable_dipole_dipole!(sys, 1.0)
    E1 = energy(sys)
    ∇E1 = Sunny.energy_grad_dipoles(sys)

    enable_dipole_dipole!(sys, 2.0)
    E2 = energy(sys)
    ∇E2 = Sunny.energy_grad_dipoles(sys)

    @test E1 ≈ E2 / 2.0
    @test ∇E1 ≈ ∇E2 / 2.0
end
