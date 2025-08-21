@testitem "Accumulated model parameters" begin
    using LinearAlgebra

    L = 6
    cryst = Sunny.chain_crystal()
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(L, 1, 1))
    b = Bond(1, 1, [1, 0, 0])

    # Two additive contributions on the same bond
    set_exchange!(sys, Diagonal([1.0, 1.0, 0.0]), b, :Jxy => 1.0)
    set_exchange!(sys, Diagonal([0.0, 0.0, 1.0]),  b, :Jz  => 2.0)

    # FM along z: only Jz contributes
    polarize_spins!(sys, [0, 0, 1])
    E1 = Sunny.energy(sys)
    @test E1 ≈ 2L * get_param(sys, :Jz)

    # Changing Jz scales energy for z-polarized state
    set_param!(sys, :Jz, 4.0)
    E2 = energy(sys)
    @test E2 ≈ 2L * 4.0

    # Changing Jxy doesn't affect energy
    set_param!(sys, :Jxy, 7.0)
    E3 = energy(sys)
    @test E2 ≈ E3

    # In [1, 0, 1] direction there a combination of both energies
    polarize_spins!(sys, [1/√2, 0, 1/√2])
    E3 = energy(sys)
    @test E3 ≈ (2L/2) * (get_param(sys, :Jxy) + get_param(sys, :Jz))

    # Bulk updates to params should also work
    set_params!(sys, [:Jxy, :Jz], [1.0, -2.0])
    E4 = energy(sys)
    @test E4 ≈ (2L/2) * (1.0 - 2.0)
end
