@testitem "Varying model parameters" begin
    using LinearAlgebra

    L = 6
    cryst = Sunny.chain_crystal()
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(1, 1, L))
    b = Bond(1, 1, [0, 0, 1])

    # Two additive contributions on the same bond
    set_exchange!(sys, Diagonal([1.0, 1.0, 0.0]), b, :Jxy => 1.0)
    set_exchange!(sys, Diagonal([0.0, 0.0, 1.0]), b, :Jz  => 2.0)

    # FM along z: only Jz contributes
    polarize_spins!(sys, [0, 0, 1])
    E1 = Sunny.energy(sys)
    @test E1 ≈ L * get_param(sys, :Jz)

    # Changing Jz scales energy for z-polarized state
    set_param!(sys, :Jz, 4.0)
    E2 = energy(sys)
    @test E2 ≈ L * 4.0

    # Changing Jxy doesn't affect energy
    set_param!(sys, :Jxy, 7.0)
    E3 = energy(sys)
    @test E2 ≈ E3

    # In [1, 0, 1] direction there a combination of both energies
    polarize_spins!(sys, [1/√2, 0, 1/√2])
    E3 = energy(sys)
    @test E3 ≈ (L/2) * sum(get_params(sys, [:Jxy, :Jz]))

    # Update to params should also work
    set_params!(sys, [:Jxy, :Jz], [1.0, -2.0])
    E4 = energy(sys)
    @test E4 ≈ (L/2) * (1.0 - 2.0)

    # A parameter without an explicit label will be implicitly labeled
    # :Unnamed1. This will be additive to named parameters.
    J0 = 1.3
    set_exchange!(sys, J0, b)
    E5 = energy(sys)
    @test E5 ≈ (L/2) * (1.0 - 2.0) + L * J0

    # Unnamed parameters can be overwritten in the traditional way
    J0 = 3.4
    msg = "Overwriting coupling for Bond(1, 1, [0, 0, 1])"
    @test_logs (:warn, msg) set_exchange!(sys, J0, b)
    E5 = energy(sys)
    @test E5 ≈ (L/2) * (1.0 - 2.0) + L * J0

    # Labeled parameters also
    msg = "Overwriting coupling :Jxy"
    @test_logs (:warn, msg) set_exchange!(sys, Diagonal([1.0, 1.0, 0.0]), b, :Jxy => 0.1)
    E6 = energy(sys)
    @test E6 ≈ (L/2) * (0.1 - 2.0) + L * J0
end


@testitem "Loss configuration" begin
    cryst = Sunny.square_crystal()
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]), :J1 => 0.0)

    # Something simple to exercise the plumbing
    loss1 = make_loss_fn(sys, [:J1], (; ϵ=2)) do sys, hp
        return get_param(sys, :J1) + hp.ϵ
    end
    @assert loss1([3.0]) == 5.0

    # A second loss with different hyperparameters
    loss2 = with_hyperparams(loss1, (; ϵ=1))
    @assert loss2([2.0]) == 3.0

    # The first loss should be unchanged
    @assert loss1([3.0]) == 5.0
end


@testitem "Code generation" begin
    using IOCapture
    cryst = Sunny.kagome_crystal()

    capt = IOCapture.capture() do
        Sunny.print_reference_exchanges(cryst, 1.0)
    end

    @test capt.output == """
        set_exchange!(sys, [1 0 0; 0 0 0; 0 0 0],  Bond(1, 2, [0, 0, 0]), :J1_A => 0)
        set_exchange!(sys, [0 0 0; 0 1 0; 0 0 0],  Bond(1, 2, [0, 0, 0]), :J1_B => 0)
        set_exchange!(sys, [0 0 0; 0 0 0; 0 0 1],  Bond(1, 2, [0, 0, 0]), :J1_C => 0)
        set_exchange!(sys, [0 1 0; -1 0 0; 0 0 0], Bond(1, 2, [0, 0, 0]), :J1_D => 0)
        set_exchange!(sys, [1 0 0; 0 0 0; 0 0 0],  Bond(1, 2, [0, 1, 0]), :J2_A => 0)
        set_exchange!(sys, [0 0 0; 0 1 0; 0 0 0],  Bond(1, 2, [0, 1, 0]), :J2_B => 0)
        set_exchange!(sys, [0 0 0; 0 0 0; 0 0 1],  Bond(1, 2, [0, 1, 0]), :J2_C => 0)
        set_exchange!(sys, [0 1 0; -1 0 0; 0 0 0], Bond(1, 2, [0, 1, 0]), :J2_D => 0)
        set_exchange!(sys, [1 0 0; 0 0 0; 0 0 0],  Bond(3, 3, [1, 0, 0]), :J3a_A => 0)
        set_exchange!(sys, [0 0 0; 0 1 0; 0 0 0],  Bond(3, 3, [1, 0, 0]), :J3a_B => 0)
        set_exchange!(sys, [0 0 0; 0 0 0; 0 0 1],  Bond(3, 3, [1, 0, 0]), :J3a_C => 0)
        set_exchange!(sys, [1 0 0; 0 0 0; 0 0 0],  Bond(1, 1, [1, 0, 0]), :J3b_A => 0)
        set_exchange!(sys, [0 0 0; 0 1 0; 0 0 0],  Bond(1, 1, [1, 0, 0]), :J3b_B => 0)
        set_exchange!(sys, [0 0 0; 0 0 0; 0 0 1],  Bond(1, 1, [1, 0, 0]), :J3b_C => 0)
        set_exchange!(sys, [0 1 0; 1 0 0; 0 0 0],  Bond(1, 1, [1, 0, 0]), :J3b_D => 0)
        """
end
