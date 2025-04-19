@testitem "Interaction lookup" begin
    using LinearAlgebra
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
    cryst = Crystal(latvecs, positions, 1)

    function run_test(mode)
        sys = System(cryst, [1 => Moment(s=1, g=2), 2 => Moment(s=1, g=2)], mode)
        randomize_spins!(sys)

        D1 = hermitianpart(randn(3, 3))
        D2 = hermitianpart(randn(3, 3))
        set_onsite_coupling!(sys, S -> S'*D1*S, 1)
        set_onsite_coupling!(sys, S -> S'*D2*S, 2)

        sys2 = reshape_supercell(sys, [1 1 0; 0 2 0; 0 0 1])
        @test energy_per_site(sys) ≈ energy_per_site(sys2)

        @test Sunny.get_quadratic_anisotropy(sys, 1) ≈ D1
        @test Sunny.get_quadratic_anisotropy(sys, 2) ≈ D2
        @test Sunny.get_quadratic_anisotropy(sys2, 1) ≈ D1
        @test Sunny.get_quadratic_anisotropy(sys2, 2) ≈ D2

        J1 = hermitianpart(randn(3, 3))
        J2 = hermitianpart(randn(3, 3))
        set_exchange!(sys2, J1, Bond(1, 2, [0, 0, 0]))
        @test Sunny.get_exchange(sys2, Bond(1, 2, [0, 0, 0])) ≈ J1

        sys3 = to_inhomogeneous(sys2)
        set_exchange_at!(sys3, J2, (1, 1, 1, 1), (1, 1, 1, 2); offset=[0, 0, 0])

        @test Sunny.get_exchange_at(sys3, (1, 1, 1, 1), (1, 1, 1, 3); offset=[0, 0, 0]) ≈ J1
        @test Sunny.get_exchange_at(sys3, (1, 1, 1, 1), (1, 1, 1, 2); offset=[0, 0, 0]) ≈ J2
    end

    foreach(run_test, (:dipole, :dipole_uncorrected, :SUN))
end
