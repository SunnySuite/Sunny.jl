@testitem "Squared error with rescaling" begin
    u = ComplexF64[1 + 2im, -3 + im, 2 - 4im]
    α = 2 - 3im
    v = α .* u
    (; error, scale) = squared_error_fitted(u, v; scale=true)
    @test error ≈ 0 atol=1e-12
    @test scale ≈ α atol=1e-12

    (; error, scale, shift) = squared_error_fitted([1, 2, 3], [1, 1, 1]; scale=true, shift=true)
    @test error ≈ 1 atol=1e-12
    @test scale ≈ 0 atol=1e-12
    @test shift ≈ 1 atol=1e-12

    (; error, scale, shift) = squared_error_fitted([2, 2, 2], [1, 1, 1]; scale=true, shift=true)
    @test error ≈ 0 atol=1e-12
    @test scale ≈ 1 atol=1e-12
    @test shift ≈ -1.0 atol=1e-12

    (; error, scale, shift) = squared_error_fitted([2, 2, 2], [1, 1, 1]; scale=true, shift=false)
    @test error ≈ 0 atol=1e-12
    @test scale ≈ 1/2 atol=1e-12
    @test shift ≈ 0 atol=1e-12
end

@testitem "Squared error bands" begin
    cryst = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0, 0, 0]])
    qpts = Sunny.QPoints([Sunny.Vec3([0.0, 0.0, 0.0]), Sunny.Vec3([0.0, 0.0, 0.5])])
    disp = [1.0 3.0; 2.0 4.0]
    data = ones(size(disp))
    res = Sunny.BandIntensities(cryst, qpts, disp, data)
    Es = [[2.1], [3.2]]

    @test squared_error_bands(Es, res) ≈ (0.1^2 + 0.2^2) / (2.1^2 + 3.2^2)
    @test Sunny.squared_error_bands_smooth(Es, res; σ=0.1) ≈ 0.0037726373192329475
end

@testitem "Wasserstein1 distance" begin
    @test Sunny.wasserstein1_distance([0.1, -0.2, 0.3], [0.1, -0.2, 0.3]) == 0
    @test Sunny.wasserstein1_distance([0, 1], [0.5, 0]) == 1
    @test Sunny.wasserstein1_distance([0, 1, 0], [0, 0, 1]) == 0.5
end

@testitem "Loss configuration" begin
    cryst = Sunny.square_crystal()
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]), :J1 => 0.0)

    # Something simple to exercise the plumbing
    loss1 = make_loss_fn(sys, [:J1]; hp=(; ϵ=2)) do sys, hp
        return get_param(sys, :J1) + hp.ϵ
    end
    @assert loss1([3.0]) == 5.0

    # A second loss with different hyperparameters
    loss2 = with_hyperparams(loss1, (; ϵ=1))
    @assert loss2([2.0]) == 3.0

    # The first loss should be unchanged
    @assert loss1([3.0]) == 5.0
end
