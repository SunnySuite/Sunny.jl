@testitem "Squared error fitted" begin
    @test squared_error([1, 3], [1, 1]; normalize=false) ≈ 4 atol=1e-12
    @test squared_error([1, 3], [1, 1]) ≈ 2/5 atol=1e-12
    @test squared_error([1, -1], [-1, 1]) ≈ 4 atol=1e-12

    r = squared_error_fitted([4, 4], [2, -1]; scale=true, weights=[0.5, 0.5])
    @test r.error ≈ 0.9 atol=1e-12
    @test r.scale ≈ 0.8 atol=1e-12

    r = squared_error_fitted([4, 4], [2, -1]; scale=true, weights=[0.5, 0.5], normalize=false)
    @test r.error ≈ 14.4 atol=1e-12
    @test r.scale ≈ 0.8 atol=1e-12

    (; error, scale, shift) = squared_error_fitted([1, 2, 3], [1, 1, 1]; shift=true)
    @test error ≈ 1 atol=1e-12
    @test scale ≈ 1 atol=1e-12
    @test shift ≈ 1 atol=1e-12

    (; error, scale, shift) = squared_error_fitted([1, 2, 3], [1, 1, 1]; shift=true, normalize=false)
    @test error ≈ 2 atol=1e-12
    @test scale ≈ 1 atol=1e-12
    @test shift ≈ 1 atol=1e-12

    (; error, scale, shift) = squared_error_fitted([3, 4, 5], [1, 1, 2]; scale=true, shift=true)
    @test error ≈ 1/4 atol=1e-12
    @test scale ≈ 3/2 atol=1e-12
    @test shift ≈ 2 atol=1e-12

    # Complex scaling factor absorbs phase
    u = ComplexF64[1 + 2im, -3 + im, 2 - 4im]
    α = 2 - 3im
    (; error, scale) = squared_error_fitted(u, u/α; scale=true)
    @test error ≈ 0 atol=1e-12
    @test scale ≈ α atol=1e-12
end

@testitem "Squared error bands" begin
    cryst = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0, 0, 0]])
    qpts = convert(Sunny.AbstractQPoints, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]])
    disp = [1.0 3.0; 2.0 4.0]
    data = ones(size(disp))
    res = Sunny.BandIntensities(cryst, qpts, disp, data)
    Es = [[2.1], [3.2]]

    # Smooth variant
    @test squared_error_bands(Es, res) ≈ (0.1^2 + 0.2^2) / (2.1^2 + 3.2^2)
    @test Sunny.squared_error_bands_smooth(Es, res; σ=0.1) ≈ 0.0037726373192329475

    # Ignore bands below intensity cutoff
    qpts = convert(Sunny.AbstractQPoints, [[0.0, 0.0, 0.0]])
    disp = [1.0, 2.0]
    data = [1.0im, 1e-6im]
    res = Sunny.BandIntensities(cryst, qpts, disp, data)
    Es = [[2.1]]
    weights = [[1.5]]
    @test squared_error_bands(Es, res; weights, intensity_cutoff=0.1, normalize=false) ≈ 1.5 * 1.1^2
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
