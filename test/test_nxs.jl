@testitem "load_nxs" begin
    import CodecZlib, IOCapture

    filename = joinpath(@__DIR__, "nxs", "1Dcut_sym_0p00T_K.nxs")

    # Discard warnings
    capt = IOCapture.capture() do
        load_nxs(filename)
    end
    params, signal = capt.value

    @test isapprox(params.binstart, [-0.05, 0.617, -0.25, -0.1]; atol=1e-7)
    @test isapprox(params.binend, [0.0, 0.667, 0.0, 0.6]; atol=1e-7)
    @test isapprox(params.binwidth, [0.1, 0.1, 0.5, 0.008]; atol=1e-7)
    @test params.covectors == [ 1.0  0.5  0.0  0.0; 0.0  1.0  0.0  0.0; 0.0  0.0  1.0  0.0; 0.0  0.0  0.0  1.0]
    @test signal[5:10] ≈ [0.0002707771, 0.0001276792, 0.0002237122, 0.0015903557, 0.0177874433, 0.0950909946]
end
