@testitem "Math Basics" begin
    import Sunny: ql_slow
    import LinearAlgebra: istril
    A = Float64[1 2 3; 4 5 6; 7 8 9]
    QF, FRF = ql_slow(A)
    @test istril(FRF)
    @test isapprox(QF * FRF, A)
end