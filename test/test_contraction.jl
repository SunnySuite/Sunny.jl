@testitem "Dipole Factor Ordering" begin
    using LinearAlgebra
    obs = Sunny.parse_observables(3; observables=nothing, correlations=nothing)
    dipoleinfo = Sunny.DipoleFactor(obs;unilateral_to_bilateral = false)
    fake_intensities = [1. 3 5; 17 7 11; 19 23 13]
    # This is the order expected by contract(...)
    fake_dipole_elements = fake_intensities[:]
    k = Sunny.Vec3(1,2,3)
    mat = Sunny.polarization_matrix(k)
    out = Sunny.contract(fake_dipole_elements,k,dipoleinfo)
    out_true = dot(fake_intensities, mat)
    @test isapprox(out,out_true)
end
