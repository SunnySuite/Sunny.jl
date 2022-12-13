@testitem "Fourier Dipole Acceleration" begin

"Tests these field-using functions give the same answer as `ewald_sum_dipole`"
function test_energy_consistency(crystal, latsize)
    sys = SpinSystem(crystal, Sunny.AbstractInteraction[], latsize)
    rand!(sys)
    (; μ0, μB) = sys.consts

    dip_real = Sunny.DipoleRealCPU(sys.crystal, sys.size, sys.site_infos, sys.consts; extent=5, η=0.5)
    dip_fourier = Sunny.DipoleFourierCPU(sys.crystal, sys.size, sys.site_infos, sys.consts; extent=5, η=0.5)

    scales = [μB*info.spin_rescaling*info.g for info in sys.site_infos]
    moments = reshape(scales, 1, 1, 1, nbasis(sys.crystal)) .* sys.dipoles
    direct_energy = (μ0/4π) * Sunny.ewald_sum_dipole(sys.crystal, moments; extent=5, η=0.5)
    real_energy = Sunny.energy(sys.dipoles, dip_real)
    fourier_energy = Sunny.energy(sys.dipoles, dip_fourier)

    @test real_energy ≈ fourier_energy
    @test direct_energy ≈ fourier_energy
end

function test_field_consistency(crystal, latsize)
    sys = SpinSystem(crystal, Sunny.AbstractInteraction[], latsize)
    rand!(sys)
    
    dip_real = Sunny.DipoleRealCPU(sys.crystal, sys.size, sys.site_infos, sys.consts; extent=4, η=0.5)
    dip_fourier = Sunny.DipoleFourierCPU(sys.crystal, sys.size, sys.site_infos, sys.consts; extent=4, η=0.5)

    H1 = zero(sys.dipoles)
    H2 = zero(sys.dipoles)
    Sunny._accum_neggrad!(H1, sys.dipoles, dip_real)
    Sunny._accum_neggrad!(H2, sys.dipoles, dip_fourier)

    @test H1 ≈ H2
end

lat_vecs = lattice_vectors(1.0, 1.0, 2.0, 90., 90., 120.)
basis_vecs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
latsize = (5, 5, 5)
crystal = Crystal(lat_vecs, basis_vecs)

test_energy_consistency(crystal, latsize)
test_field_consistency(crystal, latsize)

end