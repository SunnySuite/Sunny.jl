@testset "Fourier Dipole Acceleration" begin

"Tests these field-using functions give the same answer as `ewald_sum_dipole`"
function test_energy_consistency(crystal, latsize)
    μB = Sunny.BOHR_MAGNETON
    μ0 = Sunny.VACUUM_PERM
    sys = SpinSystem(crystal, Sunny.AbstractInteraction[], latsize; μB, μ0)
    rand!(sys)

    dipdip = dipole_dipole(; extent=5, η=0.5)
    dip_real = Sunny.DipoleRealCPU(dipdip, crystal, latsize, sys.site_infos; μB, μ0)
    dip_fourier = Sunny.DipoleFourierCPU(dipdip, crystal, latsize, sys.site_infos; μB, μ0)

    scales = [μB*info.S*info.g for info in sys.site_infos]
    moments = reshape(scales, nbasis(sys), 1, 1, 1) .* sys.sites
    direct_energy = (μ0/4π) * Sunny.ewald_sum_dipole(sys.lattice, moments; extent=5, η=0.5)
    real_energy = Sunny.energy(sys.sites, dip_real)
    fourier_energy = Sunny.energy(sys.sites, dip_fourier)

    @test real_energy ≈ fourier_energy
    @test direct_energy ≈ fourier_energy
end

function test_field_consistency(crystal, latsize)
    sys = SpinSystem(crystal, Sunny.AbstractInteraction[], latsize)
    rand!(sys)
    
    dipdip = dipole_dipole(; extent=4, η=0.5)
    dip_real = Sunny.DipoleRealCPU(dipdip, crystal, latsize, sys.site_infos)
    dip_fourier = Sunny.DipoleFourierCPU(dipdip, crystal, latsize, sys.site_infos)

    H1 = zero(sys)
    H2 = zero(sys)
    Sunny._accum_neggrad!(H1, sys.sites, dip_real)
    Sunny._accum_neggrad!(H2, sys.sites, dip_fourier)

    @test H1 ≈ H2
end

lat_vecs = lattice_vectors(1.0, 1.0, 2.0, 90., 90., 120.)
basis_vecs = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
latsize = [5, 5, 5]
crystal = Crystal(lat_vecs, basis_vecs)

test_energy_consistency(crystal, latsize)
test_field_consistency(crystal, latsize)

end