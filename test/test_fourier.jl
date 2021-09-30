"Tests these field-using functions give the same answer as `ewald_sum_dipole`"
function test_energy_consistency(sys::SpinSystem{3})
    dip_real = DipoleRealCPU(1.0, sys; extent=4, η=0.5)
    dip_fourier = DipoleFourierCPU(1.0, sys; extent=4, η=0.5)

    direct_energy = ewald_sum_dipole(sys; extent=4, η=0.5)
    real_energy = energy(sys, dip_real)
    fourier_energy = energy(sys, dip_fourier)

    @assert direct_energy ≈ real_energy      "`DipoleRealPre` energy not correct!"
    @assert direct_energy ≈ fourier_energy   "`DipoleFourier` energy not correct!"
end

function test_field_consistency(sys::SpinSystem{3})
    cryst = Crystal(sys.lattice)
    dipdip = DipoleDipole(1.0; extent=4, η=0.5)
    dip_real = DipoleRealCPU(dipdip, cryst, sys.lattice.size)
    dip_fourier = DipoleFourierCPU(dipdip, cryst, sys.lattice.size)

    H1 = zero(sys)
    H2 = zero(sys)
    _accum_field!(H1, sys.sites, dip_real)
    _accum_field!(H2, sys.sites, dip_fourier)

    @assert all(H1 .≈ H2)
end