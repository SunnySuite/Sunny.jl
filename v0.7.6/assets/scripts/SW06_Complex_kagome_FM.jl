using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.6"

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(6, 6, 8, 90, 90, 120)
cryst = Crystal(latvecs, [[1/2, 0, 0]], 147)
view_crystal(cryst; ndims=2)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
J1 = -1.0
J2 = 0.1
J3a = 0.00
J3b = 0.17
set_exchange!(sys, J1, Bond(2, 3, [0, 0, 0]))
set_exchange!(sys, J2, Bond(2, 1, [0, 0, 0]))
set_exchange!(sys, J3a, Bond(2, 2, [1, 0, 0]))
set_exchange!(sys, J3b, Bond(1, 1, [1, 0, 0]))

view_crystal(sys; ndims=2)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[-1/2,0,0], [0,0,0], [1/2,1/2,0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)

radii = range(0, 2.5, 200)
energies = range(0, 6.5, 200)
kernel = gaussian(fwhm=0.02)
res = powder_average(cryst, radii, 1000) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res; units, colorrange=(0,10))
