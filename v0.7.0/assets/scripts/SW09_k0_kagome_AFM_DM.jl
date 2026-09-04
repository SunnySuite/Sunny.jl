using Sunny, GLMakie, LinearAlgebra

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(6, 6, 8, 90, 90, 120)
cryst = Crystal(latvecs, [[1/2, 0, 0]], 147)
view_crystal(cryst; ndims=2)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
J = 1.0*I + dmvec([0, 0, -0.08])

set_exchange!(sys, J, Bond(2, 3, [0, 0, 0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[-1/2,0,0], [0,0,0], [1/2,1/2,0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)

radii = range(0, 2.5, 200)
energies = range(0, 2.5, 200)
kernel = gaussian(fwhm=0.02)
res = powder_average(cryst, radii, 1000) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res; units, colorrange=(0,20))
