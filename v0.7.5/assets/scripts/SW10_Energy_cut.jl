using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.5"

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(1.0, 1.0, 3.0, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2, 2, 1))
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)

grid = q_space_grid(cryst, [1, 0, 0], range(0, 2, 201), [0, 1, 0], range(0, 2, 201))

swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
res = intensities(swt, grid; energies=[3.75], kernel=gaussian(fwhm=0.2))
plot_intensities(res; units)

res = intensities_static(swt, grid; bounds=(3.5, 4.01))
plot_intensities(res; units)
