using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.5"

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(3.0, 3.0, 6.0, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])
view_crystal(cryst; ndims=2)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2, 2, 1))
J1 = 1.0
J2 = -0.1
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J2, Bond(1, 1, [1, 1, 0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)
