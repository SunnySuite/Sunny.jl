using Sunny, GLMakie

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(3, 8, 8, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])
view_crystal(cryst; ndims=2, ghost_radius=8)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
J1 = -1
J2 = +2 * abs(J1)
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J2, Bond(1, 1, [2, 0, 0]))

axis = [0, 0, 1]
randomize_spins!(sys)
k = minimize_spiral_energy!(sys, axis; k_guess=randn(3))

@assert k[1] ≈ 0.2300534561 || k[1] ≈ 1 - 0.2300534561
@assert spiral_energy_per_site(sys; k, axis) ≈ -33/16 * abs(J1)

sys_enlarged = repeat_periodically_as_spiral(sys, (8, 1, 1); k, axis)
plot_spins(sys_enlarged; ndims=2)

swt = SpinWaveTheorySpiral(sys; measure=ssf_perp(sys), k, axis)
qs = [[0,0,0], [1,0,0]]
path = q_space_path(cryst, qs, 401)
res = intensities_bands(swt, path)
plot_intensities(res; units)
