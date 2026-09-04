using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.5"

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(3, 8, 8, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])
view_crystal(cryst; ndims=2, ghost_radius=8)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2, 1, 1))

J = 1
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2, ghost_radius=8)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[0,0,0], [1,0,0]]
path = q_space_path(cryst, qs, 401)
res = intensities_bands(swt, path)
plot_intensities(res; units)

isapprox(res.disp[1, :], res.disp[2, :])

xs = [q[1] for q in path.qs]
ys = log10.(res.data[1, :] + res.data[2, :])
lines(xs, ys; axis=(; xlabel="[H, 0, 0]", ylabel="Log intensity (dimensionless)"))
