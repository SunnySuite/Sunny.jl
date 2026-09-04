using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.5"

units = Units(:meV, :angstrom)
a = 3.0
b = 8.0
c = 8.0
latvecs = lattice_vectors(a, b, c, 90, 90, 90)
positions = [[0, 0, 0]]
cryst = Crystal(latvecs, positions)

view_crystal(cryst; ndims=2, ghost_radius=8)

print_symmetry_table(cryst, 8)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)

J = -1
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)
energy_per_site(sys)

plot_spins(sys; ndims=2, ghost_radius=8)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))

qs = [[0,0,0], [1,0,0]]
path = q_space_path(cryst, qs, 400)

res = intensities_bands(swt, path)
plot_intensities(res; units)

radii = range(0, 2.5, 200) # 1/Å
energies = range(0, 5, 200) # meV
kernel = gaussian(fwhm=0.1)
res = powder_average(cryst, radii, 1000) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res; units)
