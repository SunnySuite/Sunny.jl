using Sunny, GLMakie

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(6, 6, 5, 90, 90, 120)
positions = [[1/2, 0, 0]]
cryst = Crystal(latvecs, positions, 147)

positions = [[1/2, 0, 0], [0, 1/2, 0], [1/2, 1/2, 0]]
cryst2 = Crystal(latvecs, positions)

view_crystal(cryst; ndims=2)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
J = -1.0
set_exchange!(sys, J, Bond(2, 3, [0, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)
energy_per_site(sys)
@assert energy_per_site(sys) ≈ 4J/2
plot_spins(sys; ndims=2)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[-1/2, 0, 0], [0, 0, 0], [1/2, 1/2, 0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)

radii = range(0, 2.5, 200)
energies = range(0, 6.5, 200)
res1 = powder_average(cryst, radii, 1000) do qs
    intensities(swt, qs; energies, kernel=gaussian(fwhm=0.02))
end
res2 = powder_average(cryst, radii, 1000) do qs
    intensities(swt, qs; energies, kernel=gaussian(fwhm=0.25))
end

fig = Figure(size=(768, 800))
plot_intensities!(fig[1, 1], res1; units, colorrange=(0,10), title="FWHM 0.02 meV")
plot_intensities!(fig[2, 1], res2; units, colorrange=(0,10), title="FWHM 0.25 meV")
fig
