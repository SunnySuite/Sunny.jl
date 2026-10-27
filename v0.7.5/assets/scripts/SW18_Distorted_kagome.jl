using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.5"

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(10.2, 5.94, 7.81, 90, 117.7, 90)
positions = [[0, 0, 0], [1/4, 1/4, 0]]
types = ["Cu1", "Cu2"]
cryst = Crystal(latvecs, positions, "C 1 2/m 1"; types)
view_crystal(cryst)

moments = [1 => Moment(s=1/2, g=2), 3 => Moment(s=1/2, g=2)]
sys = System(cryst, moments, :dipole)
J   = -2
Jp  = -1
Jab = 0.75
Ja  = -J/.66 - Jab
Jip = 0.01
set_exchange!(sys, J, Bond(1, 3, [0, 0, 0]))
set_exchange!(sys, Jp, Bond(3, 5, [0, 0, 0]))
set_exchange!(sys, Ja, Bond(3, 4, [0, 0, 0]))
set_exchange!(sys, Jab, Bond(1, 2, [0, 0, 0]))
set_exchange!(sys, Jip, Bond(3, 4, [0, 0, 1]))

axis = [0, 0, 1]
randomize_spins!(sys)
k = minimize_spiral_energy!(sys, axis; k_guess=randn(3))
plot_spins(sys; ndims=2)

k_ref = [0.785902495, 0.0, 0.107048756]
k_ref_alt = [1, 0, 1] - k_ref
@assert isapprox(k, k_ref; atol=1e-6) || isapprox(k, k_ref_alt; atol=1e-6)
@assert spiral_energy_per_site(sys; k, axis) ≈ -0.78338383838

suggest_magnetic_supercell([k_ref]; tol=1e-3)

new_shape = [14 0 1; 0 1 0; 0 0 2]
sys2 = reshape_supercell(sys, new_shape)
randomize_spins!(sys2)
minimize_energy!(sys2)
energy_per_site(sys2)

measure = ssf_perp(sys; apply_g=false)
swt = SpinWaveTheorySpiral(sys; measure, k, axis)

qs = [[0,0,0], [1,0,0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)

radii = range(0, 2, 100) # (1/Å)
energies = range(0, 6, 200)
kernel = gaussian(fwhm=0.05)
res = powder_average(cryst, radii, 400) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res; units)
