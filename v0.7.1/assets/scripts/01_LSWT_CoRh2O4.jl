using Sunny, GLMakie

units = Units(:meV, :angstrom);

a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)

positions = [[0, 0, 0]]
cryst = Crystal(latvecs, positions, 227; types=["Co"], setting="1")

view_crystal(cryst)

sys = System(cryst, [1 => Moment(s=3/2, g=2)], :dipole)

J = +0.63 # (meV)
set_exchange!(sys, J, Bond(1, 3, [0, 0, 0]))
view_crystal(sys)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles])

@assert energy_per_site(sys) ≈ -2J*(3/2)^2

shape = [0 1 1;
         1 0 1;
         1 1 0] / 2
sys_prim = reshape_supercell(sys, shape)
@assert energy_per_site(sys_prim) ≈ -2J*(3/2)^2

plot_spins(sys_prim; color=[S[3] for S in sys_prim.dipoles])

formfactors = [1 => FormFactor("Co2")]
measure = ssf_perp(sys_prim; formfactors)
swt = SpinWaveTheory(sys_prim; measure)

kernel = lorentzian(fwhm=0.8)

qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 500)

energies = range(0, 6, 300)
res = intensities(swt, path; energies, kernel)
plot_intensities(res; units, title="CoRh₂O₄ LSWT")

radii = range(0, 3, 200) # (1/Å)
res = powder_average(cryst, radii, 2000) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res; units, saturation=1.0, title="CoRh₂O₄ Powder Average")
