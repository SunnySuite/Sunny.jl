using Sunny, GLMakie

units = Units(:meV)
a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], 227, setting="1")

view_crystal(cryst)

latsize = (1, 1, 1)
S = 3/2
J = 0.63 # (meV)
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=0)
set_exchange!(sys, J, Bond(1, 3, [0,0,0]))

randomize_spins!(sys)
minimize_energy!(sys)

@assert energy_per_site(sys) ≈ -2J*S^2

s0 = sys.dipoles[1,1,1,1]
plot_spins(sys; color=[s'*s0 for s in sys.dipoles])

shape = [0 1 1;
         1 0 1;
         1 1 0] / 2
sys_prim = reshape_supercell(sys, shape)
@assert energy_per_site(sys_prim) ≈ -2J*S^2
plot_spins(sys_prim; color=[s'*s0 for s in sys_prim.dipoles])

swt = SpinWaveTheory(sys_prim; measure=ssf_perp(sys_prim))

qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 400)

kernel = Sunny.lorentzian2(fwhm=0.8)
formfactors = [FormFactor("Co2")]
energies = range(0, 6, 300)
res = intensities(swt, path; energies, kernel, formfactors)
plot_intensities(res; units)

radii = range(0, 3, 200) # (1/Å)
res = powder_average(cryst, radii, 2000) do qs
    intensities(swt, qs; energies, kernel, formfactors)
end
plot_intensities(res; units)
