using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.6"

units = Units(:meV, :angstrom)
a = b = 4.05012  # Lattice constants for triangular lattice (Å)
c = 6.75214      # Spacing between layers (Å)
latvecs = lattice_vectors(a, b, c, 90, 90, 120)
positions = [[0, 0, 0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]
types = ["Fe", "I", "I"]
cryst = Crystal(latvecs, positions; types)

cryst = subcrystal(cryst, "Fe")
view_crystal(cryst)

print_symmetry_table(cryst, 8.0)

sys = System(cryst, [1 => Moment(s=1, g=2)], :SUN; seed=2)

J1pm   = -0.236 # (meV)
J1pmpm = -0.161
J1zpm  = -0.261
J2pm   = 0.026
J3pm   = 0.166
J′0pm  = 0.037
J′1pm  = 0.013
J′2apm = 0.068

J1zz   = -0.236
J2zz   = 0.113
J3zz   = 0.211
J′0zz  = -0.036
J′1zz  = 0.051
J′2azz = 0.073

J1xx = J1pm + J1pmpm
J1yy = J1pm - J1pmpm
J1yz = J1zpm

set_exchange!(sys, [J1xx   0.0    0.0;
                    0.0    J1yy   J1yz;
                    0.0    J1yz   J1zz], Bond(1,1,[1,0,0]))
set_exchange!(sys, [J2pm   0.0    0.0;
                    0.0    J2pm   0.0;
                    0.0    0.0    J2zz], Bond(1,1,[1,2,0]))
set_exchange!(sys, [J3pm   0.0    0.0;
                    0.0    J3pm   0.0;
                    0.0    0.0    J3zz], Bond(1,1,[2,0,0]))
set_exchange!(sys, [J′0pm  0.0    0.0;
                    0.0    J′0pm  0.0;
                    0.0    0.0    J′0zz], Bond(1,1,[0,0,1]))
set_exchange!(sys, [J′1pm  0.0    0.0;
                    0.0    J′1pm  0.0;
                    0.0    0.0    J′1zz], Bond(1,1,[1,0,1]))
set_exchange!(sys, [J′2apm 0.0    0.0;
                    0.0    J′2apm 0.0;
                    0.0    0.0    J′2azz], Bond(1,1,[1,2,1]))

D = 2.165 # (meV)
set_onsite_coupling!(sys, S -> -D*S[3]^2, 1)

sys = resize_supercell(sys, (4, 4, 4))
randomize_spins!(sys)
minimize_energy!(sys)

plot_spins(sys; color=[S[3] for S in sys.dipoles])

print_wrapped_intensities(sys)

suggest_magnetic_supercell([[0, -1/4, 1/4]])

sys_min = reshape_supercell(sys, [1 0 0; 0 2 1; 0 -2 1])
randomize_spins!(sys_min)
minimize_energy!(sys_min);
plot_spins(sys_min; color=[S[3] for S in sys_min.dipoles], ghost_radius=12)

swt = SpinWaveTheory(sys_min; measure=ssf_perp(sys_min))

qs = [[0,0,0], [1,0,0], [0,1,0], [1/2,0,0], [0,1,0], [0,0,0]]
path = q_space_path(cryst, qs, 500)
res = intensities_bands(swt, path)
plot_intensities(res; units, ylims=(0, 10), title="Single Crystal Bands")

kernel = lorentzian(fwhm=0.3)
energies = range(0, 10, 300);  # 0 < ω < 10 (meV)

rotations = [([0,0,1], n*(2π/3)) for n in 0:2]
weights = [1, 1, 1]
res = domain_average(cryst, path; rotations, weights) do path_rotated
    intensities(swt, path_rotated; energies, kernel)
end
plot_intensities(res; units, colormap=:viridis, title="Domain Averaged Intensities")
