using Sunny, GLMakie

a = b = 4.05012  # Lattice constants for triangular lattice
c = 6.75214      # Spacing in the z-direction

latvecs = lattice_vectors(a, b, c, 90, 90, 120) # A 3x3 matrix of lattice vectors that
                                                # define the conventional unit cell
positions = [[0, 0, 0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]  # Positions of atoms in fractions
                                                           # of lattice vectors
types = ["Fe", "I", "I"]
FeI2 = Crystal(latvecs, positions; types)

cryst = subcrystal(FeI2, "Fe")

view_crystal(cryst)

print_symmetry_table(cryst, 8.0)

sys = System(cryst, (4, 4, 4), [SpinInfo(1, S=1, g=2)], :SUN, seed=2)

J1pm   = -0.236
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

D = 2.165
set_onsite_coupling!(sys, S -> -D*S[3]^2, 1)

randomize_spins!(sys)
minimize_energy!(sys)

plot_spins(sys; color=[s[3] for s in sys.dipoles])

print_wrapped_intensities(sys)

suggest_magnetic_supercell([[0, -1/4, 1/4]])

sys_min = reshape_supercell(sys, [1 0 0; 0 2 1; 0 -2 1])
randomize_spins!(sys_min)
minimize_energy!(sys_min);

plot_spins(sys_min; color=[s[3] for s in sys_min.dipoles], ghost_radius=12)

swt = SpinWaveTheory(sys_min)

q_points = [[0,0,0], [1,0,0], [0,1,0], [1/2,0,0], [0,1,0], [0,0,0]];

density = 50
path, xticks = reciprocal_space_path(cryst, q_points, density);

disp = dispersion(swt, path);

formula = intensity_formula(swt, :perp; kernel=delta_function_kernel)

disp, intensity = intensities_bands(swt, path, formula);

fig = Figure()
ax = Axis(fig[1,1]; xlabel="Momentum (r.l.u.)", ylabel="Energy (meV)", xticks, xticklabelrotation=π/6)
ylims!(ax, 0.0, 7.5)
xlims!(ax, 1, size(disp, 1))
colorrange = extrema(intensity)
for i in axes(disp, 2)
    lines!(ax, 1:length(disp[:,i]), disp[:,i]; color=intensity[:,i], colorrange)
end
fig

fwhm = 0.3 # full width at half maximum (meV)
broadened_formula = intensity_formula(swt, :perp; kernel=lorentzian(; fwhm))

energies = collect(0:0.01:10)  # 0 < ω < 10 (meV).
is1 = intensities_broadened(swt, path, energies, broadened_formula);

R = rotation_in_rlu(cryst, [0, 0, 1], 2π/3)
is2 = intensities_broadened(swt, [R*q for q in path], energies, broadened_formula)
is3 = intensities_broadened(swt, [R*R*q for q in path], energies, broadened_formula)
is_averaged = (is1 + is2 + is3) / 3

fig = Figure()
ax = Axis(fig[1,1]; xlabel="Momentum (r.l.u.)", ylabel="Energy (meV)", xticks, xticklabelrotation=π/6)
heatmap!(ax, 1:size(is_averaged, 1), energies, is_averaged)
fig

disp, is = dssf(swt, path);
