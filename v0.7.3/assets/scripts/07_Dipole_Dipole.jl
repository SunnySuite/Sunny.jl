using Sunny, GLMakie

units = Units(:K, :angstrom)
latvecs = lattice_vectors(10.19, 10.19, 10.19, 90, 90, 90)
positions = [[0, 0, 0]]
cryst = Crystal(latvecs, positions, 227)
view_crystal(cryst)

sys = System(cryst, [1 => Moment(s=7/2, g=2)], :dipole; seed=0)
sys = reshape_supercell(sys, primitive_cell(cryst))
J1 = 0.304 # (K)
set_exchange!(sys, J1, Bond(1, 2, [0,0,0]))

sys_lr = clone_system(sys)
enable_dipole_dipole!(sys_lr, units.vacuum_permeability)

sys_lr_cut = clone_system(sys)
modify_exchange_with_truncated_dipole_dipole!(sys_lr_cut, 5.0, units.vacuum_permeability)

randomize_spins!(sys_lr)
minimize_energy!(sys_lr)
plot_spins(sys_lr; ghost_radius=8, color=[:red, :blue, :yellow, :purple])

sys.dipoles .= sys_lr.dipoles
sys_lr_cut.dipoles .= sys_lr.dipoles;

qs = [[0,0,0], [0,1,0], [1,1/2,0], [1/2,1/2,1/2], [3/4,3/4,0], [0,0,0]]
labels = ["Γ", "X", "W", "L", "K", "Γ"]
path = q_space_path(cryst, qs, 500; labels)

measure = ssf_trace(sys)
swt = SpinWaveTheory(sys; measure)
res1 = intensities_bands(swt, path)

swt = SpinWaveTheory(sys_lr; measure)
res2 = intensities_bands(swt, path)

swt = SpinWaveTheory(sys_lr_cut; measure)
res3 = intensities_bands(swt, path);

fig = Figure(size=(768, 300))
plot_intensities!(fig[1, 1], res1; units, title="Without long-range dipole")
ax = plot_intensities!(fig[1, 2], res2; units, title="With long-range dipole")
for c in eachrow(res3.disp)
    lines!(ax, eachindex(c), c; linestyle=:dash, color=:black)
end
fig
