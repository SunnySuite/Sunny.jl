using Sunny, GLMakie

units = Units(:K, :angstrom)
latvecs = lattice_vectors(10.19, 10.19, 10.19, 90, 90, 90)
positions = [[1/8, 1/8, 1/8]]
cryst = Crystal(latvecs, positions, 227, setting="1")
view_crystal(cryst)

sys = System(cryst, [1 => Moment(s=7/2, g=2)], :dipole)
J1 = 0.304 # (K)
set_exchange!(sys, J1, Bond(1, 2, [0,0,0]))

shape = [1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2]
sys_prim = reshape_supercell(sys, shape)

set_dipole!(sys_prim, [+1, -1, 0], position_to_site(sys_prim, [1/8, 1/8, 1/8]))
set_dipole!(sys_prim, [-1, +1, 0], position_to_site(sys_prim, [3/8, 3/8, 1/8]))
set_dipole!(sys_prim, [+1, +1, 0], position_to_site(sys_prim, [3/8, 1/8, 3/8]))
set_dipole!(sys_prim, [-1, -1, 0], position_to_site(sys_prim, [1/8, 3/8, 3/8]))

plot_spins(sys_prim; ghost_radius=8, color=[:red, :blue, :yellow, :purple])

qs = [[0,0,0], [0,1,0], [1,1/2,0], [1/2,1/2,1/2], [3/4,3/4,0], [0,0,0]]
labels = ["Γ", "X", "W", "L", "K", "Γ"]
path = q_space_path(cryst, qs, 500; labels)

measure = ssf_trace(sys_prim)
swt = SpinWaveTheory(sys_prim; measure)
res1 = intensities_bands(swt, path)

sys_prim_dd = clone_system(sys_prim)
enable_dipole_dipole!(sys_prim_dd, units.vacuum_permeability)
swt = SpinWaveTheory(sys_prim_dd; measure)
res2 = intensities_bands(swt, path)

sys_prim_tdd = clone_system(sys_prim)
modify_exchange_with_truncated_dipole_dipole!(sys_prim_tdd, 5.0, units.vacuum_permeability)
swt = SpinWaveTheory(sys_prim_tdd; measure)
res3 = intensities_bands(swt, path)

fig = Figure(size=(768, 300))
plot_intensities!(fig[1, 1], res1; units, title="Local Exchange Only")
ax = plot_intensities!(fig[1, 2], res2; units, title="Local Exchange and Dipole-Dipole")
for c in eachrow(res3.disp)
    lines!(ax, eachindex(c), c; linestyle=:dash, color=:black)
end
fig
