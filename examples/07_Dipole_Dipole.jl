# # 7. Long-range dipole interactions
#
# This example shows how long-range dipole-dipole interactions can affect a spin
# wave calculation. These interactions can be included two ways: Ewald summation
# or in real-space with a distance cutoff. The study follows [Del Maestro and
# Gingras, J. Phys.: Cond. Matter, **16**, 3339
# (2004)](https://arxiv.org/abs/cond-mat/0403494).

using Sunny, GLMakie

# Create a pyrochlore crystal from spacegroup 227.

units = Units(:K, :angstrom)
latvecs = lattice_vectors(10.19, 10.19, 10.19, 90, 90, 90)
positions = [[1/8, 1/8, 1/8]]
cryst = Crystal(latvecs, positions, 227, setting="1")
view_crystal(cryst)

# Create a system with antiferromagnetic nearest neighbor exchange.

sys = System(cryst, [1 => Moment(s=7/2, g=2)], :dipole)
J1 = 0.304 # (K)
set_exchange!(sys, J1, Bond(1, 2, [0,0,0]))

# Reshape to the primitive cell, which contains four atoms. To facilitate
# indexing, the function [`position_to_site`](@ref) accepts positions with
# respect to the original (cubic) cell.

shape = [1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2]
sys_prim = reshape_supercell(sys, shape)

set_dipole!(sys_prim, [+1, -1, 0], position_to_site(sys_prim, [1/8, 1/8, 1/8]))
set_dipole!(sys_prim, [-1, +1, 0], position_to_site(sys_prim, [3/8, 3/8, 1/8]))
set_dipole!(sys_prim, [+1, +1, 0], position_to_site(sys_prim, [3/8, 1/8, 3/8]))
set_dipole!(sys_prim, [-1, -1, 0], position_to_site(sys_prim, [1/8, 3/8, 3/8]))

plot_spins(sys_prim; ghost_radius=8, color=[:red, :blue, :yellow, :purple])

# Calculate dispersions with and without long-range dipole interactions. The
# high-symmetry ``𝐪``-points are specified with respect to the conventional
# cubic cell.

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

# Create a panel that qualitatively reproduces Fig. 2 of [Del Maestro and
# Gingras](https://arxiv.org/abs/cond-mat/0403494). That previous work had two
# errors: Its energy scale is too small by a factor of 2 and, in addition,
# slight corrections are needed for the third dispersion band.

fig = Figure(size=(768, 300))
plot_intensities!(fig[1, 1], res1; units, title="Without long-range dipole")
ax = plot_intensities!(fig[1, 2], res2; units, title="With long-range dipole")
for c in eachrow(res3.disp)
    lines!(ax, eachindex(c), c; linestyle=:dash, color=:black)
end
fig
