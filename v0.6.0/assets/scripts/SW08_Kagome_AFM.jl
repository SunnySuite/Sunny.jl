using Sunny, GLMakie

a = 1
latvecs = lattice_vectors(a, a, 10a, 90, 90, 120)
cryst = Crystal(latvecs, [[1/2,0,0]], 147)

S = 1
sys = System(cryst, (3, 3, 1), [SpinInfo(1; S, g=2)], :dipole)
J = 1.0
set_exchange!(sys, J, Bond(2, 3, [0, 0, 0]))

k = -[1/3, 1/3, 0]
axis = [0,0,1]
set_spiral_order_on_sublattice!(sys, 1; k, axis, S0=[cos(0),sin(0),0])
set_spiral_order_on_sublattice!(sys, 2; k, axis, S0=[cos(0),sin(0),0])
set_spiral_order_on_sublattice!(sys, 3; k, axis, S0=[cos(2π/3),sin(2π/3),0])
plot_spins(sys; dims=2)

@assert energy_per_site(sys) ≈ (4/2)*J*S^2*cos(2π/3)

points_rlu = [[-1/2, 0, 0], [0, 0, 0], [1/2, 1/2, 0]]
density = 100
path, xticks = reciprocal_space_path(cryst, points_rlu, density);

swt = SpinWaveTheory(sys)
formula = intensity_formula(swt, :perp; kernel=delta_function_kernel)
disp, intensity = intensities_bands(swt, path, formula);

fig = Figure()
ax = Axis(fig[1,1]; xlabel="Momentum (r.l.u.)", ylabel="Energy (meV)",
          xticks, xticklabelrotation=π/6)
ylims!(ax, -1e-1, 2.3)
for i in axes(disp, 2)
    lines!(ax, 1:length(disp[:,i]), disp[:,i]; color=intensity[:,i], colorrange=(0,1e-2))
end
fig
