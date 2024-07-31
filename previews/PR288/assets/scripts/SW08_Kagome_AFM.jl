using Sunny, GLMakie

a = 1.0
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
path = q_space_path(cryst, points_rlu, 400)

swt = SpinWaveTheory(sys; corrspec=ssf_perp(sys))
res = intensities_bands(swt, path)

plot_intensities(res; saturation=0.45)
