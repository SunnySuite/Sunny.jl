using Sunny, GLMakie

latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], "P1")

sys = System(cryst, (1, 1, 25), [SpinInfo(1; S=1, g=2)], :dipole; seed=0)
D = 0.1 # meV
B = 5D # ~ 8.64T
set_exchange!(sys, dmvec([0, 0, D]), Bond(1, 1, [0, 0, 1]))
set_field!(sys, [0, 0, B])

randomize_spins!(sys)
minimize_energy!(sys)
@assert energy_per_site(sys) ≈ -10D

dt = 0.1
kT = 0.02
damping = 0.1
langevin = Langevin(dt; kT, damping)
suggest_timestep(sys, langevin; tol=1e-2)
for _ in 1:10_000
    step!(sys, langevin)
end

plot_spins(sys)

sc = dynamical_correlations(sys; dt, nω=100, ωmax=15D)
add_sample!(sc, sys)
nsamples = 100
for _ in 1:nsamples
    for _ in 1:1_000
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end
density = 100
path, xticks = reciprocal_space_path(cryst, [[0,0,-1/2], [0,0,+1/2]], density)
data = intensities_interpolated(sc, path, intensity_formula(sc, :trace; kT));

sys_small = resize_supercell(sys, (1,1,1))
minimize_energy!(sys_small)
swt = SpinWaveTheory(sys_small)
formula = intensity_formula(swt, :trace; kernel=delta_function_kernel)
disp_swt, intens_swt = intensities_bands(swt, path, formula);

fig = Figure()
ax = Axis(fig[1,1]; aspect=1.4, ylabel="ω (meV)", xlabel="𝐪 (r.l.u.)",
          xticks, xticklabelrotation=π/10)
heatmap!(ax, 1:size(data, 1), available_energies(sc), data;  colorrange=(0.0, 50.0))
lines!(ax, disp_swt[:,1]; color=:magenta, linewidth=2)
fig
