using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.4"

latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], "P1")

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(1, 1, 25))
D = 0.1
B = 5D
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

sc = SampledCorrelations(sys; dt, energies=range(0, 15D, 100), measure=ssf_trace(sys))
add_sample!(sc, sys)
nsamples = 100
for _ in 1:nsamples
    for _ in 1:1_000
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end
path = q_space_path(cryst, [[0,0,-1/2], [0,0,+1/2]], 400)
res1 = intensities(sc, path; energies=:available, kT)

sys_small = resize_supercell(sys, (1,1,1))
minimize_energy!(sys_small)
swt = SpinWaveTheory(sys_small; measure=ssf_trace(sys_small))
res2 = intensities_bands(swt, path)

fig = Figure(size=(768, 300))
plot_intensities!(fig[1, 1], res1; title="Classical dynamics")
plot_intensities!(fig[1, 2], res2; title="Spin wave theory")
fig
