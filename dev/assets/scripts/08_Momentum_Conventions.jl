using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.8.0"

latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], "P1")

s = 3/2
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
J = 1.0
D = 0.2
z = [0, 0, 1]
set_exchange!(sys, dmvec(D * z) - J * z * z', Bond(1, 1, [0, 0, 1]))

polarize_spins!(sys, [0, 0, 1])
@assert energy(sys) ≈ - s^2
plot_spins(sys)

path = q_space_path(cryst, [[0, 0, -1/2], [0, 0, 0], [0, 0, +1/2]], 400)
swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
res = intensities_bands(swt, path)
plot_intensities(res; ylims=(0, 5))

sys2 = resize_supercell(sys, (1, 1, 32))

dt = 0.03 / J
kT = 0.03 * J
damping = 0.1
langevin = Langevin(dt; kT, damping)
suggest_timestep(sys2, langevin; tol=1e-2)
for _ in 1:10_000
    step!(sys, langevin)
end

sc = SampledCorrelations(sys2; dt, energies=range(0, 5, 100), measure=ssf_trace(sys2))
add_sample!(sc, sys2)
nsamples = 100
for _ in 1:nsamples
    for _ in 1:1000
        step!(sys2, langevin)
    end
    add_sample!(sc, sys2)
end

res2 = intensities(sc, path; energies=:available, kT)
plot_intensities(res2)
