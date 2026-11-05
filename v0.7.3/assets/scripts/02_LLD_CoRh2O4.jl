using Sunny, GLMakie

units = Units(:meV, :angstrom)
a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)
cryst = Crystal(latvecs, [[1/8, 1/8, 1/8]], 227)

sys = System(cryst, [1 => Moment(s=3/2, g=2)], :dipole)
J = 0.63 # (meV)
set_exchange!(sys, J, Bond(2, 3, [0,0,0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles])

sys = repeat_periodically(sys, (10, 10, 10))
plot_spins(sys; color=[S[3] for S in sys.dipoles])

langevin = Langevin(; damping=0.2, kT=16*units.K)

suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.025;

energies = [energy_per_site(sys)]
for _ in 1:1000
    step!(sys, langevin)
    push!(energies, energy_per_site(sys))
end

suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.042;

lines(energies, color=:blue, figure=(size=(600,300),), axis=(xlabel="Timesteps", ylabel="Energy (meV)"))

S0 = sys.dipoles[1,1,1,1]
plot_spins(sys; color=[S'*S0 for S in sys.dipoles])

formfactors = [1 => FormFactor("Co2")]
measure = ssf_perp(sys; formfactors)
sc = SampledCorrelationsStatic(sys; measure)
add_sample!(sc, sys)    # Accumulate the newly sampled structure factor into `sf`

for _ in 1:20
    for _ in 1:100
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end

grid = q_space_grid(cryst, [1, 0, 0], range(-10, 10, 200), [0, 1, 0], (-10, 10))

res = intensities_static(sc, grid)
plot_intensities(res; saturation=1.0, title="Static Intensities at T = 16 K")

dt = 2*langevin.dt
energies = range(0, 6, 50)
sc = SampledCorrelations(sys; dt, energies, measure)

for _ in 1:5
    for _ in 1:100
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end

qs = [[3/4, 3/4,   0],
      [  0,   0,   0],
      [  0, 1/2, 1/2],
      [1/2,   1,   0],
      [  0,   1,   0],
      [1/4,   1, 1/4],
      [  0,   1,   0],
      [  0,  -4,   0]]
path = q_space_path(cryst, qs, 500)

res = intensities(sc, path; energies, langevin.kT)
plot_intensities(res; units, title="Intensities at 16 K")

radii = range(0, 3.5, 200) # (1/Å)
res = powder_average(cryst, radii, 350) do qs
    intensities(sc, qs; energies, langevin.kT)
end
plot_intensities(res; units, title="Powder Average at 16 K")
