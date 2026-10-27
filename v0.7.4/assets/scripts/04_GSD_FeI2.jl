using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.4"

units = Units(:meV, :angstrom)
a = b = 4.05012
c = 6.75214
latvecs = lattice_vectors(a, b, c, 90, 90, 120)
cryst = Crystal(latvecs, [[0,0,0]], 164; types=["Fe"])

sys = System(cryst, [1 => Moment(s=1, g=2)], :SUN)
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
set_exchange!(sys, [J1xx 0.0 0.0; 0.0 J1yy J1yz; 0.0 J1yz J1zz], Bond(1,1,[1,0,0]))
set_exchange!(sys, [J2pm 0.0 0.0; 0.0 J2pm 0.0; 0.0 0.0 J2zz], Bond(1,1,[1,2,0]))
set_exchange!(sys, [J3pm 0.0 0.0; 0.0 J3pm 0.0; 0.0 0.0 J3zz], Bond(1,1,[2,0,0]))
set_exchange!(sys, [J′0pm 0.0 0.0; 0.0 J′0pm 0.0; 0.0 0.0 J′0zz], Bond(1,1,[0,0,1]))
set_exchange!(sys, [J′1pm 0.0 0.0; 0.0 J′1pm 0.0; 0.0 0.0 J′1zz], Bond(1,1,[1,0,1]))
set_exchange!(sys, [J′2apm 0.0 0.0; 0.0 J′2apm 0.0; 0.0 0.0 J′2azz], Bond(1,1,[1,2,1]))
D = 2.165
set_onsite_coupling!(sys, S -> -D*S[3]^2, 1)

sys = resize_supercell(sys, (16, 16, 4))

langevin = Langevin(; damping=0.2, kT=2.3*units.K)

randomize_spins!(sys)
minimize_energy!(sys; maxiters=10)
suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.03;

for _ in 1:10_000
    step!(sys, langevin)
end
plot_spins(sys; color=[S[3] for S in sys.dipoles])

print_wrapped_intensities(sys)

suggest_timestep(sys, langevin; tol=1e-2)

langevin.kT = 3.5 * units.K
for _ in 1:10_000
    step!(sys, langevin)
end

suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.040;

dt = 2*langevin.dt
energies = range(0, 7.5, 120)
formfactors = [1 => FormFactor("Fe2"; g_lande=3/2)]
sc = SampledCorrelations(sys; dt, energies, measure=ssf_perp(sys; formfactors))

add_sample!(sc, sys)

for _ in 1:2
    for _ in 1:1000               # Enough steps to decorrelate spins
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end

res = intensities(sc, [[0, 0, 0], [0.5, 0.5, 0.5]]; energies, langevin.kT)
fig = lines(res.energies, res.data[:, 1]; axis=(xlabel="meV", ylabel="Intensity"), label="(0,0,0)")
lines!(res.energies, res.data[:, 2]; label="(π,π,π)")
axislegend()
fig

qs = [[0,   0, 0],  # List of wave vectors that define a path
      [1,   0, 0],
      [0,   1, 0],
      [1/2, 0, 0],
      [0,   1, 0],
      [0,   0, 0]]
qpath = q_space_path(cryst, qs, 500)
res = intensities(sc, qpath; energies, langevin.kT)
plot_intensities(res; colorrange=(0.0, 1.0), title="Intensities at T = 2.3 K")

grid = q_space_grid(cryst, [1, 0, 0], range(-1.5, 1.5, 300), [0, 1, 0], (-1.5, 1.5); orthogonalize=true)
res = intensities(sc, grid; energies=[3.5], langevin.kT)
plot_intensities(res; title="Intensity slice at ω = 3.5 meV")
