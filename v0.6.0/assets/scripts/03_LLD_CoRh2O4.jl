using Sunny, GLMakie, Statistics

a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], 227, setting="1")
latsize = (2, 2, 2)
S = 3/2
J = 0.63 # (meV)
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=0)
set_exchange!(sys, J, Bond(1, 3, [0,0,0]))
randomize_spins!(sys)
minimize_energy!(sys)

sys = resize_supercell(sys, (10, 10, 10))
@assert energy_per_site(sys) ≈ -2J*S^2

kT = 16 * meV_per_K  # 16K, a temperature slightly below ordering
langevin = Langevin(; damping=0.2, kT)

suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.025;

energies = [energy_per_site(sys)]
for _ in 1:1000
    step!(sys, langevin)
    push!(energies, energy_per_site(sys))
end

suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.042;

plot(energies, color=:blue, figure=(size=(600,300),), axis=(xlabel="Timesteps", ylabel="Energy (meV)"))

S_ref = sys.dipoles[1,1,1,1]
plot_spins(sys; color=[s'*S_ref for s in sys.dipoles])

sc = instant_correlations(sys)
add_sample!(sc, sys)    # Accumulate the newly sampled structure factor into `sf`

for _ in 1:20
    for _ in 1:100
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end

q1s = -10:0.1:10
q2s = -10:0.1:10
qs = [[q1, q2, 0.0] for q1 in q1s, q2 in q2s];

formfactors = [FormFactor("Co2")]
instant_formula = intensity_formula(sc, :perp; formfactors)
iq = instant_intensities_interpolated(sc, qs, instant_formula);

heatmap(q1s, q2s, iq;
    colorrange = (0, maximum(iq)/2),
    axis = (
        xlabel="Momentum Transfer Qx (r.l.u)", xlabelsize=16,
        ylabel="Momentum Transfer Qy (r.l.u)", ylabelsize=16,
        aspect=true,
    )
)

dt = 2*langevin.dt
ωmax = 6.0  # Maximum energy to resolve (meV)
nω = 50     # Number of energies to resolve
sc = dynamical_correlations(sys; dt, nω, ωmax)

for _ in 1:5
    for _ in 1:100
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end

points = [[3/4, 3/4,   0],
          [  0,   0,   0],
          [  0, 1/2, 1/2],
          [1/2,   1,   0],
          [  0,   1,   0],
          [1/4,   1, 1/4],
          [  0,   1,   0],
          [  0,  -4,   0]]
density = 50 # (Å)
path, xticks = reciprocal_space_path(cryst, points, density);

formula = intensity_formula(sc, :perp; formfactors, kT=langevin.kT)
fwhm = 0.2
iqw = intensities_interpolated(sc, path, formula)
iqwc = broaden_energy(sc, iqw, (ω, ω₀) -> lorentzian(; fwhm)(ω-ω₀));

ωs = available_energies(sc)
heatmap(1:size(iqwc, 1), ωs, iqwc;
    colorrange = (0, maximum(iqwc)/50),
    axis = (;
        xlabel="Momentum Transfer (r.l.u)",
        ylabel="Energy Transfer (meV)",
        xticks,
        xticklabelrotation=π/5,
        aspect = 1.4,
    )
)

radii = 0:0.05:3.5 # (1/Å)
output = zeros(Float64, length(radii), length(ωs))
for (i, radius) in enumerate(radii)
    pts = reciprocal_space_shell(sc.crystal, radius, 100)
    is = intensities_interpolated(sc, pts, formula)
    is = broaden_energy(sc, is, (ω,ω₀) -> lorentzian(; fwhm)(ω-ω₀))
    output[i, :] = mean(is , dims=1)[1,:]
end

heatmap(radii, ωs, output;
    axis = (
        xlabel="|Q| (Å⁻¹)",
        ylabel="Energy Transfer (meV)",
        aspect = 1.4,
    ),
    colorrange = (0, 20.0)
)
