# # 3. Landau-Lifshitz dynamics of CoRh₂O₄ at finite *T*

using Sunny, GLMakie, Statistics

# ### System construction

# Construct the system as in the previous [CoRh₂O₄ tutorial](@ref "2. Spin wave
# simulations of CoRh₂O₄"). After optimization, the system will be in an
# unfrustrated antiferromagnetic ground state.

a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], 227; setting="1")
latsize = (2, 2, 2)
S = 3 / 2
J = 0.63 # (meV)
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=0)
set_exchange!(sys, J, Bond(1, 3, [0, 0, 0]))
randomize_spins!(sys)
minimize_energy!(sys)

# Use [`resize_supercell`](@ref) to build a new system with a lattice of
# 10×10×10 unit cells. The desired Néel order is retained.

sys = resize_supercell(sys, (10, 10, 10))
@assert energy_per_site(sys) ≈ -2J * S^2

# Use the stochastic Landau-Lifshitz dynamics to thermalize system into
# equilibrium at finite temperature. This is a [`Langevin`](@ref) equation,
# which includes damping and noise terms. The strength of the noise term is
# automatically fixed according to the damping time scale `λ` and the target
# temperature, according to a fluctuation-dissipation theorem.

Δt = 0.05 / abs(J * S)   # Time step
λ = 0.1             # Dimensionless damping time-scale
kT = 16 * meV_per_K  # 16K, a temperature slightly below ordering
langevin = Langevin(Δt; λ, kT);

# Because the magnetic order has been initialized correctly, relatively few
# additional Langevin time-steps are required to reach thermal equilibrium.

energies = [energy_per_site(sys)]
for _ in 1:1000
    step!(sys, langevin)
    push!(energies, energy_per_site(sys))
end
plot(
    energies;
    color=:blue,
    figure=(size=(600, 300),),
    axis=(xlabel="Time steps", ylabel="Energy (meV)"),
)

# Thermal fluctuations are apparent in the spin configuration.

S_ref = sys.dipoles[1, 1, 1, 1]
plot_spins(sys; color=[s' * S_ref for s in sys.dipoles])

# ### Instantaneous structure factor

# To visualize the instantaneous (equal-time) structure factor, create an object
# [`instant_correlations`](@ref) and use [`add_sample!`](@ref) to accumulated
# data for the equilibrated spin configuration.

sc = instant_correlations(sys)
add_sample!(sc, sys)    # Accumulate the newly sampled structure factor into `sf`

# Collect 20 additional decorrelated samples. For each sample, about 100
# Langevin time-steps is sufficient to collect approximately uncorrelated
# statistics.

for _ in 1:20
    for _ in 1:100
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end

# Define a slice of momentum space. Wavevectors are specified in reciprocal
# lattice units (RLU). The notation `q1s in -10:0.1:10` indicates that the first
# ``q``-component ranges from -10 to 10 in intervals of 0.1. That is, ``q``
# spans over 20 Brillouin zones. To convert to absolute momentum units, each
# component of ``q`` would need to be scaled by a reciprocal lattice vector.

q1s = -10:0.1:10
q2s = -10:0.1:10
qs = [[q1, q2, 0.0] for q1 in q1s, q2 in q2s];

# Plot the instantaneous structure factor for the given ``q``-slice. We employ
# the appropriate [`FormFactor`](@ref) for Co2⁺. An [`intensity_formula`](@ref)
# defines how dynamical correlations correspond to the observable structure
# factor. The function [`instant_intensities_interpolated`](@ref) calculates
# intensities at the target `qs` by interpolating over the data available at
# discrete reciprocal-space lattice points.

formfactors = [FormFactor("Co2")]
instant_formula = intensity_formula(sc, :perp; formfactors)
iq = instant_intensities_interpolated(sc, qs, instant_formula);

# Plot the resulting intensity data ``I(𝐪)``. The color scale is clipped to 50%
# of the maximum intensity.

heatmap(
    q1s,
    q2s,
    iq;
    colorrange=(0, maximum(iq) / 2),
    axis=(
        xlabel="Momentum Transfer Qx (r.l.u)",
        xlabelsize=16,
        ylabel="Momentum Transfer Qy (r.l.u)",
        ylabelsize=16,
        aspect=true,
    ),
)

# ### Dynamical structure factor

# To collect statistics for the dynamical structure factor intensities
# ``I(𝐪,ω)`` at finite temperature, use [`dynamical_correlations`](@ref). Now,
# each call to `add_sample!` will run a classical spin dynamics trajectory.
# Longer-time trajectories will be required to achieve greater energy
# resolution, as controlled by `nω`. Here, we pick a moderate number of
# energies, `nω = 50`, which will make the simulation run quickly.

ωmax = 6.0  # Maximum  energy to resolve (meV)
nω = 50     # Number of energies to resolve
sc = dynamical_correlations(sys; Δt, nω, ωmax, process_trajectory=:symmetrize)

# Each sample requires running a full dynamical trajectory to measure
# correlations, so we here restrict to 5 samples.

for _ in 1:5
    for _ in 1:100
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end

# Define points that define a piecewise-linear path through reciprocal space,
# and a sampling density.

points = [
    [3 / 4, 3 / 4, 0],
    [0, 0, 0],
    [0, 1 / 2, 1 / 2],
    [1 / 2, 1, 0],
    [0, 1, 0],
    [1 / 4, 1, 1 / 4],
    [0, 1, 0],
    [0, -4, 0],
]
density = 50 # (Å)
path, xticks = reciprocal_space_path(cryst, points, density);

# Calculate ``I(𝐪, ω)`` intensities along this path with Lorentzian broadening
# on the scale of 0.1 meV.

formula = intensity_formula(sc, :perp; formfactors, kT=langevin.kT)
η = 0.1
iqw = intensities_interpolated(sc, path, formula)
iqwc = broaden_energy(sc, iqw, (ω, ω₀) -> lorentzian(ω - ω₀, η));

# Plot the intensity data on a clipped color scale

ωs = available_energies(sc)
heatmap(
    1:size(iqwc, 1),
    ωs,
    iqwc;
    colorrange=(0, maximum(iqwc) / 50),
    axis=(;
        xlabel="Momentum Transfer (r.l.u)",
        ylabel="Energy Transfer (meV)",
        xticks,
        xticklabelrotation=π / 5,
        aspect=1.4,
    ),
)

# ### Powder averaged intensity

# Define spherical shells in reciprocal space via their radii, in absolute units
# of 1/Å. For each shell, calculate and average the intensities at 100
# ``𝐪``-points, sampled approximately uniformly.

radii = 0:0.05:3.5 # (1/Å)
output = zeros(Float64, length(radii), length(ωs))
for (i, radius) in enumerate(radii)
    pts = reciprocal_space_shell(sc.crystal, radius, 100)
    is = intensities_interpolated(sc, pts, formula)
    is = broaden_energy(sc, is, (ω, ω₀) -> lorentzian(ω - ω₀, η))
    output[i, :] = mean(is; dims=1)[1, :]
end

# Plot resulting powder-averaged structure factor

heatmap(
    radii,
    ωs,
    output;
    axis=(xlabel="|Q| (Å⁻¹)", ylabel="Energy Transfer (meV)", aspect=1.4),
    colorrange=(0, 20.0),
)
