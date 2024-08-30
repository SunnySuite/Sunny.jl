# # 8. Momentum transfer conventions
#
# This example illustrates Sunny's conventions for dynamical structure factor
# intensities, ``\mathcal{S}(𝐪,ω)``, as documented in the page [Structure
# Factor Conventions](@ref). The variables ``𝐪`` and ``ω`` describe momentum
# and energy transfer _to_ the sample.
#
# For systems without inversion-symmetry, the structure factor intensities at
# ``± 𝐪`` may be inequivalent. To highlight this, consider a simple 1D chain
# that includes only Dzyaloshinskii–Moriya interactions between neighboring
# sites. Coupling to an external field then breaks time-reversal symmetry,
# giving rise to an inequivalence of intensities ``\mathcal{S}(±𝐪,ω)``

using Sunny, GLMakie

# Selecting the P1 spacegroup will effectively disable all symmetry analysis.
# This can be a convenient way to avoid symmetry-imposed constraints on the
# couplings. A disadvantage is that all bonds are treated as inequivalent, and
# Sunny will therefore not be able to propagate any couplings between bonds.

latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], "P1")

# Construct a 1D chain system that extends along the global Cartesian ``ẑ``
# axis. The Hamiltonian includes DM and Zeeman coupling terms, ``ℋ = ∑_j D ẑ ⋅
# (𝐒_j × 𝐒_{j+1}) - ∑_j 𝐁 ⋅ μ_j``, where ``μ_j = - g 𝐒_j`` is the
# [`magnetic_moment`](@ref) and ``𝐁 ∝ ẑ``.

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(1, 1, 25))
D = 0.1
B = 5D
set_exchange!(sys, dmvec([0, 0, D]), Bond(1, 1, [0, 0, 1]))
set_field!(sys, [0, 0, B])

# The large external field fully polarizes the system. Here, the DM coupling
# contributes nothing, leaving only Zeeman coupling.

randomize_spins!(sys)
minimize_energy!(sys)
@assert energy_per_site(sys) ≈ -10D

# Sample from the classical Boltzmann distribution at a low temperature.

dt = 0.1
kT = 0.02
damping = 0.1
langevin = Langevin(dt; kT, damping)
suggest_timestep(sys, langevin; tol=1e-2)
for _ in 1:10_000
    step!(sys, langevin)
end

# The Zeeman coupling polarizes the magnetic moments in the ``𝐁 ∝ ẑ``
# direction. The spin dipoles, however, are anti-aligned with the magnetic
# moments, and therefore point towards ``-ẑ``. This is shown below.

plot_spins(sys)

# Estimate the dynamical structure factor using classical dynamics.

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

# Calculate the same quantity with linear spin wave theory at ``T = 0``. Because
# the ground state is fully polarized, the required magnetic cell is a ``1×1×1``
# grid of chemical unit cells.

sys_small = resize_supercell(sys, (1,1,1))
minimize_energy!(sys_small)
swt = SpinWaveTheory(sys_small; measure=ssf_trace(sys_small))
res2 = intensities_bands(swt, path)

# This model system has a single magnon band with dispersion ``ϵ(𝐪) = 1 - D/B
# \sin(2πq₃)`` and uniform intensity. Both calculation methods reproduce this
# analytical solution. Observe that ``𝐪`` and ``-𝐪`` are inequivalent. The
# structure factor calculated from classical dynamics additionally shows an
# elastic peak at ``𝐪 = [0,0,0]``, reflecting the ferromagnetic ground state.

fig = Figure(size=(768, 300))
plot_intensities!(fig[1, 1], res1; title="Classical dynamics")
plot_intensities!(fig[1, 2], res2; title="Spin wave theory")
fig
