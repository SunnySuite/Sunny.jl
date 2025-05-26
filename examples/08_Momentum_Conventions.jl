# # 8. Momentum transfer conventions
#
# This example illustrates Sunny's conventions for dynamical structure factor
# intensities, ``\mathcal{S}(𝐪,ω)``, as documented in the page [Structure
# Factor Conventions](@ref). The variables ``𝐪`` and ``ω`` describe momentum
# and energy transfer _to_ the sample.
#
# The structure factor intensities at ``± 𝐪`` may be inequivalent if the model
# lacks inversion symmetry. This tutorial considers a 1D chain with
# Dzyaloshinskii–Moriya coupling between neighboring sites. An external field
# then breaks time-reversal symmetry, giving rise to an inequivalence of
# intensities ``\mathcal{S}(±𝐪,ω)``

using Sunny, GLMakie

# Selecting the P1 spacegroup will effectively disable all symmetry analysis.
# This can be a convenient way to avoid symmetry-imposed constraints on the
# couplings. A disadvantage is that all bonds are treated as
# symmetry-inequivalent, such that each coupling within the unit cell must be
# specified independently.

latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], "P1")

# Construct a 1D chain system that extends along the global Cartesian ``ẑ``
# axis. The Hamiltonian includes DM and Zeeman coupling terms, ``ℋ = ∑_j D ẑ ⋅
# (𝐒_j × 𝐒_{j+1}) - ∑_j 𝐁 ⋅ μ_j``, where ``μ_j = - 2 𝐒_j`` is the
# [`magnetic_moment`](@ref) and ``𝐁 ∝ ẑ``.

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
B = 0.8
D = 0.2
set_exchange!(sys, dmvec([0, 0, D]), Bond(1, 1, [0, 0, 1]))
set_field!(sys, [0, 0, B])

# The large external field fully polarizes the system. Here, the DM coupling
# contributes nothing, leaving only Zeeman coupling.

randomize_spins!(sys)
minimize_energy!(sys)
@assert magnetic_moment(sys, (1, 1, 1, 1)) ≈ [0, 0, 2 * sign(B)]
@assert energy(sys) ≈ - 2 * abs(B)

# ### Calculation using linear spin wave theory

# The [`SpinWaveTheory`](@ref) calculation shows a single band with dispersion
# ``ϵ(𝐪) = 2 |B| (1 - \sin(2πq_3) D / B)`` and uniform intensity. Notice the
# different excitation energies at ``𝐪`` and ``-𝐪``.

path = q_space_path(cryst, [[0,0,-1/2], [0,0,+1/2]], 400)
swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
res = intensities_bands(swt, path)
plot_intensities(res; ylims=(0, 3))

# ### Calculation using classical spin dynamics

# Enlarge the system with [`resize_supercell`](@ref)

sys2 = resize_supercell(sys, (1, 1, 32))

# Use [`Langevin`](@ref) dynamics to sample from the classical Boltzmann
# distribution at a low temperature.

dt = 0.05 / abs(B)
kT = 0.05 * abs(B)
damping = 0.1
langevin = Langevin(dt; kT, damping)
suggest_timestep(sys2, langevin; tol=1e-2)
for _ in 1:10_000
    step!(sys, langevin)
end

# Estimate intensities with [`SampledCorrelations`](@ref), which employs
# classical spin dynamics.

sc = SampledCorrelations(sys2; dt, energies=range(0, 3, 100), measure=ssf_trace(sys2))
add_sample!(sc, sys2)
nsamples = 100
for _ in 1:nsamples
    for _ in 1:1000
        step!(sys2, langevin)
    end
    add_sample!(sc, sys2)
end

# In the limit ``T → 0``, the excitation band matches linear spin wave theory up
# to finite-size effects and statistical error. Additionally, there is an
# elastic peak at ``𝐪 = [0,0,0]`` associated with the ferromagnetic ground
# state.

res2 = intensities(sc, path; energies=:available, kT)
plot_intensities(res2)
