# # 8. Momentum transfer conventions
#
# Sunny defines the dynamical spin structure factor following conventions such
# as in [Squire](https://doi.org/10.1017/CBO9781139107808) and
# [Boothroyd](https://groups.physics.ox.ac.uk/Boothroyd/PNS/).
#
# ```math
# \mathcal{S}^{α, β}(𝐪, ω) ≡ \frac{1}{2π} \int_{-∞}^{∞} e^{-iωt} ⟨\hat{M}^{†α}_𝐪(0) \hat{M}^β_𝐪(t)⟩ dt.
# ```
# 
# Momentum space operators ``\hat{𝐌}_𝐪`` are obtained from the real-space
# density ``\hat{𝐌}(𝐫)`` using the Fourier transform convention,
#
# ```math
# \hat{𝐌}_𝐪 ≡ \int_V e^{+ i 𝐪⋅𝐫} \hat{𝐌}(𝐫) d𝐫.
# ```
#
# The structure factor, integrated over a finite ``𝐪``-region, is extensive in
# sample volume ``V``. Sunny will report it as an intensive quantity by dividing
# by the number of chemical cells in the sample. For full details, see the
# documentation page [Structure Factor Conventions](@ref).
#
# With appropriate contraction of spin components, ``\mathcal{S}^{α, β}(𝐪, ω)``
# can be directly related to the neutron scattering cross-section, where ``𝐪``
# and ``ω`` denote momentum and energy transfer _to_ the sample. For models that
# lack inversion symmetry, the intensities at ``± 𝐪`` may become inequivalent.
# We illustrate such a case using a 1D chain with competing Ising and
# Dzyaloshinskii–Moriya couplings between neighboring sites.

using Sunny, GLMakie

# Selecting the P1 spacegroup will effectively disable all symmetry analysis.
# This can be a convenient way to avoid symmetry-imposed constraints on the
# couplings. A disadvantage is that all bonds are treated as
# symmetry-inequivalent, such that each coupling within the unit cell must be
# specified independently.

latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], "P1")

# Consider a 1D chain oriented along the third lattice vector, such that site
# ``j`` has position ``𝐫_j = j 𝐚_3``. The Hamiltonian includes DM and
# Ising-like couplings between nearest neighbors on the chain,
#
# ```math
# ℋ = 𝐃 ⋅ ∑_j 𝐒_j × 𝐒_{j+1} - J ∑_j ∑_j S^z_j S^z_{j+1}.
# ```
# Select the DM vector ``𝐃 = D ẑ``.

s = 3/2
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
J = 1.0
D = 0.2
z = [0, 0, 1]
set_exchange!(sys, dmvec(D * z) - J * z * z', Bond(1, 1, [0, 0, 1]))

# The relatively large Ising coupling favors one of two polarized states,
# ``±ẑ``. Align spins ``𝐒`` in the ``+ẑ`` direction to break the symmetry by
# hand.

polarize_spins!(sys, [0, 0, 1])
@assert energy(sys) ≈ - s^2
plot_spins(sys)

# ### Calculation using linear spin wave theory

# The [`SpinWaveTheory`](@ref) calculation shows a single band with dispersion
# ``ϵ(𝐪) = 2 s [J ± D \sin(2πq_3)]`` for the polarization state ``𝐒 = ± s ẑ``.
# Note the sensitivity to the sign of ``q_3``.

path = q_space_path(cryst, [[0, 0, -1/2], [0, 0, +1/2]], 400)
swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
res = intensities_bands(swt, path)
plot_intensities(res; ylims=(0, 5))

# ### Calculation using classical spin dynamics

# Enlarge the system with [`resize_supercell`](@ref)

sys2 = resize_supercell(sys, (1, 1, 32))

# Use [`Langevin`](@ref) dynamics to sample from the classical Boltzmann
# distribution at a low temperature.

dt = 0.03 / abs(J)
kT = 0.03 * abs(J)
damping = 0.1
langevin = Langevin(dt; kT, damping)
suggest_timestep(sys2, langevin; tol=1e-2)
for _ in 1:10_000
    step!(sys, langevin)
end

# Estimate intensities with [`SampledCorrelations`](@ref), which employs
# classical spin dynamics.

sc = SampledCorrelations(sys2; dt, energies=range(0, 5, 100), measure=ssf_trace(sys2))
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
