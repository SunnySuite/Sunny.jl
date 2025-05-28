# # 8. Momentum transfer conventions
#
# Sunny defines the dynamical spin structure factor following the sign
# conventions in [Squire](https://doi.org/10.1017/CBO9781139107808) and
# [Boothroyd](https://groups.physics.ox.ac.uk/Boothroyd/PNS/) textbooks,
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
# and ``ω`` denote momentum and energy transfer _to_ the sample. Special care is
# required for models that lack inversion symmetry, in which case ``± 𝐪``
# become inequivalent. We illustrate such a case using a 1D chain with
# Dzyaloshinskii–Moriya coupling between neighboring sites. With an external
# field, the inequivalence of ``\mathcal{S}(±𝐪,ω)`` becomes apparent.

using Sunny, GLMakie

# Selecting the P1 spacegroup will effectively disable all symmetry analysis.
# This can be a convenient way to avoid symmetry-imposed constraints on the
# couplings. A disadvantage is that all bonds are treated as
# symmetry-inequivalent, such that each coupling within the unit cell must be
# specified independently.

latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], "P1")

# Consider a 1D chain oriented along the third lattice vector, such that site
# ``j`` has position ``𝐫_j = j 𝐚_3``. The Hamiltonian includes DM coupling
# between nearest neighbors and Zeeman coupling,
#
# ```math
# ℋ = 𝐃 ⋅ ∑_j 𝐒_j × 𝐒_{j+1} - 𝐁 ⋅ ∑_j μ_j.
# ```
# Select DM vector ``𝐃 = D ê`` and external field ``𝐁 = B ê`` in an
# arbitrary direction ``ê``. This field couples to each dimensionless
# [`magnetic_moment`](@ref) ``μ_j = - g 𝐒_j``. The convention of Sunny is to
# absorb the ``μ_B`` factor into ``𝐁``, giving it units of energy. Fix ``g =
# -2`` for simplicity; other values could be viewed as a rescaling of ``𝐁``.

s = 2
sys = System(cryst, [1 => Moment(; s, g=-2)], :dipole)
B = 1.0
D = 0.1
e = [1, 0, 0]
set_exchange!(sys, dmvec(D * e), Bond(1, 1, [0, 0, 1]))
set_field!(sys, B * e)

# The large external field fully polarizes the system. Here, the DM coupling
# contributes nothing, leaving only Zeeman coupling.

randomize_spins!(sys)
minimize_energy!(sys)
@assert magnetic_moment(sys, (1, 1, 1, 1)) ≈ 2s * sign(B) * e
@assert energy(sys) ≈ - 2s * abs(B)

# ### Calculation using linear spin wave theory

# The [`SpinWaveTheory`](@ref) calculation shows a single band with dispersion
# ``ϵ(𝐪) = |B| + sgn(B) 2 s D \sin(2πq_3)`` and uniform intensity. Notice
# the different excitation energies at ``𝐪`` and ``-𝐪``.

path = q_space_path(cryst, [[0, 0, -1/2], [0, 0, +1/2]], 400)
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
