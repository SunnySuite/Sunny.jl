# # 8. Momentum transfer conventions
#
# Sunny defines the dynamical structure factor following conventions as in
# [Squire](https://doi.org/10.1017/CBO9781139107808) and
# [Boothroyd](https://groups.physics.ox.ac.uk/Boothroyd/PNS/),
#
# ```math
# \mathcal{S}^{αβ}(𝐪, ω) ≡ \frac{1}{2π} \int_{-∞}^{∞} e^{-iωt} ⟨\hat{M}^{†α}_𝐪(0) \hat{M}^β_𝐪(t)⟩ dt.
# ```
# 
# The magnetic moment in momentum space ``\hat{𝐌}_𝐪`` is obtained from the
# real-space density ``\hat{𝐌}(𝐫)`` using the Fourier transform convention,
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
# With appropriate contraction of spin components, ``\mathcal{S}^{αβ}(𝐪, ω)``
# can be directly related to the neutron scattering cross-section with ``𝐪``
# and ``ω`` denoting momentum and energy transfer **to** the sample. For models
# that lack inversion symmetry, the intensities at ``±𝐪`` may be inequivalent.
# We illustrate such a case using a 1D chain with competing Ising and
# Dzyaloshinskii–Moriya couplings between neighboring sites.
#
# Be aware that other codes, e.g. [SpinW](https://spinw.org/), may employ an
# alternate structure factor convention that effectively reverses the direction
# of momentum transfer, ``𝐪 → -𝐪``.

# ### 1D model lacking reflection symmetry

using Sunny, GLMakie

# The model will live on a 1D chain. Select the spacegroup P1 to effectively
# disable all symmetry analysis. This eliminates all symmetry-imposed
# constraints on the couplings. Furthermore, because all bonds become
# symmetry-inequivalent, each coupling must be specified independently (there is
# no automatic propagation to "equivalent" bonds).

latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], "P1")

# Include DM and Ising couplings between nearest neighbors on the chain. Letting
# ``j`` denote the site position along axis ``𝐚_3``, the full Hamiltonian is
# ```math
# ℋ = 𝐃 ⋅ ∑_j 𝐒_j × 𝐒_{j+1} - J ∑_j ∑_j S^z_j S^z_{j+1},
# ```
# with DM vector ``𝐃 = D ẑ``.

s = 3/2
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
J = 1.0
D = 0.2
z = [0, 0, 1]
set_exchange!(sys, dmvec(D * z) - J * z * z', Bond(1, 1, [0, 0, 1]))

# The relatively large Ising coupling favors one of two polarized states,
# ``±ẑ``. Use [`polarize_spins!`](@ref) to align ``𝐒`` with ``+ẑ``. This
# polarization could also be realized via an applied ``𝐁`` field in the
# _opposite_ direction, as documented in [`set_field!`](@ref).

polarize_spins!(sys, [0, 0, 1])
@assert energy(sys) ≈ - s^2
plot_spins(sys)

# ### Intensities from linear spin wave theory

# The [`SpinWaveTheory`](@ref) calculation shows a single band with dispersion
# ``ϵ(𝐪) = 2 s [J ± D \sin(2πq_3)]`` for the polarization state ``𝐒 = ± s
# ẑ``. There is a clear dependence on the sign of ``q_3``.

path = q_space_path(cryst, [[0, 0, -1/2], [0, 0, 0], [0, 0, +1/2]], 400)
swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
res = intensities_bands(swt, path)
plot_intensities(res; ylims=(0, 5))

# ### Intensities from classical spin dynamics

# Classical dynamics at low temperatures produces, in principle, the same
# excitation spectrum. Here, finite-size effects will limit ``𝐪``-space
# resolution. Use [`resize_supercell`](@ref) to study a chain of 32 sites.

sys2 = resize_supercell(sys, (1, 1, 32))

# Use [`Langevin`](@ref) dynamics to sample from the classical Boltzmann
# distribution at a relatively low temperature, ``k_B T / J = 0.03``

dt = 0.03 / J
kT = 0.03 * J
damping = 0.1
langevin = Langevin(dt; kT, damping)
suggest_timestep(sys2, langevin; tol=1e-2)
for _ in 1:10_000
    step!(sys, langevin)
end

# Use [`SampledCorrelations`](@ref) to collect statistics from classical spin
# dynamic trajectories.

sc = SampledCorrelations(sys2; dt, energies=range(0, 5, 100), measure=ssf_trace(sys2))
add_sample!(sc, sys2)
nsamples = 100
for _ in 1:nsamples
    for _ in 1:1000
        step!(sys2, langevin)
    end
    add_sample!(sc, sys2)
end

# In the limit ``T → 0``, the [`intensities`](@ref) match linear spin wave
# theory up to finite-size effects and statistical error. Additionally, the
# classical dynamics calculation shows the elastic peak at ``𝐪 = [0,0,0]``,
# associated with the ferromagnetic ground state.

res2 = intensities(sc, path; energies=:available, kT)
plot_intensities(res2)
