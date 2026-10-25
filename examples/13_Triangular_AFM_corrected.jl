# # 13. Perturbative corrections to spin wave theory
#
# This tutorial calculates the spin wave spectrum of a triangular
# antiferromagnet, including all perturbative corrections at order ``1/s``. The
# study follows [Mourigal, Fuhrman, Chernyshev and Zhitomirsky, Phys. Rev. B
# **88**, 094407 (2013)](https://doi.org/10.1103/PhysRevB.88.094407).

using Sunny, LinearAlgebra, Printf, Statistics
using GLMakie
load_blas_for_threading()

# Triangular lattice with antiferromagnetic nearest-neighbor interactions.

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
(s, s_str) = (1/2, "1/2")
## (s, s_str) = (3/2, "3/2")
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))

# The ground state is a 120° spiral in the three-site magnetic cell.

sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)

# Calculate spin-wave intensities, including all corrections to order ``1/s``.
# The regulator ``η > 0`` defines an energy resolution. The numerical cost of
# the momentum integrals scales like ``η^{-D}`` in effective dimension ``D``.
# The function `corrected_intensities` is currently experimental and subject to
# change.

swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
qpts = [[2/3, -1/3, 0], [0, 0, 0], [1/2, 0, 0], [1/6, 1/6, 0], [0, 1/4, 0]]
labels=["K", "Γ", "M", "Y₁", "Y"]
path = q_space_path(cryst, qpts, 200; labels)
η = 0.03s
energies = 0:(η/2):(20s/3)
res = Sunny.corrected_intensities(swt, path; energies, η, threaded=true, verbose=true)
;#hide

# Plot the corrected intensities. Gray pixels indicate regions where the theory
# is uncontrolled, as detected by resummed magnon poles that fall below half of
# the lowest harmonic energy.
#
# A visual comparison with Fig. 4 of Mourigal et al. highlights significant
# deviations, especially in the ``s = 1/2`` case at lower energies. The two
# schemes agree at order ``1/s``, but differ in their resummation procedure.
# Sunny solves the Dyson equation in the full Nambu (particle/hole) space,
# whereas Mourigal et al. project to the particle space alone [concretely, see
# Eq. (12) for their resummation procedure]. Sunny's choice is the standard form
# of the Dyson equation for bosons with broken symmetry, and is the one that
# respects Goldstone's theorem. It can also produce renormalized poles with
# imaginary frequency and these signal where perturbation theory becomes
# uncontrolled. To more closely reproduce Fig. 4 of Mourigal et al., recalculate
# `corrected_intensities` with the hidden option `resummation=:particle`.

plot_intensities(res; colormap=:jet, colorrange=(0, s+3/2),
                 title="s = $s_str", axis=(; xlabel="", ylabel="ω / J"))
