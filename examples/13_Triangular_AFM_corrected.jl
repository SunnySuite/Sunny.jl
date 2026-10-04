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
(s, s_str, ωmax) = (1/2, "1/2", 3.5)
## (s, s_str, ωmax) = (3/2, "3/2", 10)
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))

# The ground state is a 120° spiral in the three-site magnetic cell.

sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)

# Calculate spin-wave intensities with loop corrections at order ``1/s`` using
# `corrected_intensities` (an experimental feature). The option
# `dyson=:particle` applies the particle-sector Dyson equation as in Eq. (12) of
# Mourigal et al. Sunny's default, `dyson=:nambu`, would instead apply the Dyson
# equation in the full particle-hole space, which may be more accurate. The
# regulator ``η > 0`` defines an effective energy resolution. Numerical costs
# grow as ``η^{-D}`` in the effective model dimension ``D``.

swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
qpts = [[2/3, -1/3, 0], [0, 0, 0], [1/2, 0, 0], [1/6, 1/6, 0], [0, 1/4, 0]]
labels=["K", "Γ", "M", "Y₁", "Y"]
path = q_space_path(cryst, qpts, 200; labels)
η = 0.03s
energies = 0:(η/2):ωmax
res = Sunny.corrected_intensities(swt, path; energies, η, threaded=true, verbose=true, dyson=:particle)
;#hide

# The intensities can be compared with Fig. 4 of Mourigal et al. (corrected in
# the
# [erratum](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.099901)).

plot_intensities(res; colormap=:jet, colorrange=(0, s+3/2),
                 title="s = $s_str", axis=(; xlabel="", ylabel="ω / J"))
