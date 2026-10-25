# Reproduce Fig. 4a of Mourigal, Fuhrman, Chernyshev and Zhitomirsky, PRB 88,
# 094407 (2013) [arXiv:1306.1231]. Calculates the dynamical structure factor of
# the triangular-lattice Heisenberg antiferromagnet, including ``1/s``
# corrections.
#
# The reference calculation omitted interference between transverse and
# longitudinal channels. Sunny includes this additional ``1/s`` correction term,
# which is found to redistribute the two-magnon continuum.

using Sunny, LinearAlgebra, Printf, Statistics
using GLMakie

# `threaded=true` parallelizes over 𝐪, but the OpenBLAS that Julia ships does not merely
# fail to scale on matrices this small, it anti-scales: called from 18 threads,
# `bogoliubov!` slows by a factor of 30, so total throughput saturates at two cores'
# worth. The contention is inside the library and `set_num_threads` does not reach it.
# Apple's Accelerate is reentrant, and takes `bogoliubov!` from 0.6× to 12.6×, at the
# price of a slightly slower serial eigensolve.
BLAS.set_num_threads(1)
Sys.isapple() && @eval using AppleAccelerate

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])

(s, s_str) = (1/2, "1/2")
## (s, s_str) = (3/2, "3/2")

sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))

# 120° spiral on the triangular lattice, in the three-site magnetic cell
sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
randomize_spins!(sys)
minimize_energy!(sys)

swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))

qpts = [[2/3, -1/3, 0], [0, 0, 0], [1/2, 0, 0], [1/6, 1/6, 0], [0, 1/4, 0]]
labels=["K", "Γ", "M", "Y₁", "Y"]
path = q_space_path(cryst, qpts, 200; labels)

η = 0.03s
energies = 0:(η/2):(20s/3)
res = Sunny.corrected_intensities(swt, path; energies, η, threaded=true, verbose=true)

# Compare with Fig. 4a of Mourigal et al.

plot_intensities(res; colormap=:jet, colorrange=(0, s+3/2),
                 title="s = $s_str", axis=(; xlabel="", ylabel="ω / J"))
