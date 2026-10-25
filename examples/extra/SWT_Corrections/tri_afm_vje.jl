# Self-consistent internal lines in the scheme of Veillette, James and Essler,
# PRB 72, 134429 (2005) (https://arxiv.org/abs/cond-mat/0506667), here for the
# triangular-lattice Heisenberg antiferromagnet. The magnon energies in the loop
# denominators are replaced by the renormalized ones, keeping the vertices and
# Bogoliubov coefficients harmonic, and the replacement is iterated. The bare
# propagator of the Dyson equation stays that of LSWT, as in their Eq. (40).
#
# The scheme is uncontrolled. At s = 4 it converges in a few iterations, but at
# s ≲ 5/2 the iteration drifts away from a fixed point.

using Sunny, LinearAlgebra
using GLMakie
load_blas_for_threading()

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
s = 4
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
randomize_spins!(sys)
minimize_energy!(sys)
swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
η = 0.03s

# Each iteration dresses the internal lines with the on-shell poles computed
# from the previous ones. The poles are found on an 8×8 grid of the magnetic
# zone and interpolated between.

vacuum = Sunny.MagnonVacuum(swt)
qM = [[1/2, 0, 0]]
for iter in 1:4
    global vacuum = Sunny.dressed_vacuum(Sunny.OneLoop(swt; η, vacuum); grid=(8, 8, 1), threaded=true)
    bands = Sunny.corrected_intensities_bands(swt, qM; η, vacuum)
    println("Iteration $iter: bands at M = ", round.(sort(bands.disp[:]); digits=4))
end

# Compare intensities with harmonic and self-consistent internal lines. The
# latter moves the threshold of the two-magnon continuum with the magnons.

qpts = [[2/3, -1/3, 0], [0, 0, 0], [1/2, 0, 0], [1/6, 1/6, 0], [0, 1/4, 0]]
path = q_space_path(cryst, qpts, 200; labels=["K", "Γ", "M", "Y₁", "Y"])
energies = 0:(η/2):2.5s
res0 = Sunny.corrected_intensities(swt, path; energies, η, threaded=true)
res1 = Sunny.corrected_intensities(swt, path; energies, η, threaded=true, vacuum)

fig = Figure(size=(1000, 400))
plot_intensities!(fig[1, 1], res0; colorrange=(0, 2s), title="Harmonic internal lines")
plot_intensities!(fig[1, 2], res1; colorrange=(0, 2s), title="Self-consistent internal lines")
fig
