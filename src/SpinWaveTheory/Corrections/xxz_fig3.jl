# Reproduce Fig. 3a of Maksimov, Zhitomirsky and Chernyshev, PRB 94, 140407(R)
# (2016) [arXiv:1607.08238]. Calculates the dynamical structure factor of the
# easy-plane XXZ triangular-lattice antiferromagnet in an out-of-plane field,
# including ``1/s`` corrections.
#
# The reference calculation omitted interference between transverse and
# longitudinal channels. Sunny includes this additional ``1/s`` correction term,
# which is found to redistribute the two-magnon continuum.

using Sunny, LinearAlgebra, Printf
using GLMakie

# See the comment in fig4b.jl: Julia's OpenBLAS anti-scales when `bogoliubov!` is
# called from many threads, and Apple's Accelerate does not.
BLAS.set_num_threads(1)
Sys.isapple() && @eval using AppleAccelerate

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])

Δ = 0.9
hfrac = 0.2

sys = System(cryst, [1 => Moment(; s=1/2, g=1)], :dipole)
set_exchange!(sys, diagm([1.0, 1.0, Δ]), Bond(1, 1, [1, 0, 0]))
Hs = 3 * (Δ + 1/2)
set_field!(sys, [0, 0, -hfrac * Hs])

# Minimize to the umbrella state in the three-site magnetic cell.

sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
randomize_spins!(sys)
minimize_energy!(sys)

# The ground state breaks chiral symmetry in one of two ways. Reflect in the
# ``y`` plane as necessary to ensure a consistent chiral orientation.

if (sys.dipoles[1] × sys.dipoles[2])[3] < 0
    for site in eachsite(sys)
        set_dipole!(sys, diagm([+1, -1, 1]) * sys.dipoles[site], site)
    end
end
plot_spins(sys; ndims=2)

# The relevant path in Fourier space.

M = [1/2, 0, 0]
K′ = [1/3, 1/3, 0]
Γ = [0, 0, 0]
K = [2/3, -1/3, 0]
path = q_space_path(cryst, [M, K′, Γ, K, M, Γ], 400;
                    labels=["M", "K′", "Γ", "K", "M", "Γ"])

# Calculate intensities ``S^{αα}(𝐪, ω)`` to order 1/s. Enlarging the artificial
# broadening η by 4× accelerates the calculation by 16×.

η = 0.02 # vs. 0.005 used in Maksimov et al.
energies = 0:(η/2):2.0
swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
res = Sunny.corrected_intensities(swt, path; energies, η, threaded=true, verbose=true)

# Compare with Fig. 3a of Maksimov et al.

plot_intensities(res; colormap=:jet, colorrange=(0, 21.0),
                 title="XXZ triangular AFM, s = 1/2, Δ = $Δ, H = $hfrac Hs",
                 axis=(; xlabel="", ylabel="ω / J"))
