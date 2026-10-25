# Reproduce Fig. 3(a) of Maksimov, Zhitomirsky and Chernyshev, PRB 94, 140407(R)
# (2016) [arXiv:1607.08238]. Calculates the dynamical structure factor of the
# easy-plane XXZ triangular-lattice antiferromagnet in an out-of-plane field,
# including 1/s corrections.
#
# Sunny makes one small correction to the reference calculation: it retains
# interference between transverse and longitudinal channels. For this model, the
# correction removes a small fraction of the two-magnon weight, which itself is
# small compared to the quasiparticle branch intensities.

using Sunny, LinearAlgebra, Printf
using GLMakie

# See the comment in fig4b.jl: Julia's OpenBLAS anti-scales when `bogoliubov!` is
# called from many threads, and Apple's Accelerate does not.
BLAS.set_num_threads(1)
Sys.isapple() && @eval using AppleAccelerate

# The paper's labeled wavevectors, in r.l.u. of the one-site hexagonal cell. K is
# the ordering wavevector (4π/3, 0); K′ is the adjacent corner (2π/3, 2π/√3);
# M = (π, π/√3) is the midpoint of the zone edge that joins them.
Kpt = [2/3, -1/3, 0]
Kppt = [1/3, 1/3, 0]
Mpt = [1/2, 0, 0]

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])

(s, Δ, hfrac) = (1/2, 0.9, 0.2)

sys = System(cryst, [1 => Moment(; s, g=1)], :dipole)
set_exchange!(sys, diagm([1.0, 1.0, Δ]), Bond(1, 1, [1, 0, 0]))
Hs = 6s * (Δ + 1/2)
set_field!(sys, [0, 0, -hfrac * Hs])

# Umbrella state in the three-site magnetic cell. The winding is taken along
# 𝐐 = K rather than -K, which selects the chirality domain in which K is the
# corner that the field pushes down, as in the paper.
sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
θ = asin(hfrac)
Q = cryst.recipvecs * Kpt
for site in eachsite(sys)
    φ = dot(Q, global_positions(sys)[site])
    set_dipole!(sys, [cos(θ)cos(φ), cos(θ)sin(φ), sin(θ)], site)
end

swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))

# Harmonic cross-check against their Fig. 1. Folded into the three-site cell, the
# bands at Γ are the extended-zone ε_𝐤 at 𝐤 = 0, ±K, i.e. the Goldstone mode
# together with the two corner energies, which the field splits.
ε0 = sort(vec(dispersion(swt, [[0, 0, 0]])))
@printf("harmonic: ε(Γ) = %.4f, ε(K) = %.4f, ε(K′) = %.4f  (their Fig. 1: 0, 0.27, 1.18)\n",
        ε0[1], ε0[2], ε0[3])
@printf("umbrella: θ = %.2f°, Hs = %.3f J, H = %.3f J\n", rad2deg(θ), Hs, hfrac*Hs)

path = q_space_path(cryst, [Mpt, Kppt, [0, 0, 0], Kpt, Mpt, [0, 0, 0]], 200;
                    labels=["M", "K′", "Γ", "K", "M", "Γ"])

η = 0.01
energies = 0:(η/2):2.0
@time res = Sunny.corrected_intensities(swt, path; energies, η, tol=0.01, threaded=true, verbose=true)

# Their colour scale is a rainbow from zero, cut off at ≈21 so that the Bragg
# divergence at the ordering wavevector is allowed to saturate.
plot_intensities(res; colormap=:jet, colorrange=(0, 21.0),
                 title="XXZ triangular AFM, s = 1/2, Δ = $Δ, H = $hfrac Hs",
                 axis=(; xlabel="", ylabel="ω / J"))
