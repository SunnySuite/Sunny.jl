# Reproduce Fig. 2 of Mourigal, Fuhrman, Chernyshev and Zhitomirsky, PRB 88,
# 094407 (2013) (https://arxiv.org/abs/1306.1231). Calculates the magnon
# spectral function A₁₁(𝐪, ω) of the triangular-lattice Heisenberg
# antiferromagnet, including ``1/s`` corrections.
#
# Unlike the structure factor, A₁₁ is not an observable. It is the particle
# block of the magnon Green's function, available from the lower-level function
# `corrected_channels` with `spectral=true`.

using Sunny, LinearAlgebra
using GLMakie
load_blas_for_threading()

# Triangular lattice with antiferromagnetic nearest-neighbor interactions.

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
(s, s_str, ωmax) = (1/2, "1/2", 3.5)
# (s, s_str, ωmax) = (3/2, "3/2", 10)
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))

# The ground state is a 120° spiral in the three-site magnetic cell.

sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
randomize_spins!(sys)
minimize_energy!(sys)

# The path K-Γ-M-Y₁-Y of Mourigal et al.

qpts = [[2/3, -1/3, 0], [0, 0, 0], [1/2, 0, 0], [1/6, 1/6, 0], [0, 1/4, 0]]
labels = ["K", "Γ", "M", "Y₁", "Y"]
path = q_space_path(cryst, qpts, 200; labels)

# Calculate the spectral matrix at order ``1/s``, applying the particle-sector
# Dyson equation as in Eq. (12) of Mourigal et al.

η = 0.03s
energies = 0:(η/2):ωmax
swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
res = Sunny.corrected_channels(swt, path; energies, η, dyson=:particle, spectral=true,
                               threaded=true, verbose=true)

# The three-site magnetic cell folds 𝐪 together with 𝐪 ± 𝐐, where 𝐐 is the
# ordering wavevector K. The band belonging to 𝐪 itself is the one whose
# harmonic energy matches Eq. (11) of Mourigal et al., for the unfolded one-site
# cell.

γ(q) = (cos(2π*q[1]) + cos(2π*q[2]) + cos(2π*(q[1] + q[2]))) / 3
ε(q) = 3s * sqrt((1 - γ(q)) * (1 + 2γ(q)))
bands = [argmin(abs.(res.disp[:, iq] .- ε(q))) for (iq, q) in enumerate(path.qs)]
A11 = [real(res.specfunc[n, n, iω, iq]) for iω in eachindex(energies), (iq, n) in enumerate(bands)]

# Compare with Fig. 2 of Mourigal et al. The sharp lower edges along K-Γ and
# Y₁-Y are thresholds of the two-magnon continuum. Peak heights scale as
# ``1/η ∝ 1/s``, and so does the color range.

fig = Figure(size=(600, 400))
plot_intensities!(fig[1, 1], Sunny.Intensities(cryst, res.qpts, res.energies, A11);
                  colormap=:jet, colorrange=(0, 3/2s), title="A₁₁(𝐪, ω), s = $s_str",
                  axis=(; xlabel="", ylabel="ω / J"))
lines!(eachindex(path.qs), ε.(path.qs); color=:white, linestyle=:dash)
fig
