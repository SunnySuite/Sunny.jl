using Sunny, LinearAlgebra, Printf, Statistics
@assert pkgversion(Sunny) >= v"0.10.0"
using GLMakie
load_fast_blas()

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
(s, s_str) = (1/2, "1/2")
# (s, s_str) = (3/2, "3/2")
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))

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

plot_intensities(res; colormap=:jet, colorrange=(0, s+3/2),
                 title="s = $s_str", axis=(; xlabel="", ylabel="ω / J"))
