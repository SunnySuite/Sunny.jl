# Reproduce Fig. 4 of Zhang et al. (https://arxiv.org/abs/2508.21142) using
# diagrammatic ladder resummation. Calculates the longitudinal structure factor
# Sᶻᶻ of the up-up-down (UUD) plateau phase of the triangular XXZ
# antiferromagnet in a field.
#
# The UUD state is collinear, so Sᶻᶻ is purely two-magnon. The `:one_loop`
# correction scheme gives only a faint continuum. The `:two_magnon_ladder`
# scheme allows magnon pairs to interact and bind into sharp branches.

using Sunny, LinearAlgebra
using GLMakie
load_blas_for_threading()

# XXZ exchange with anisotropy Δ, and the field 1.905 J of the paper.

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole)
set_exchange!(sys, diagm([1, 1, 0]), Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, diagm([0, 0, 1]), Bond(1, 1, [1, 0, 0]), :Δ => 5.0)
set_field!(sys, [0, 0, 1.905])

# The UUD state in the three-site magnetic cell. It is the classical ground
# state at Δ = 5, and remains a stationary point as Δ is lowered below.

sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
randomize_spins!(sys)
minimize_energy!(sys)

# The energy scale η is a regulator. Additional Gaussian broadening σ is applied
# for consistency with the published figure. Color ranges tuned to match the
# paper's arbitrary units.

η = 0.05
kernel = gaussian(; σ=0.1)
path = q_space_path(cryst, [[0, 0, 0], [1/3, 1/3, 0], [1/2, 1/2, 0]], 120; labels=["Γ", "K", "M"])
panels = [
    (; Δ=5, energies=range(8, 12, 321), colorrange=(0, 0.022)),
    (; Δ=2, energies=range(2, 6, 321), colorrange=(0, 0.14)),
    (; Δ=1, energies=range(0, 4, 321), colorrange=(0, 0.19)),
]

fig = Figure(size=(1000, 1100))
for (row, (; Δ, energies, colorrange)) in enumerate(panels)
    set_param!(sys, :Δ, Δ)
    measure = ssf_custom((q, ssf) -> real(ssf[3, 3]), sys; apply_g=false)
    swt = SpinWaveTheory(sys; measure)
    swt_hf = Sunny.to_hartree_fock_corrected_vacuum(swt)
    for (col, scheme) in enumerate((:one_loop, :two_magnon_ladder))
        res = Sunny.corrected_intensities(swt_hf, path, scheme; energies, η,
                                          kernel, threaded=true, verbose=true)
        plot_intensities!(fig[row, col], res; colormap=:viridis, colorrange,
                          title="Sᶻᶻ, Δ = $Δ, scheme = :$scheme")
    end
end
fig

# Zhang et al. used truncated Hilbert space exact diagonalization (THED), which
# builds the bosonic Hamiltonian for states with at most two magnons, and
# diagonalizes it exactly. Sunny's `:two_magnon_ladder` scheme is a diagrammatic
# resummation that captures the same two-magnon bound states. It additionally
# includes processes that create or annihilate magnon pairs out of the ground
# state, which reflect its quantum fluctuations.
#
# At Δ = 5, the ladder scheme gives the four flat bound-state branches below the
# continuum. At Δ = 2, the upper two branches enter the continuum and broaden.
# At Δ = 1, the UUD state is classically unstable, but is stabilized by quantum
# fluctuations. The calculation captures this by starting from the Hartree-Fock
# corrected ground state `swt_hf`, which gives the soft mode a gap. The lower
# branch then has its K point minimum at 0.69 J (slighly lower than 0.75 J from
# THED). The paper's matrix product state calculation also finds a weak mode
# near 1.5 J, attributed to four-magnon states, which both THED and the ladder
# scheme omit.
