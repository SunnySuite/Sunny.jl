# Reproduce Fig. 4 of Zhang et al. (https://arxiv.org/abs/2508.21142).
# Calculates the longitudinal structure factor Sᶻᶻ of the up-up-down (UUD)
# plateau phase of the triangular XXZ antiferromagnet in a field.
#
# The UUD state is collinear, so Sᶻᶻ is purely two-magnon. A one-loop
# correction (`dyson=:nambu`) gives only a faint continuum. The ladder
# resummation (`dyson=:ladder`) includes the interaction between the two
# magnons, which binds them into sharp branches. This reproduces the truncated
# Hilbert space calculation (THED) of the paper.

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
    vacuum = Sunny.MagnonVacuum(swt, Sunny.hartree_fock_correction(swt; tol=1e-4).terms2)
    for (col, dyson) in enumerate((:nambu, :ladder))
        res = Sunny.corrected_intensities(swt, path; energies, η, kernel, dyson, vacuum, threaded=true, verbose=true)
        plot_intensities!(fig[row, col], res; colormap=:viridis, colorrange,
                          title="Sᶻᶻ, Δ = $Δ, dyson = :$dyson")
    end
end
fig

# At Δ = 5, the ladder gives the four flat bound-state branches of Fig. 4(a),
# below the continuum. At Δ = 2, the upper two branches enter the continuum and
# broaden, as in Fig. 4(b). At Δ = 1, the UUD state is classically unstable, but
# is stabilized by quantum fluctuations. The Hartree-Fock mean field of
# `MagnonVacuum` captures this, gapping the soft mode. The lower branch then has
# minimum 0.69 J at the K point, compared with 0.75 J from THED. The weak mode
# near 1.5 J that MPS finds is attributed to four-magnon states, which both THED
# and the ladder omit.
