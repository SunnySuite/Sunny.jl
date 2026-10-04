# Two-magnon bound states of the up-up-down (UUD) phase of the triangular XXZ
# antiferromagnet, after Fig. 4 of Zhang et al., arXiv:2508.21142. That work
# diagonalizes the Hamiltonian in the truncated space of one- and two-magnon
# states (THED). Here `dyson=:ladder` solves the same two-magnon problem as the
# resolvent of a quadratic model of magnons and interacting pairs, in the full
# Nambu space.
#
# The UUD state is collinear, so the cubic vertex vanishes and the longitudinal
# structure factor S^zz is purely two-magnon. Without the interaction between
# the two magnons, `:nambu` gives only a faint continuum. The ladder binds the
# pairs into sharp branches that carry most of the weight. Their energies agree
# with THED at Δ = 5, where the paper finds good agreement with MPS.

using Sunny, LinearAlgebra
using GLMakie
load_blas_for_threading()

# XXZ exchange with a field along z, in units of J. The field g μB B/J = 1.905
# of the paper lies inside the classical plateau for both anisotropies below.
# At Δ = 1 the UUD state is classically degenerate with canted states, so that
# LSWT about it is unstable without the order-by-disorder that THED adds by
# hand, and that case is omitted.

function uud_system(Δ)
    cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole)
    set_exchange!(sys, diagm([1, 1, Δ]), Bond(1, 1, [1, 0, 0]))
    set_field!(sys, [0, 0, 1.905])
    sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
    for (i, sz) in enumerate((-1, -1, +1))
        set_dipole!(sys, [0, 0, sz], (1, 1, 1, i))
    end
    return sys
end

# The magnon lines are dressed by the Hartree-Fock mean field, which is the
# leading 1/s correction that THED includes in its one-magnon energies. With
# the counterterm that `MagnonVacuum` subtracts, the result differs from a 1/s
# expansion about LSWT only at higher order.

path = q_space_path(Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]]),
                    [[0, 0, 0], [1/3, 1/3, 0], [1/2, 1/2, 0]], 120)
η = 0.05

fig = Figure(size=(1000, 750))
for (row, (Δ, energies)) in enumerate(((5, range(7, 12, 401)), (2, range(2, 6, 321))))
    sys = uud_system(Δ)
    measure = ssf_custom((q, ssf) -> real(ssf[3, 3]), sys; apply_g=false)
    swt = SpinWaveTheory(sys; measure)
    vacuum = Sunny.MagnonVacuum(swt, Sunny.hartree_fock_correction(swt; tol=1e-4).terms2)
    for (col, dyson) in enumerate((:nambu, :ladder))
        res = Sunny.corrected_intensities(swt, path; energies, η, dyson, vacuum, threaded=true)
        plot_intensities!(fig[row, col], res; colorrange=(0, Δ == 5 ? 0.05 : 0.3),
                          title="S^zz, Δ = $Δ, dyson = :$dyson")
    end
end
fig

# At Δ = 5 the ladder gives four flat two-magnon branches near ω ≈ 9-11 J,
# split off below the continuum: the four branches of the paper's Fig. 4(a).
# At Δ = 2 the lower two branches stay sharp below the continuum, while the
# upper two lie within it and hybridize into broader resonances, as in Fig.
# 4(b). Those that remain sharp inside the continuum are bound states that
# cannot decay, being in a sector of different total Sᶻ.
