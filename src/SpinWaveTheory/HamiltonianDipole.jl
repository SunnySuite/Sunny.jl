# Adds the long-range dipole-dipole interaction to the dynamical matrix `H` in
# dipole mode. Everything else is compiled into the quadratic table.
function accum_ewald_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_supercell::Vec3)
    (; sys, data) = swt
    (; local_rotations, sqrtS) = data
    (; gs) = sys
    L = nbands(swt)

    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    H22 = view(H, L+1:2L, L+1:2L)
    Rs = local_rotations

    # Interaction matrices between sites at wavevector q and at zero
    Aq = ewald_interaction_tensor(sys, q_supercell)
    A0 = ewald_interaction_tensor(sys, zero(Vec3))

    # Loop over sublattice pairs
    for i in 1:L, j in 1:L
        # An ordered pair of magnetic moments contribute (μᵢ A μⱼ)/2 to the
        # energy. A symmetric contribution will appear for the bond reversal
        # (i, j) → (j, i).  Note that μ = -μB g S.
        J = gs[i]' * ewald_pair(Aq, sys, q_supercell, i, j) * gs[j] / 2
        J0 = gs[i]' * ewald_pair(A0, sys, zero(Vec3), i, j) * gs[j] / 2

        # Perform same transformation as appears in usual bilinear exchange.
        # Rⱼ denotes a rotation from ẑ to the ground state dipole Sⱼ.
        J = sqrtS[i]*sqrtS[j] * Rs[i]' * J * Rs[j]
        J0 = sqrtS[i]*sqrtS[j] * Rs[i]' * J0 * Rs[j]

        # Interactions for Jˣˣ, Jʸʸ, Jˣʸ, and Jʸˣ at wavevector q.
        Q⁻ = 0.5 * (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1]))
        Q⁺ = 0.5 * (J[1, 1] + J[2, 2] + im*(J[1, 2] - J[2, 1]))
        H11[i, j] += Q⁻
        H11[j, i] += conj(Q⁻)
        H22[i, j] += Q⁺
        H22[j, i] += conj(Q⁺)

        P⁻ = 0.5 * (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1]))
        P⁺ = 0.5 * (J[1, 1] - J[2, 2] + im*(J[1, 2] + J[2, 1]))
        H21[i, j] += P⁻
        H21[j, i] += conj(P⁺)
        H12[i, j] += P⁺
        H12[j, i] += conj(P⁻)

        # Interactions for Jᶻᶻ at wavevector (0,0,0).
        H11[i, i] -= J0[3, 3]
        H11[j, j] -= J0[3, 3]
        H22[i, i] -= J0[3, 3]
        H22[j, j] -= J0[3, 3]
    end
end
