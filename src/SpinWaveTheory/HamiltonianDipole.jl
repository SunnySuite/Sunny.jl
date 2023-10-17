###########################################################################
# Below are the implementations of the SU(N) linear spin-wave calculations #
###########################################################################

# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; R_mat, c_coef) = data
    H .= 0.0

    N = sys.Ns[1]            # Dimension of SU(N) coherent states
    S = (N-1)/2              # Spin magnitude
    L  = natoms(sys.crystal) # Number of quasiparticle bands

    @assert size(H) == (2L, 2L)

    # Zeeman contributions
    (; extfield, gs, units) = sys
    for atom = 1:L
        effB = units.μB * (gs[1, 1, 1, atom]' * extfield[1, 1, 1, atom])
        res = dot(effB, R_mat[atom][:, 3]) / 2
        H[atom, atom]     += res
        H[atom+L, atom+L] += res
    end

    # pairexchange interactions
    for atom = 1:L
        ints = sys.interactions_union[atom]

        # Quadratic exchange
        for coupling in ints.pair
            (; isculled, bond) = coupling
            isculled && break

            J = Mat3(coupling.bilin*I)
            sub_i, sub_j = bond.i, bond.j
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            R_mat_i = R_mat[sub_i]
            R_mat_j = R_mat[sub_j]
            Rij = S * (R_mat_i' * J * R_mat_j)

            P = 0.25 * (Rij[1, 1] - Rij[2, 2] - im*Rij[1, 2] - im*Rij[2, 1])
            Q = 0.25 * (Rij[1, 1] + Rij[2, 2] - im*Rij[1, 2] + im*Rij[2, 1])

            H[sub_i, sub_j] += Q  * phase
            H[sub_j, sub_i] += conj(Q) * conj(phase)
            H[sub_i+L, sub_j+L] += conj(Q) * phase
            H[sub_j+L, sub_i+L] += Q  * conj(phase)

            H[sub_i+L, sub_j] += P * phase
            H[sub_j+L, sub_i] += P * conj(phase)
            H[sub_i, sub_j+L] += conj(P) * phase
            H[sub_j, sub_i+L] += conj(P) * conj(phase)

            H[sub_i, sub_i] -= 0.5 * Rij[3, 3]
            H[sub_j, sub_j] -= 0.5 * Rij[3, 3]
            H[sub_i+L, sub_i+L] -= 0.5 * Rij[3, 3]
            H[sub_j+L, sub_j+L] -= 0.5 * Rij[3, 3]

            ### Biquadratic exchange

            (coupling.biquad isa Number) || error("General biquadratic interactions not yet implemented in LSWT.")
            J = coupling.biquad::Float64

            # ⟨Ω₂, Ω₁|[(𝐒₁⋅𝐒₂)^2 + 𝐒₁⋅𝐒₂/2]|Ω₁, Ω₂⟩ = (Ω₁⋅Ω₂)^2
            # The biquadratic part
            Ri = R_mat[sub_i]
            Rj = R_mat[sub_j]
            Rʳ = Ri' * Rj
            C0 = Rʳ[3, 3]*S^2
            C1 = S*√S/2*(Rʳ[1, 3] + im*Rʳ[2, 3])
            C2 = S*√S/2*(Rʳ[3, 1] + im*Rʳ[3, 2])
            A11 = -Rʳ[3, 3]*S
            A22 = -Rʳ[3, 3]*S
            A21 = S/2*(Rʳ[1, 1] - im*Rʳ[1, 2] - im*Rʳ[2, 1] + Rʳ[2, 2])
            A12 = S/2*(Rʳ[1, 1] + im*Rʳ[1, 2] + im*Rʳ[2, 1] + Rʳ[2, 2])
            B21 = S/4*(Rʳ[1, 1] + im*Rʳ[1, 2] + im*Rʳ[2, 1] - Rʳ[2, 2])
            B12 = B21

            H[sub_i, sub_i] += J* (C0*A11 + C1*conj(C1))
            H[sub_j, sub_j] += J* (C0*A22 + C2*conj(C2))
            H[sub_i, sub_j] += J* ((C0*A12 + C1*conj(C2)) * phase)
            H[sub_j, sub_i] += J* ((C0*A21 + C2*conj(C1)) * conj(phase))
            H[sub_i+L, sub_i+L] += J* (C0*A11 + C1*conj(C1))
            H[sub_j+L, sub_j+L] += J* (C0*A22 + C2*conj(C2))
            H[sub_j+L, sub_i+L] += J* ((C0*A12 + C1*conj(C2)) * conj(phase))
            H[sub_i+L, sub_j+L] += J* ((C0*A21 + C2*conj(C1)) * phase)

            H[sub_i, sub_i+L] += J* (C1*conj(C1))
            H[sub_j, sub_j+L] += J* (C2*conj(C2))
            H[sub_i+L, sub_i] += J* (C1*conj(C1))
            H[sub_j+L, sub_j] += J* (C2*conj(C2))

            H[sub_i, sub_j+L] += J* ((2C0*B12 + C1*C2) * phase)
            H[sub_j, sub_i+L] += J* ((2C0*B21 + C2*C1) * conj(phase))
            H[sub_i+L, sub_j] += J* (conj(2C0*B12 + C1*C2) * phase)
            H[sub_j+L, sub_i] += J* (conj(2C0*B21 + C2*C1) * conj(phase))
        end
    end

    # single-ion anisotropy
    for atom = 1:L
        (; c2, c4, c6) = c_coef[atom]
        H[atom, atom]     += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[atom+L, atom+L] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[atom, atom+L]   += -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        H[atom+L, atom]   +=  im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
    end

    # H must be hermitian up to round-off errors
    if hermiticity_norm(H) > 1e-12
        println("norm(H-H')= ", norm(H-H'))
        throw("H is not hermitian!")
    end
    
    # make H exactly hermitian for cholesky decomposition.
    make_hermitian!(H) 

    # add tiny part to the diagonal elements for cholesky decomposition.
    for i = 1:2L
        H[i, i] += swt.energy_ϵ
    end
end