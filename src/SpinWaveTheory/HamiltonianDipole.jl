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
    L = natoms(sys.crystal)  # Number of quasiparticle bands

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

            if !iszero(coupling.bilin)
                J = Mat3(coupling.bilin*I)
                sub_i, sub_j = bond.i, bond.j
                phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

                # Transform couplings according to the local quantization basis
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
            end

            ### Biquadratic exchange

            if !iszero(coupling.biquad)
                J = coupling.biquad
                J = Mat5(J isa Number ? J * diagm(scalar_biquad_metric) : J)

                # Transform couplings according to the local quantization basis
                Vi = operator_for_stevens_rotation(2, R_mat_i)
                Vj = operator_for_stevens_rotation(2, R_mat_j)
                J = S^3 * Mat5(Vi' * J * Vj)
            
                H[sub_i, sub_i] += -6J[3, 3]
                H[sub_j, sub_j] += -6J[3, 3]
                H[sub_i+L, sub_i+L] += -6J[3, 3]
                H[sub_j+L, sub_j+L] += -6J[3, 3]
                H[sub_i+L, sub_i] += 12*(J[1, 3] - 1im*J[5, 3])
                H[sub_i, sub_i+L] += 12*(J[1, 3] + 1im*J[5, 3])
                H[sub_j+L, sub_j] += 12*(J[3, 1] - 1im*J[3, 5])
                H[sub_j, sub_j+L] += 12*(J[3, 1] + 1im*J[3, 5])

                P = 0.25 * (-J[4, 4]+J[2, 2] - 1im*( J[4, 2]+J[2, 4]))
                Q = 0.25 * ( J[4, 4]+J[2, 2] - 1im*(-J[4, 2]+J[2, 4]))

                H[sub_i, sub_j] += Q * phase
                H[sub_j, sub_i] += conj(Q) * conj(phase)
                H[sub_i+L, sub_j+L] += conj(Q) * phase
                H[sub_j+L, sub_i+L] += Q  * conj(phase)

                H[sub_i+L, sub_j] += P * phase
                H[sub_j+L, sub_i] += P * conj(phase)
                H[sub_i, sub_j+L] += conj(P) * phase
                H[sub_j, sub_i+L] += conj(P) * conj(phase)
            end
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
    @assert hermiticity_norm(H) < 1e-12
    
    # Make H exactly hermitian
    hermitianpart!(H) 

    # Add small constant shift for positive-definiteness
    for i = 1:2L
        H[i, i] += swt.energy_ϵ
    end
end
