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
    for i in 1:L
        effB = units.μB * (gs[1, 1, 1, i]' * extfield[1, 1, 1, i])
        res = dot(effB, R_mat[i][:, 3]) / 2
        H[i, i]     += res
        H[i+L, i+L] += res
    end

    # pairexchange interactions
    for i in 1:L
        ints = sys.interactions_union[i]

        # Quadratic exchange
        for coupling in ints.pair
            (; isculled, bond) = coupling
            isculled && break

            if !iszero(coupling.bilin)
                J = Mat3(coupling.bilin*I)
                i, j = bond.i, bond.j
                phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

                # Transform couplings according to the local quantization basis
                R_mat_i = R_mat[i]
                R_mat_j = R_mat[j]
                Rij = S * (R_mat_i' * J * R_mat_j)

                P = 0.25 * (Rij[1, 1] - Rij[2, 2] - im*Rij[1, 2] - im*Rij[2, 1])
                Q = 0.25 * (Rij[1, 1] + Rij[2, 2] - im*Rij[1, 2] + im*Rij[2, 1])

                H[i, j]     += Q  * phase
                H[j, i]     += conj(Q) * conj(phase)
                H[i+L, j+L] += conj(Q) * phase
                H[j+L, i+L] += Q  * conj(phase)

                H[i+L, j] += P * phase
                H[j+L, i] += P * conj(phase)
                H[i, j+L] += conj(P) * phase
                H[j, i+L] += conj(P) * conj(phase)

                H[i, i]     -= 0.5 * Rij[3, 3]
                H[j, j]     -= 0.5 * Rij[3, 3]
                H[i+L, i+L] -= 0.5 * Rij[3, 3]
                H[j+L, j+L] -= 0.5 * Rij[3, 3]
            end

            ### Biquadratic exchange

            if !iszero(coupling.biquad)
                J = coupling.biquad
                J = Mat5(J isa Number ? J * diagm(scalar_biquad_metric) : J)

                # Transform couplings according to the local quantization basis
                Vi = operator_for_stevens_rotation(2, R_mat_i)
                Vj = operator_for_stevens_rotation(2, R_mat_j)
                J = S^3 * Mat5(Vi' * J * Vj)
            
                H[i, i] += -6J[3, 3]
                H[j, j] += -6J[3, 3]
                H[i+L, i+L] += -6J[3, 3]
                H[j+L, j+L] += -6J[3, 3]
                H[i+L, i] += 12*(J[1, 3] - im*J[5, 3])
                H[i, i+L] += 12*(J[1, 3] + im*J[5, 3])
                H[j+L, j] += 12*(J[3, 1] - im*J[3, 5])
                H[j, j+L] += 12*(J[3, 1] + im*J[3, 5])

                P = 0.25 * (-J[4, 4]+J[2, 2] - im*( J[4, 2]+J[2, 4]))
                Q = 0.25 * ( J[4, 4]+J[2, 2] - im*(-J[4, 2]+J[2, 4]))

                H[i, j] += Q * phase
                H[j, i] += conj(Q) * conj(phase)
                H[i+L, j+L] += conj(Q) * phase
                H[j+L, i+L] += Q  * conj(phase)

                H[i+L, j] += P * phase
                H[j+L, i] += P * conj(phase)
                H[i, j+L] += conj(P) * phase
                H[j, i+L] += conj(P) * conj(phase)
            end
        end
    end

    # single-ion anisotropy
    for i in 1:L
        (; c2, c4, c6) = c_coef[i]
        H[i, i]     += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[i+L, i+L] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[i, i+L]   += -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        H[i+L, i]   +=  im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
    end

    # H must be hermitian up to round-off errors
    @assert hermiticity_norm(H) < 1e-12
    
    # Make H exactly hermitian
    hermitianpart!(H) 

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i, i] += swt.energy_ϵ
    end
end
