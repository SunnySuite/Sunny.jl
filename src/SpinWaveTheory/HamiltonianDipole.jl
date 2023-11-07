# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; onsite_hamiltonian) = data

    L = natoms(sys.crystal)  # Number of quasiparticle bands
    @assert size(H) == (2L, 2L)

    # Initialize Hamiltonian with onsite contributions (Zeeman and single-ion)
    H .= onsite_hamiltonian

    # Add pairwise terms 
    for ints in sys.interactions_union

        # Bilinear exchange
        for coupling in ints.pair
            (; isculled, bond) = coupling
            isculled && break
            i, j = bond.i, bond.j
            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            if !iszero(coupling.bilin)
                J = coupling.bilin  # This is Rij in previous notation (transformed exchange matrix)

                P = 0.25 * (J[1, 1] - J[2, 2] - im*J[1, 2] - im*J[2, 1])
                Q = 0.25 * (J[1, 1] + J[2, 2] - im*J[1, 2] + im*J[2, 1])

                H[i, j]     += Q * phase
                H[j, i]     += conj(Q) * conj(phase)
                H[i+L, j+L] += conj(Q) * phase
                H[j+L, i+L] += Q  * conj(phase)

                H[i+L, j] += P * phase
                H[j+L, i] += P * conj(phase)
                H[i, j+L] += conj(P) * phase
                H[j, i+L] += conj(P) * conj(phase)

                H[i, i]     -= 0.5 * J[3, 3]
                H[j, j]     -= 0.5 * J[3, 3]
                H[i+L, i+L] -= 0.5 * J[3, 3]
                H[j+L, j+L] -= 0.5 * J[3, 3]
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                J = coupling.biquad  # Transformed quadrupole exchange matrix
            
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

    # H must be hermitian up to round-off errors
    @assert hermiticity_norm(H) < 1e-12
    
    # Make H exactly hermitian
    hermitianpart!(H) 

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i, i] += swt.energy_ϵ
    end
end