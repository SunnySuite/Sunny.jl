# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; local_rotations, stevens_coefs) = data

    N = swt.sys.Ns[1]
    S = (N-1)/2
    L = num_bands(swt) 
    @assert size(H) == (2L, 2L)

    # Initialize Hamiltonian buffer 
    H .= 0.0 
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    H22 = view(H, L+1:2L, L+1:2L)

    # Add Zeeman term
    (; extfield, gs, units) = sys
    for i in 1:L
        B = units.μB * (gs[1, 1, 1, i]' * extfield[1, 1, 1, i]) 
        B′ = dot(B, local_rotations[i][:, 3]) / 2 
        H11[i, i] += B′
        H22[i, i] += B′
    end

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

                H11[i, j] += Q * phase
                H11[j, i] += conj(Q) * conj(phase)
                H22[i, j] += conj(Q) * phase
                H22[j, i] += Q  * conj(phase)

                H21[i, j] += P * phase
                H21[j, i] += P * conj(phase)
                H12[i, j] += conj(P) * phase
                H12[j, i] += conj(P) * conj(phase)

                H11[i, i] -= 0.5 * J[3, 3]
                H11[j, j] -= 0.5 * J[3, 3]
                H22[i, i] -= 0.5 * J[3, 3]
                H22[j, j] -= 0.5 * J[3, 3]
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                J = coupling.biquad  # Transformed quadrupole exchange matrix
            
                H11[i, i] += -6J[3, 3]
                H11[j, j] += -6J[3, 3]
                H22[i, i] += -6J[3, 3]
                H22[j, j] += -6J[3, 3]
                H21[i, i] += 12*(J[1, 3] - im*J[5, 3])
                H12[i, i] += 12*(J[1, 3] + im*J[5, 3])
                H21[j, j] += 12*(J[3, 1] - im*J[3, 5])
                H12[j, j] += 12*(J[3, 1] + im*J[3, 5])

                P = 0.25 * (-J[4, 4]+J[2, 2] - im*( J[4, 2]+J[2, 4]))
                Q = 0.25 * ( J[4, 4]+J[2, 2] - im*(-J[4, 2]+J[2, 4]))

                H11[i, j] += Q * phase
                H11[j, i] += conj(Q) * conj(phase)
                H22[i, j] += conj(Q) * phase
                H22[j, i] += Q  * conj(phase)

                H21[i, j] += P * phase
                H21[j, i] += P * conj(phase)
                H12[i, j] += conj(P) * phase
                H12[j, i] += conj(P) * conj(phase)
            end
        end
    end

    # Add single-ion anisotropy
    for i in 1:L
        (; c2, c4, c6) = stevens_coefs[i]
        H11[i, i] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H22[i, i] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H12[i, i] += -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        H21[i, i] +=  im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
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


# Calculate y = H*x, where H is the quadratic Hamiltonian matrix (dynamical matrix).
function multiply_by_hamiltonian_dipole!(y, x, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; local_rotations, stevens_coefs) = data

    N = swt.sys.Ns[1]
    S = (N-1)/2
    L = num_bands(swt)

    x = reshape(x, natoms(sys.crystal), 2)
    y = reshape(y, natoms(sys.crystal), 2)

    # Add Zeeman term
    (; extfield, gs, units) = sys
    for i in 1:L
        B = units.μB * (gs[1, 1, 1, i]' * extfield[1, 1, 1, i]) 
        B′ = dot(B, local_rotations[i][:, 3]) / 2 

        y[i, 1] += B′ * x[i, 1]
        y[i, 2] += B′ * x[i, 2]
    end

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

                y[i, 1] += Q * phase * x[j, 1]
                y[j, 1] += conj(Q) * conj(phase) * x[i, 1]
                y[i, 2] += conj(Q) * phase * x[j, 2]
                y[j, 2] += Q * conj(phase) * x[i, 2]

                y[i, 2] += P * phase * x[j, 1]
                y[j, 2] += P * conj(phase) * x[i, 1]
                y[i, 1] += conj(P) * phase * x[j, 2]
                y[j, 1] += conj(P) * conj(phase) * x[i, 2]

                y[i, 1] -= 0.5 * J[3, 3] * x[i, 1]
                y[j, 1] -= 0.5 * J[3, 3] * x[j, 1]
                y[i, 2] -= 0.5 * J[3, 3] * x[i, 2]
                y[j, 2] -= 0.5 * J[3, 3] * x[j, 2]
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                J = coupling.biquad  # Transformed quadrupole exchange matrix
            
                y[i, 1] += -6J[3, 3] * x[i, 1]
                y[j, 1] += -6J[3, 3] * x[j, 1]

                y[i, 2] += -6J[3, 3] * x[i, 2]
                y[j, 2] += -6J[3, 3] * x[j, 2]

                y[i, 2] += 12*(J[1, 3] - im*J[5, 3]) * x[i, 1]
                y[i, 1] += 12*(J[1, 3] + im*J[5, 3]) * x[i, 2]
                y[j, 2] += 12*(J[3, 1] - im*J[3, 5]) * x[j, 1]
                y[j, 1] += 12*(J[3, 1] + im*J[3, 5]) * x[j, 2]

                P = 0.25 * (-J[4, 4]+J[2, 2] - im*( J[4, 2]+J[2, 4]))
                Q = 0.25 * ( J[4, 4]+J[2, 2] - im*(-J[4, 2]+J[2, 4]))

                y[i, 1] += Q * phase * x[j, 1]
                y[j, 1] += conj(Q) * conj(phase) * x[i, 1]
                y[i, 2] += conj(Q) * phase * x[j, 2]
                y[j, 2] += Q * conj(phase) * x[i, 2]

                y[i, 2] += P * phase * x[j, 1]
                y[j, 2] += P * conj(phase) * x[i, 1]
                y[i, 1] += conj(P) * phase * x[j, 2]
                y[j, 1] += conj(P) * conj(phase) * x[i, 2]
            end
        end
    end

    # Add single-ion anisotropy
    for i in 1:L
        (; c2, c4, c6) = stevens_coefs[i]
        y[i, 1] += (-3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]) * x[i, 1]
        y[i, 2] += (-3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]) * x[i, 2]
        y[i, 1] += (-im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])) * x[i, 2]
        y[i, 2] += (im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])) * x[i, 1]
    end

    # Add small constant shift for positive-definiteness
    @. y += swt.energy_ϵ * x

    nothing
end