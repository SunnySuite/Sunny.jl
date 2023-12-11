# Set the dynamical quadratic Hamiltonian matrix in dipole mode. 
function swt_hamiltonian_dipole!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; local_rotations, stevens_coefs) = data

    N = swt.sys.Ns[1]
    S = (N-1)/2
    L = nbands(swt) 
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
# Delete this before mergining. For benchmarking.
function multiply_by_hamiltonian_dipole_ref!(y, x, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; local_rotations, stevens_coefs) = data

    N = swt.sys.Ns[1]
    S = (N-1)/2
    L = nbands(swt)

    y .= 0

    # Add Zeeman term (can be precalculated)
    (; extfield, gs, units) = sys
    for i in 1:L
        B = units.μB * (gs[1, 1, 1, i]' * extfield[1, 1, 1, i]) 
        B′ = dot(B, local_rotations[i][:, 3]) / 2 

        y[i] += B′ * x[i]
        y[i+L] += B′ * x[i+L]
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

                y[i] += Q * phase * x[j]
                y[j] += conj(Q) * conj(phase) * x[i]
                y[i+L] += conj(Q) * phase * x[j+L]
                y[j+L] += Q * conj(phase) * x[i+L]

                y[i+L] += P * phase * x[j]
                y[j+L] += P * conj(phase) * x[i]
                y[i] += conj(P) * phase * x[j+L]
                y[j] += conj(P) * conj(phase) * x[i+L]

                y[i] -= 0.5 * J[3, 3] * x[i]
                y[j] -= 0.5 * J[3, 3] * x[j]
                y[i+L] -= 0.5 * J[3, 3] * x[i+L]
                y[j+L] -= 0.5 * J[3, 3] * x[j+L]
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                J = coupling.biquad  # Transformed quadrupole exchange matrix
            
                y[i] += -6J[3, 3] * x[i]
                y[j] += -6J[3, 3] * x[j]

                y[i+L] += -6J[3, 3] * x[i+L]
                y[j+L] += -6J[3, 3] * x[j+L]

                y[i+L] += 12*(J[1, 3] - im*J[5, 3]) * x[i]
                y[i] += 12*(J[1, 3] + im*J[5, 3]) * x[i+L]
                y[j+L] += 12*(J[3, 1] - im*J[3, 5]) * x[j]
                y[j] += 12*(J[3, 1] + im*J[3, 5]) * x[j+L]

                P = 0.25 * (-J[4, 4]+J[2, 2] - im*( J[4, 2]+J[2, 4]))
                Q = 0.25 * ( J[4, 4]+J[2, 2] - im*(-J[4, 2]+J[2, 4]))

                y[i] += Q * phase * x[j]
                y[j] += conj(Q) * conj(phase) * x[i]
                y[i+L] += conj(Q) * phase * x[j+L]
                y[j+L] += Q * conj(phase) * x[i+L]

                y[i+L] += P * phase * x[j]
                y[j+L] += P * conj(phase) * x[i]
                y[i] += conj(P) * phase * x[j+L]
                y[j] += conj(P) * conj(phase) * x[i+L]
            end
        end
    end

    # Add single-ion anisotropy
    for i in 1:L
        (; c2, c4, c6) = stevens_coefs[i]
        y[i] += (-3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]) * x[i]
        y[i+L] += (-3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]) * x[i+L]
        y[i] += (-im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])) * x[i+L]
        y[i+L] += (im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])) * x[i]
    end

    # Add small constant shift for positive-definiteness
    @. y += swt.energy_ϵ * x

    nothing
end

function multiply_by_hamiltonian_dipole(x::Array{ComplexF64, 2}, swt::SpinWaveTheory, qs_reshaped::Array{Vec3})
    # Preallocate buffers
    y = zeros(ComplexF64, (size(qs_reshaped)..., size(x, 2)))
    phasebuf = zeros(ComplexF64, length(qs_reshaped))

    # Precompute e^{2πq_α} components
    qphase = map(qs_reshaped) do q  
        (exp(2π*im*q[1]), exp(2π*im*q[2]), exp(2π*im*q[3]))
    end

    # Perform batched matrix-vector multiply
    multiply_by_hamiltonian_dipole_aux!(reshape(y, (length(qs_reshaped), length(x))), x, phasebuf, qphase, swt)

    return y 
end

# function multiply_by_hamiltonian_dipole_aux!(y::Array{ComplexF64, 2}, x::Array{ComplexF64, 2}, phasebuf::Vector{ComplexF64}, qphase, swt::SpinWaveTheory)
function multiply_by_hamiltonian_dipole_aux!(y, x, phasebuf, qphase, swt)
    (; sys, data) = swt
    (; stevens_coefs, local_rotations) = data

    N = swt.sys.Ns[1]
    S = (N-1)/2
    L = natoms(sys.crystal) 

    y .= 0
    nq = size(y, 1) 
    x = Base.ReshapedArray(x, (nq, L, 2), ())
    y = Base.ReshapedArray(y, (nq, L, 2), ())

    # Add single-site terms, both Zeeman and single-ion anisotropy.
    # These entries are q-independent and could be precomputed. I tried this --
    # it was actually very marginally slower. I.e., memory lookup actually costs
    # more than online computation of coefs and B′.
    (; units, extfield, gs) = sys
    for i in 1:L
        (; c2, c4, c6) = stevens_coefs[i]

        coef1 = -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        coef2 = -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        coef3 = -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        coef4 = im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])

        B = units.μB * (gs[1, 1, 1, i]' * extfield[1, 1, 1, i]) 
        B′ = dot(B, view(local_rotations[i], :, 3)) / 2 

        # Seems to be no benefit to breaking this into two loops acting on different final indices.
        for q in 1:nq
            y[q, i, 1] += (B′ + coef1) * x[q, i, 1] + coef3 * x[q, i, 2]
            y[q, i, 2] += (B′ + coef2) * x[q, i, 2] + coef4 * x[q, i, 1]
        end
    end

    # Add pairwise terms 
    for ints in sys.interactions_union

        # Bilinear exchange
        for coupling in ints.pair
            (; isculled, bond) = coupling
            isculled && break
            i, j = bond.i, bond.j

            # Calculate phase associated with periodic wrapping
            n1, n2, n3 = bond.n
            map!(qp -> (qp[1]^n1)*(qp[2]^n2)*(qp[3]^n3), phasebuf, qphase)

            if !iszero(coupling.bilin)
                J = coupling.bilin  # This is Rij in previous notation (transformed exchange matrix)

                P = 0.25 * (J[1, 1] - J[2, 2] - im*J[1, 2] - im*J[2, 1])
                Q = 0.25 * (J[1, 1] + J[2, 2] - im*J[1, 2] + im*J[2, 1])

                @inbounds begin
                    for q in axes(y, 1) 
                        y[q, i, 1] += Q * phasebuf[q] * x[q, j, 1]
                        y[q, i, 1] += conj(P) * phasebuf[q] * x[q, j, 2]
                        y[q, i, 1] -= 0.5 * J[3, 3] * x[q, i, 1]
                    end
                    for q in axes(y, 1) 
                        y[q, i, 2] += conj(Q) * phasebuf[q] * x[q, j, 2]
                        y[q, i, 2] += P * phasebuf[q] * x[q, j, 1]
                        y[q, i, 2] -= 0.5 * J[3, 3] * x[q, i, 2]
                    end
                    for q in axes(y, 1) 
                        y[q, j, 1] += conj(P) * conj(phasebuf[q]) * x[q, i, 2]
                        y[q, j, 1] += conj(Q) * conj(phasebuf[q]) * x[q, i, 1]
                        y[q, j, 1] -= 0.5 * J[3, 3] * x[q, j, 1]
                    end
                    for q in axes(y, 1) 
                        y[q, j, 2] += Q * conj(phasebuf[q]) * x[q, i, 2]
                        y[q, j, 2] += P * conj(phasebuf[q]) * x[q, i, 1]
                        y[q, j, 2] -= 0.5 * J[3, 3] * x[q, j, 2]
                    end
                end

                # for q in axes(y, 1) 
                #     y[q, i, 1] += Q * phasebuf[q] * x[q, j, 1]
                #     y[q, j, 1] += conj(Q) * conj(phasebuf[q]) * x[q, i, 1]
                #     y[q, i, 2] += conj(Q) * phasebuf[q] * x[q, j, 2]
                #     y[q, j, 2] += Q * conj(phasebuf[q]) * x[q, i, 2]

                #     y[q, i, 2] += P * phasebuf[q] * x[q, j, 1]
                #     y[q, j, 2] += P * conj(phasebuf[q]) * x[q, i, 1]
                #     y[q, i, 1] += conj(P) * phasebuf[q] * x[q, j, 2]
                #     y[q, j, 1] += conj(P) * conj(phasebuf[q]) * x[q, i, 2]

                #     y[q, i, 1] -= 0.5 * J[3, 3] * x[q, i, 1]
                #     y[q, j, 1] -= 0.5 * J[3, 3] * x[q, j, 1]
                #     y[q, i, 2] -= 0.5 * J[3, 3] * x[q, i, 2]
                #     y[q, j, 2] -= 0.5 * J[3, 3] * x[q, j, 2]
                # end
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                J = coupling.biquad  # Transformed quadrupole exchange matrix

                P = 0.25 * (-J[4, 4]+J[2, 2] - im*( J[4, 2]+J[2, 4]))
                Q = 0.25 * ( J[4, 4]+J[2, 2] - im*(-J[4, 2]+J[2, 4]))

                @inbounds for q in 1:nq
                    y[q, i, 1] += -6J[3, 3] * x[q, i, 1]
                    y[q, i, 1] += 12*(J[1, 3] + im*J[5, 3]) * x[q, i, 2]
                    y[q, i, 1] += Q * phasebuf[q] * x[q, j, 1]
                    y[q, i, 1] += conj(P) * phasebuf[q] * x[q, j, 2]
                end
                @inbounds for q in 1:nq
                    y[q, i, 2] += -6J[3, 3] * x[q, i, 2]
                    y[q, i, 2] += 12*(J[1, 3] - im*J[5, 3]) * x[q, i, 1]
                    y[q, i, 2] += conj(Q) * phasebuf[q] * x[q, j, 2]
                    y[q, i, 2] += P * phasebuf[q] * x[q, j, 1]
                end
                @inbounds for q in 1:nq
                    y[q, j, 1] += -6J[3, 3] * x[q, j, 1]
                    y[q, j, 1] += 12*(J[3, 1] + im*J[3, 5]) * x[q, j, 2]
                    y[q, j, 1] += conj(Q) * conj(phasebuf[q]) * x[q, i, 1]
                    y[q, j, 1] += conj(P) * conj(phasebuf[q]) * x[q, i, 2]
                end
                @inbounds for q in 1:nq
                    y[q, j, 2] += -6J[3, 3] * x[q, j, 2]
                    y[q, j, 2] += 12*(J[3, 1] - im*J[3, 5]) * x[q, j, 1]
                    y[q, j, 2] += Q * conj(phasebuf[q]) * x[q, i, 2]
                    y[q, j, 2] += P * conj(phasebuf[q]) * x[q, i, 1]
                end

                # P = 0.25 * (-J[4, 4]+J[2, 2] - im*( J[4, 2]+J[2, 4]))
                # Q = 0.25 * ( J[4, 4]+J[2, 2] - im*(-J[4, 2]+J[2, 4]))
                # for q in 1:nq
                #     y[q, i, 1] += -6J[3, 3] * x[q, i, 1]
                #     y[q, j, 1] += -6J[3, 3] * x[q, j, 1]

                #     y[q, i, 2] += -6J[3, 3] * x[q, i, 2]
                #     y[q, j, 2] += -6J[3, 3] * x[q, j, 2]

                #     y[q, i, 2] += 12*(J[1, 3] - im*J[5, 3]) * x[q, i, 1]
                #     y[q, i, 1] += 12*(J[1, 3] + im*J[5, 3]) * x[q, i, 2]
                #     y[q, j, 2] += 12*(J[3, 1] - im*J[3, 5]) * x[q, j, 1]
                #     y[q, j, 1] += 12*(J[3, 1] + im*J[3, 5]) * x[q, j, 2]

                #     y[q, i, 1] += Q * phasebuf[q] * x[q, j, 1]
                #     y[q, j, 1] += conj(Q) * conj(phasebuf[q]) * x[q, i, 1]
                #     y[q, i, 2] += conj(Q) * phasebuf[q] * x[q, j, 2]
                #     y[q, j, 2] += Q * conj(phasebuf[q]) * x[q, i, 2]

                #     y[q, i, 2] += P * phasebuf[q] * x[q, j, 1]
                #     y[q, j, 2] += P * conj(phasebuf[q]) * x[q, i, 1]
                #     y[q, i, 1] += conj(P) * phasebuf[q] * x[q, j, 2]
                #     y[q, j, 1] += conj(P) * conj(phasebuf[q]) * x[q, i, 2]
                # end
            end
        end
    end

    # Add small constant shift for positive-definiteness. 
    @. y += swt.energy_ϵ * x

    nothing
end