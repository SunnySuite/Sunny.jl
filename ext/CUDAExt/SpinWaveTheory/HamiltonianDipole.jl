using LinearAlgebra

function _dot(a, b)
    return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]    
end

function fill_matrix(H, swt, qs_reshaped, qs, L)
    iq = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if iq > size(H, 3)
        return
    end

    Hq = @view H[:, :, iq]
    H11 = @view Hq[1:L, 1:L]
    H12 = @view Hq[1:L, L+1:2L]
    H21 = @view Hq[L+1:2L, 1:L]
    H22 = @view Hq[L+1:2L, L+1:2L]

    q_reshaped = qs_reshaped * qs[iq]

    (; sys, data) = swt
    (; local_rotations, stevens_coefs, sqrtS) = data
    (; extfield, pairs, gs) = sys

    for (i, int) in enumerate(sys.interactions_union)
        # Zeeman term
        B = gs[1, 1, 1, i]' * extfield[1, 1, 1, i]
        B′ = - dot(B, local_rotations[i][:, 3])
        H11[i, i] += B′
        H22[i, i] += B′

        # Single-ion anisotropy
        (; c2, c4, c6) = stevens_coefs[i]
        s = sqrtS[i]^2
        A1 = -6s*c2[3] - 80*s^3*c4[5] - 336*s^5*c6[7]
        A2 = 2s*(c2[1]+im*c2[5]) + 12s^3*(c4[3]+im*c4[7]) + 32s^5*(c6[5]+im*c6[9])
        H11[i, i] += A1
        H22[i, i] += A1
        H12[i, i] += A2
        H21[i, i] += conj(A2)
        
        # Pair interactions
        for idx in int.pair[1]:int.pair[2]
            coupling = pairs[idx]
            (; isculled, bond) = coupling
            isculled && break
            @assert i == bond.i
            j = bond.j

            phase = exp(2π*im * _dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            si = sqrtS[i]^2
            sj = sqrtS[j]^2
            sij = sqrtS[i] * sqrtS[j]
            # Bilinear exchange
            if !iszero(coupling.bilin)
                J = coupling.bilin  # Transformed exchange matrix

                Q = 0.5 * sij * (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1]))
                H11[i, j] += Q * phase
                H11[j, i] += conj(Q) * conj(phase)
                H22[i, j] += conj(Q) * phase
                H22[j, i] += Q  * conj(phase)

                P = 0.5 * sij * (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1]))
                H21[i, j] += P * phase
                H21[j, i] += P * conj(phase)
                H12[i, j] += conj(P) * phase
                H12[j, i] += conj(P) * conj(phase)

                H11[i, i] -= sj * J[3, 3]
                H11[j, j] -= si * J[3, 3]
                H22[i, i] -= sj * J[3, 3]
                H22[j, j] -= si * J[3, 3]
            end

            # Biquadratic exchange
            if !iszero(coupling.biquad)
                K = coupling.biquad  # Transformed quadrupole exchange matrix

                Sj2Si = sj^2 * si
                Si2Sj = si^2 * sj
                H11[i, i] += -12 * Sj2Si * K[3, 3]
                H22[i, i] += -12 * Sj2Si * K[3, 3]
                H11[j, j] += -12 * Si2Sj * K[3, 3]
                H22[j, j] += -12 * Si2Sj * K[3, 3]
                H21[i, i] += 4 * Sj2Si * (K[1, 3] - im*K[5, 3])
                H12[i, i] += 4 * Sj2Si * (K[1, 3] + im*K[5, 3])
                H21[j, j] += 4 * Si2Sj * (K[3, 1] - im*K[3, 5])
                H12[j, j] += 4 * Si2Sj * (K[3, 1] + im*K[3, 5])

                Q = 0.5 * sij^3 * ( K[4, 4]+K[2, 2] - im*(-K[4, 2]+K[2, 4]))
                H11[i, j] += Q * phase
                H11[j, i] += conj(Q * phase)
                H22[i, j] += conj(Q) * phase
                H22[j, i] += Q  * conj(phase)

                P = 0.5 * sij^3 * (-K[4, 4]+K[2, 2] - im*( K[4, 2]+K[2, 4]))
                H21[i, j] += P * phase
                H12[j, i] += conj(P * phase)
                H21[j, i] += P * conj(phase)
                H12[i, j] += conj(P) * phase
            end
        end
    end

    # H must be hermitian up to round-off errors
    @assert Sunny.diffnorm2(Hq, Hq') < 1e-12

    # Make H exactly hermitian
    hermitianpart!(Hq)

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i, i, iq] += swt.regularization
    end
    return
end

function swt_hamiltonian_dipole!(H::CUDA.CuArray{ComplexF64, 3}, swt::SpinWaveTheoryDevice, qs_reshaped, qs::CUDA.CuArray{Sunny.Vec3})
    L = Sunny.nbands(swt)
    Nq = size(qs, 1)
    @assert size(H, 3) == Nq
    @assert size(view(H,:,:,1)) == (2L, 2L)

    H .= 0.0

    kernel = CUDA.@cuda launch=false fill_matrix(H, swt, qs_reshaped, qs, L)
    config = launch_configuration(kernel.fun)
    threads = Base.min(Nq, config.threads)
    blocks = cld(Nq, threads)
    kernel(H, swt, qs_reshaped, qs, L; threads=threads, blocks=blocks)
end