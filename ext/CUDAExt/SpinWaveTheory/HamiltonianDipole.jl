using LinearAlgebra

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
    (; extfield, pairs, indices, gs) = sys

    @inbounds begin
        for (i, int) in enumerate(sys.interactions_union)
            H11_ii = ComplexF64(0.)
            H22_ii = ComplexF64(0.)
            H12_ii = ComplexF64(0.)
            H21_ii = ComplexF64(0.)

            # Zeeman term
            B = gs[1, 1, 1, i]' * extfield[1, 1, 1, i]
            B′ = - dot(B, local_rotations[i][:, 3])
            H11_ii += B′
            H22_ii += B′

            # Single-ion anisotropy
            (; c2, c4, c6) = stevens_coefs[i]
            s = sqrtS[i]^2
            A1 = -6s*c2[3] - 80*s^3*c4[5] - 336*s^5*c6[7]
            A2 = 2s*(c2[1]+im*c2[5]) + 12s^3*(c4[3]+im*c4[7]) + 32s^5*(c6[5]+im*c6[9])
            H11_ii += A1
            H22_ii += A1
            H12_ii += A2
            H21_ii += conj(A2)
        
            # Pair interactions
            for idx in indices[i]:indices[i+1] - 1
                coupling = pairs[idx]
                (; isculled, bond) = coupling
                isculled && break
                @assert i == bond.i
                j = bond.j

                phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

                si = sqrtS[i]^2
                sj = sqrtS[j]^2
                sij = sqrtS[i] * sqrtS[j]

                # Bilinear exchange
                if !iszero(coupling.bilin)
                    J = coupling.bilin  # Transformed exchange matrix

                    Q = 0.5 * sij * (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1]))
                    H11[i,j] += Q * phase
                    H11[j,i] += conj(Q) * conj(phase)
                    H22[i,j] += conj(Q) * phase
                    H22[j,i] += Q  * conj(phase)

                    P = 0.5 * sij * (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1]))
                    H21[i,j] += P * phase
                    H21[j,i] += P * conj(phase)
                    H12[i,j] += conj(P) * phase
                    H12[j,i] += conj(P) * conj(phase)

                    H11_ii -= sj * J[3, 3]
                    H11[j,j] -= si * J[3, 3]
                    H22_ii -= sj * J[3, 3]
                    H22[j,j] -= si * J[3, 3]
                end

                # Biquadratic exchange
                if !iszero(coupling.biquad)
                    K = coupling.biquad  # Transformed quadrupole exchange matrix

                    Sj2Si = sj^2 * si
                    Si2Sj = si^2 * sj
                    H11_ii += -12 * Sj2Si * K[3, 3]
                    H22_ii += -12 * Sj2Si * K[3, 3]
                    H11[j, j] += -12 * Si2Sj * K[3, 3]
                    H22[j, j] += -12 * Si2Sj * K[3, 3]
                    H21_ii += 4 * Sj2Si * (K[1, 3] - im*K[5, 3])
                    H12_ii += 4 * Sj2Si * (K[1, 3] + im*K[5, 3])
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
            H11[i, i] += H11_ii
            H22[i, i] += H22_ii
            H12[i, i] += H12_ii
            H21[i, i] += H21_ii
        end
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

    kernel = CUDA.@cuda launch=false matrix_cleanup(H, swt, L)

    config = launch_configuration(kernel.fun, shmem=threads->get_shmem_matrix_cleanup(threads))
    optimal_threads_1d = config.threads
    threads_x = Base.min(L, optimal_threads_1d)
    threads_y = Base.min(Nq, optimal_threads_1d ÷ threads_x)
    threads = (threads_x, threads_y)

    blocks_x = cld(L, threads_x)
    blocks_y = cld(Nq, threads_y)
    blocks = (blocks_x, blocks_y)
    kernel(H, swt, L; threads=threads, blocks=blocks, shmem=get_shmem_matrix_cleanup(threads))
end
