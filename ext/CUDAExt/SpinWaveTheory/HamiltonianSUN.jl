# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 

@inline δ(x, y) = (x==y)

function fill_matrix(H11, H12, H21, H22, swt, qs_reshaped, qs)
    iq = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if iq > size(H11, 5)
        return
    end

    (; sys) = swt
    (; pairs, onsite, general) = sys
    N = sys.Ns
    q_reshaped = qs_reshaped * qs[iq]
    for (i, int) in enumerate(sys.interactions_union)
        # Onsite coupling, including Zeeman. Note that op has already been
        # transformed according to the local frame of sublattice i.
        op = view(onsite,:,:,i)
        for m in 1:N-1
            for n in 1:N-1
                c = op[m, n] - δ(m, n) * op[N, N]
                H11[m, i, n, i, iq] += c
                H22[n, i, m, i, iq] += c
            end
        end
        for idx in int.pair[1]:int.pair[2]
            coupling = pairs[idx]
            (; isculled, bond) = coupling
            isculled && break

            @assert i == bond.i
            j = bond.j

            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            # Set "general" pair interactions of the form Aᵢ⊗Bⱼ. Note that Aᵢ
            # and Bᵢ have already been transformed according to the local frames
            # of sublattice i and j, respectively.
            for jdx in 1:size(general,4)
                Ai = view(general,:,:,1,jdx,idx)
                Bj = view(general,:,:,2,jdx,idx)
                for m in 1:N-1, n in 1:N-1
                    c = (Ai[m,n] - δ(m,n)*Ai[N,N]) * (Bj[N,N])
                    H11[m, i, n, i, iq] += c
                    H22[n, i, m, i, iq] += c

                    c = Ai[N,N] * (Bj[m,n] - δ(m,n)*Bj[N,N])
                    H11[m, j, n, j, iq] += c
                    H22[n, j, m, j, iq] += c

                    c = Ai[m,N] * Bj[N,n]
                    H11[m, i, n, j, iq] += c * phase
                    H22[n, j, m, i, iq] += c * conj(phase)

                    c = Ai[N,m] * Bj[n,N]
                    H11[n, j, m, i, iq] += c * conj(phase)
                    H22[m, i, n, j, iq] += c * phase

                    c = Ai[m,N] * Bj[n,N]
                    H12[m, i, n, j, iq] += c * phase
                    H12[n, j, m, i, iq] += c * conj(phase)
                    H21[n, j, m, i, iq] += conj(c) * conj(phase)
                    H21[m, i, n, j, iq] += conj(c) * phase
                end
            end
        end
    end
    return
end


function matrix_cleanup(H, swt, L)
    iq = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if iq > size(H, 3)
        return
    end

    Hq = view(H,:,:,iq)

    # H must be hermitian up to round-off errors
    @assert Sunny.diffnorm2(Hq, Hq') < 1e-12

    # Make H exactly hermitian
    hermitianpart!(Hq)

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        Hq[i,i] += swt.regularization
    end
    return
end

function swt_hamiltonian_SUN!(H::CUDA.CuArray{ComplexF64,3}, swt::SpinWaveTheoryDevice, qs_reshaped, qs::CUDA.CuArray{Sunny.Vec3})
    (; sys) = swt
    N = sys.Ns
    Na = Sunny.natoms(sys.crystal)
    L = (N - 1) * Na

    Nq = size(qs, 1)
    @assert size(H, 3) == Nq
    @assert size(view(H, :, :, 1)) == (2L, 2L)

    H .= 0.0
    blockdims = (N - 1, Na, N - 1, Na, Nq)

    H11 = reshape(view(H, 1:L, 1:L, :), blockdims)
    H12 = reshape(view(H, 1:L, L+1:2L, :), blockdims)
    H21 = reshape(view(H, L+1:2L, 1:L, :), blockdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L, :), blockdims)

    kernel = CUDA.@cuda launch=false fill_matrix(H11, H12, H21, H22, swt, qs_reshaped, qs)
    config = launch_configuration(kernel.fun)
    threads = Base.min(Nq, config.threads)
    blocks = cld(Nq, threads)
    kernel(H11, H12, H21, H22, swt, qs_reshaped, qs; threads=threads, blocks=blocks)

    kernel = CUDA.@cuda launch=false matrix_cleanup(H, swt, L)
    config = launch_configuration(kernel.fun)
    threads = Base.min(Nq, config.threads)
    blocks = cld(Nq, threads)
    kernel(H, swt, L; threads=threads, blocks=blocks)
end
