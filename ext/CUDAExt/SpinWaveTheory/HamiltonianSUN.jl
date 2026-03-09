# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode. 

@inline δ(x, y) = (x==y)

function fill_matrix(H11, H12, H21, H22, swt, qs_reshaped, qs)
    iq = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if iq > size(H11, 5)
        return
    end

    (; sys) = swt
    (; pairs, indices, onsite, general) = sys
    N = sys.Ns
    q_reshaped = qs_reshaped * qs[iq]

    @inbounds begin
        for i in 1:length(indices) - 1 
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
            for idx in indices[i]:indices[i+1] - 1    
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
    end
    return
end

function fill_matrix_ewald(H11, H12, H21, H22, A0, Aqs, swt, Na)
    iq = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if iq > size(H, 3)
        return
    end

    (; sys, data) = swt
    (; spins_localized) = data
    (; gs) = sys

    @assert !isnothing(sys.ewald)

    N = sys.Ns

    Aq = @view Aqs[:,:,iq]

    for i in 1:Na, j in 1:Na
        # An ordered pair of magnetic moments contribute (μᵢ A μⱼ)/2 to the
        # energy, where μ = - g S. A symmetric contribution will appear for
        # the bond reversal (i, j) → (j, i).
        J = gs[i]' * Aq[i, j] * gs[j] / 2
        J0 = gs[i]' * A0[i, j] * gs[j] / 2
        for α in 1:3, β in 1:3
            Ai = spins_localized[α, i]
            Bj = spins_localized[β, j]

            for m in 1:N-1, n in 1:N-1
                c = (Ai[m,n] - δ(m,n)*Ai[N,N]) * (Bj[N,N])
                H11[m, i, n, i] += c * J0[α, β]
                H22[n, i, m, i] += c * J0[α, β]

                c = Ai[N,N] * (Bj[m,n] - δ(m,n)*Bj[N,N])
                H11[m, j, n, j] += c * J0[α, β]
                H22[n, j, m, j] += c * J0[α, β]

                c = Ai[m,N] * Bj[N,n]
                H11[m, i, n, j] += c * J[α, β]
                H22[n, j, m, i] += c * conj(J[α, β])

                c = Ai[N,m] * Bj[n,N]
                H11[n, j, m, i] += c * conj(J[α, β])
                H22[m, i, n, j] += c * J[α, β]

                c = Ai[m,N] * Bj[n,N]
                H12[m, i, n, j] += c * J[α, β]
                H12[n, j, m, i] += c * conj(J[α, β])
                H21[n, j, m, i] += conj(c) * conj(J[α, β])
                H21[m, i, n, j] += conj(c) * J[α, β]
            end
        end
    end
    return
end

function swt_hamiltonian_ewald!(H11, H12, H21, H22, swt::SpinWaveTheoryDevice, qs_reshaped, qs::CUDA.CuArray{Sunny.Vec3})
    # Interaction matrix for wavevector (0,0,0). It could be recalculated as:
    # precompute_dipole_ewald(sys.crystal, (1,1,1), demag) * μ0_μB²

    (; demag, μ0_μB², A) = swt.sys.ewald
    Na = Sunny.natoms(swt.sys.crystal)
    Nq = size(qs, 1)
    A0 = reshape(A, Na, Na)

    # Interaction matrix for wavevector q
    A_qs_d = CUDA.zeros(Sunny.SMatrix{3,3,ComplexF64,9}, 1, 1, 1, Na, Na, Nq)
    kernel = CUDA.@cuda launch=false precompute_dipole_ewald_at_wavevector_kernel(A_qs_d, swt.sys.crystal, (1,1,1), demag, qs_reshaped, qs, μ0_μB²)
    config = launch_configuration(kernel.fun)
    threads = Base.min(Nq, config.threads)
    blocks = cld(Nq, threads)
    kernel(A_qs_d, swt.sys.crystal, (1,1,1), demag, qs_reshaped, qs, μ0_μB²; threads=threads, blocks=blocks)
    A_qs_reshape = reshape(A_qs_d, Na, Na, Nq)

    kernel = CUDA.@cuda launch=false fill_matrix_ewald(H11, H12, H21, H22, A0, A_qs_reshape, swt, Na)
    config = launch_configuration(kernel.fun)
    threads = Base.min(Nq, config.threads)
    blocks = cld(Nq, threads)
    kernel(H, A0, A_qs_reshape, swt, Na; threads=threads, blocks=blocks)
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

    if !isnothing(swt.sys.ewald)
        swt_hamiltonian_ewald!(H11, H12, H21, H22, swt, qs_reshaped, qs)
    end

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
