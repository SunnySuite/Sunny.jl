@kernel inbounds=true function multiply_by_hamiltonian_SUN_batched_kernel_ka!(
    Y4, X4,
    onsite,
    general,
    pairs,
    indices,
    incoming_pairs,
    incoming_indices,
    pair_map,
    q_reshaped_chain,
    regularization,
    Nf, N, L,
)
    tid      = @index(Global, Linear)
    N_chains = size(Y4, 1)

    chain = ((tid - 1) % N_chains) + 1
    rest  =  (tid - 1) ÷ N_chains
    l     =  (rest % L) + 1
    sigma =  (rest ÷ L) + 1
    f     = (l - 1) % Nf + 1
    me    = (l - 1) ÷ Nf + 1

    CT = eltype(Y4)
    T  = real(CT)

    q   = q_reshaped_chain[chain]
    Ngd = size(general, 4)

    acc = CT(0)

    OSUM = onsite[N, N, me]
    if sigma == 1
        for n in 1:Nf
            c = onsite[f, n, me] - (f == n ? OSUM : CT(0))
            acc += c * X4[chain, n, me, 1]
        end
    else
        for n in 1:Nf
            c = onsite[n, f, me] - (n == f ? OSUM : CT(0))
            acc += c * X4[chain, n, me, 2]
        end
    end

    for pidx in indices[me] : indices[me+1] - 1
        coupling = pairs[pidx]
        coupling.isculled && break
        j     = coupling.bond.j
        phase = exp(T(2π) * im * dot(q, coupling.bond.n))

        if sigma == 1
            for k in 1:Ngd
                ANN = general[N, N, 1, k, pidx]
                BNN = general[N, N, 2, k, pidx]
                AfN = general[f, N, 1, k, pidx]

                for n in 1:Nf
                    c1 = (general[f, n, 1, k, pidx] - (f == n ? ANN : CT(0))) * BNN
                    acc += c1 * X4[chain, n, me, 1]
                end

                for n in 1:Nf
                    c3 = AfN * general[N, n, 2, k, pidx]
                    acc += c3 * phase * X4[chain, n, j, 1]
                end

                for n in 1:Nf
                    c5 = AfN * general[n, N, 2, k, pidx]
                    acc += c5 * phase * X4[chain, n, j, 2]
                end
            end
        else
            for k in 1:Ngd
                ANN = general[N, N, 1, k, pidx]
                BNN = general[N, N, 2, k, pidx]
                ANf = general[N, f, 1, k, pidx]
                AfN = general[f, N, 1, k, pidx]

                for n in 1:Nf
                    c1 = (general[n, f, 1, k, pidx] - (n == f ? ANN : CT(0))) * BNN
                    acc += c1 * X4[chain, n, me, 2]
                end

                for n in 1:Nf
                    c4 = ANf * general[n, N, 2, k, pidx]
                    acc += c4 * phase * X4[chain, n, j, 2]
                end

                for n in 1:Nf
                    c5 = AfN * general[n, N, 2, k, pidx]
                    acc += conj(c5) * phase * X4[chain, n, j, 1]
                end
            end
        end
    end

    for iidx in incoming_indices[me] : incoming_indices[me+1] - 1
        coupling = incoming_pairs[iidx]
        nbr   = coupling.bond.i
        phase = exp(T(2π) * im * dot(q, coupling.bond.n))
        oidx  = pair_map[iidx]

        if sigma == 1
            for k in 1:Ngd
                ANN = general[N, N, 1, k, oidx]
                BNN = general[N, N, 2, k, oidx]
                BfN = general[f, N, 2, k, oidx]

                for n in 1:Nf
                    c2 = ANN * (general[f, n, 2, k, oidx] - (f == n ? BNN : CT(0)))
                    acc += c2 * X4[chain, n, me, 1]
                end

                for m in 1:Nf
                    c4 = general[N, m, 1, k, oidx] * BfN
                    acc += c4 * conj(phase) * X4[chain, m, nbr, 1]
                end

                for m in 1:Nf
                    c5 = general[m, N, 1, k, oidx] * BfN
                    acc += c5 * conj(phase) * X4[chain, m, nbr, 2]
                end
            end
        else
            for k in 1:Ngd
                ANN = general[N, N, 1, k, oidx]
                BNN = general[N, N, 2, k, oidx]
                BNf = general[N, f, 2, k, oidx]
                BfN = general[f, N, 2, k, oidx]

                for n in 1:Nf
                    c2 = ANN * (general[n, f, 2, k, oidx] - (n == f ? BNN : CT(0)))
                    acc += c2 * X4[chain, n, me, 2]
                end

                for m in 1:Nf
                    c3 = general[m, N, 1, k, oidx] * BNf
                    acc += c3 * conj(phase) * X4[chain, m, nbr, 2]
                end

                for m in 1:Nf
                    c5 = general[m, N, 1, k, oidx] * BfN
                    acc += conj(c5) * conj(phase) * X4[chain, m, nbr, 1]
                end
            end
        end
    end

    acc += T(regularization) * X4[chain, f, me, sigma]

    Y4[chain, f, me, sigma] = acc
end

function multiply_by_hamiltonian_SUN_batched_ka!(
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    swt::SpinWaveTheoryDeviceKA{TSys, TData},
    incoming::IncomingBondCSRKA,
    q_reshaped_chain::AbstractVector,
) where {TSys <: SystemDeviceSUNKA, TData}
    sys      = swt.sys
    Na       = Sunny.natoms(sys.crystal)
    Nf       = sys.Ns - 1
    N        = sys.Ns
    L        = Nf * Na
    N_chains = size(Y, 1)

    @assert size(Y) == size(X) == (N_chains, 2L) "Bogoliubov vector size mismatch: expected ($N_chains, $(2L)), got $(size(Y)) and $(size(X))"
    @assert length(q_reshaped_chain) == N_chains "q_reshaped_chain length must equal N_chains"
    @assert eltype(Y) == eltype(X) "Y and X must have the same element type"

    fill!(Y, zero(eltype(Y)))

    Y4 = reshape(Y, (N_chains, Nf, Na, 2))
    X4 = reshape(X, (N_chains, Nf, Na, 2))

    total     = N_chains * L * 2
    backend   = KernelAbstractions.get_backend(Y)
    workgroup = max(min(total, 256), 1)
    kernel!   = multiply_by_hamiltonian_SUN_batched_kernel_ka!(backend, workgroup)
    kernel!(Y4, X4,
            sys.onsite, sys.general,
            sys.pairs, sys.indices,
            incoming.incoming_pairs, incoming.incoming_indices, incoming.pair_map,
            q_reshaped_chain,
            Float64(swt.regularization),
            Nf, N, L;
            ndrange=total)
    return Y
end

function mul_dynamical_matrix_batched_ka!(
    swt::SpinWaveTheoryDeviceKA{TSys, TData},
    incoming::IncomingBondCSRKA,
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    q_reshaped_chain::AbstractVector,
) where {TSys <: SystemDeviceSUNKA, TData}
    multiply_by_hamiltonian_SUN_batched_ka!(Y, X, swt, incoming, q_reshaped_chain)
    return Y
end
