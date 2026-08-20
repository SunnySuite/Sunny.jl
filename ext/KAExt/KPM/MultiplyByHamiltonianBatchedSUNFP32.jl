@kernel inbounds=true function multiply_by_hamiltonian_SUN_f32_batched_kernel_ka!(
    Y4, X4,
    W_onsite,
    W_out_self, W_out_cross, W_out_anom,
    bond_j_out, bond_n_out, indices_out,
    W_in_self, W_in_cross, W_in_anom,
    bond_i_in, bond_n_in, indices_in,
    q_reshaped_chain,
    regularization,
    Nf, L,
)
    tid      = @index(Global, Linear)
    N_chains = size(Y4, 1)

    chain = ((tid - 1) % N_chains) + 1
    rest  =  (tid - 1) ÷ N_chains
    l     =  (rest % L) + 1
    sigma =  (rest ÷ L) + 1
    f     = (l - 1) % Nf + 1
    me    = (l - 1) ÷ Nf + 1

    CT  = @uniform eltype(Y4)
    pi2 = Float32(2π)
    q   = q_reshaped_chain[chain]

    acc = CT(0)

    if sigma == 1
        for n in 1:Nf
            acc += W_onsite[f, n, me] * X4[chain, n, me, 1]
        end
    else
        for n in 1:Nf
            acc += conj(W_onsite[f, n, me]) * X4[chain, n, me, 2]
        end
    end

    for pidx in indices_out[me] : indices_out[me + 1] - 1
        j     = Int(bond_j_out[pidx])
        nv    = bond_n_out[pidx]
        phase = exp(CT(0, pi2 * (q[1]*nv[1] + q[2]*nv[2] + q[3]*nv[3])))

        if sigma == 1

            for n in 1:Nf
                acc += W_out_self[f, n, pidx] * X4[chain, n, me, 1]
            end

            for n in 1:Nf
                acc += phase * W_out_cross[f, n, pidx] * X4[chain, n, j, 1]
            end

            for n in 1:Nf
                acc += phase * W_out_anom[f, n, pidx] * X4[chain, n, j, 2]
            end
        else

            for n in 1:Nf
                acc += conj(W_out_self[f, n, pidx]) * X4[chain, n, me, 2]
            end

            for n in 1:Nf
                acc += phase * conj(W_out_cross[f, n, pidx]) * X4[chain, n, j, 2]
            end

            for n in 1:Nf
                acc += phase * conj(W_out_anom[f, n, pidx]) * X4[chain, n, j, 1]
            end
        end
    end

    for iidx in indices_in[me] : indices_in[me + 1] - 1
        nbr   = Int(bond_i_in[iidx])
        nv    = bond_n_in[iidx]
        cph   = conj(exp(CT(0, pi2 * (q[1]*nv[1] + q[2]*nv[2] + q[3]*nv[3]))))

        if sigma == 1

            for n in 1:Nf
                acc += W_in_self[f, n, iidx] * X4[chain, n, me, 1]
            end

            for m in 1:Nf
                acc += cph * W_in_cross[f, m, iidx] * X4[chain, m, nbr, 1]
            end

            for m in 1:Nf
                acc += cph * W_in_anom[f, m, iidx] * X4[chain, m, nbr, 2]
            end
        else

            for n in 1:Nf
                acc += conj(W_in_self[f, n, iidx]) * X4[chain, n, me, 2]
            end

            for m in 1:Nf
                acc += cph * conj(W_in_cross[f, m, iidx]) * X4[chain, m, nbr, 2]
            end

            for m in 1:Nf
                acc += cph * conj(W_in_anom[f, m, iidx]) * X4[chain, m, nbr, 1]
            end
        end
    end

    acc += CT(regularization) * X4[chain, f, me, sigma]

    Y4[chain, f, me, sigma] = acc
end

function multiply_by_hamiltonian_SUN_f32_batched_ka!(
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    sys::SystemDeviceSUNF32KA,
    q_reshaped_chain::AbstractVector,
    regularization::Float32,
)
    Nf       = sys.Ns - 1
    Na       = sys.Na
    L        = Nf * Na
    N_chains = size(Y, 1)

    @assert size(Y) == size(X) == (N_chains, 2L) "Bogoliubov vector size mismatch: expected ($N_chains, $(2L))"
    @assert length(q_reshaped_chain) == N_chains  "q_reshaped_chain length must equal N_chains"

    fill!(Y, zero(eltype(Y)))

    Y4 = reshape(Y, (N_chains, Nf, Na, 2))
    X4 = reshape(X, (N_chains, Nf, Na, 2))

    total     = N_chains * L * 2
    backend   = KernelAbstractions.get_backend(Y)
    workgroup = max(min(total, 256), 1)
    kernel!   = multiply_by_hamiltonian_SUN_f32_batched_kernel_ka!(backend, workgroup)
    kernel!(Y4, X4,
            sys.W_onsite,
            sys.W_out_self, sys.W_out_cross, sys.W_out_anom,
            sys.bond_j_out, sys.bond_n_out, sys.indices_out,
            sys.W_in_self, sys.W_in_cross, sys.W_in_anom,
            sys.bond_i_in, sys.bond_n_in, sys.indices_in,
            q_reshaped_chain,
            regularization,
            Nf, L;
            ndrange = total)
    return Y
end

function mul_dynamical_matrix_batched_ka!(
    swt::SpinWaveTheoryDeviceKA{TSys, TData},
    _incoming::IncomingBondCSRKA,
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    q_reshaped_chain::AbstractVector,
) where {TSys <: SystemDeviceSUNF32KA, TData}
    sys_f32 = swt.sys :: SystemDeviceSUNF32KA
    multiply_by_hamiltonian_SUN_f32_batched_ka!(Y, X, sys_f32, q_reshaped_chain,
                                                Float32(swt.regularization))
    return Y
end
