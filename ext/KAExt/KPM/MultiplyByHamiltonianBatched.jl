@inline pow5_ka(x::Float32) = x*x*x*x*x
@inline pow5_ka(x) = x^5

@kernel inbounds=true function multiply_by_hamiltonian_dipole_batched_kernel_ka!(
    Y, X, sys, incoming, data,
    q_reshaped_chain, regularization, L,
)
    tid = @index(Global, Linear)
    N_chains = size(Y, 1)

    chain_idx = ((tid - 1) % N_chains) + 1
    rest      =  (tid - 1) ÷ N_chains
    me        =  (rest % L) + 1
    sigma     =  (rest ÷ L) + 1

    q = q_reshaped_chain[chain_idx]

    (; local_rotations, stevens_coefs, sqrtS) = data
    (; extfield, pairs, indices, gs) = sys

    T = real(eltype(Y))
    CT = Complex{T}

    acc = CT(0)

    s_me = sqrtS[me]^2
    c2_me = stevens_coefs[me].c2
    c4_me = stevens_coefs[me].c4
    c6_me = stevens_coefs[me].c6
    A1 = -6*s_me*c2_me[3] - 80*s_me^3*c4_me[5] - 336*pow5_ka(s_me)*c6_me[7]
    A2 = 2*s_me*(c2_me[1] + im*c2_me[5]) + 12*s_me^3*(c4_me[3] + im*c4_me[7]) + 32*pow5_ka(s_me)*(c6_me[5] + im*c6_me[9])

    B = gs[1, 1, 1, me]' * extfield[1, 1, 1, me]
    Bprime = -dot(B, local_rotations[me][:, 3])

    if sigma == 1
        acc += T(Bprime + A1) * X[chain_idx, me, 1]
        acc += CT(A2) * X[chain_idx, me, 2]
    else
        acc += T(Bprime + A1) * X[chain_idx, me, 2]
        acc += CT(conj(A2)) * X[chain_idx, me, 1]
    end

    for idx in indices[me]:indices[me + 1] - 1
        coupling = pairs[idx]
        (; isculled, bond) = coupling
        isculled && break
        j = bond.j

        phase = exp(T(2π)*im * dot(q, bond.n))

        si_local  = sqrtS[me]^2
        sj_local  = sqrtS[j]^2
        sij_local = sqrtS[me] * sqrtS[j]

        if !iszero(coupling.bilin)
            J = coupling.bilin
            P = T(0.5) * sij_local * (J[1, 1] - J[2, 2] - im*J[1, 2] - im*J[2, 1])
            Q = T(0.5) * sij_local * (J[1, 1] + J[2, 2] - im*J[1, 2] + im*J[2, 1])

            if sigma == 1
                acc += CT(Q * phase) * X[chain_idx, j, 1]
                acc += CT(conj(P) * phase) * X[chain_idx, j, 2]
                acc -= T(sj_local * J[3, 3]) * X[chain_idx, me, 1]
            else
                acc += CT(conj(Q) * phase) * X[chain_idx, j, 2]
                acc += CT(P * phase) * X[chain_idx, j, 1]
                acc -= T(sj_local * J[3, 3]) * X[chain_idx, me, 2]
            end
        end

        if !iszero(coupling.biquad)
            K = coupling.biquad
            Sj2Si = sj_local^2 * si_local
            Q = T(0.5) * sij_local^3 * ( K[4, 4] + K[2, 2] - im*(-K[4, 2] + K[2, 4]))
            P = T(0.5) * sij_local^3 * (-K[4, 4] + K[2, 2] - im*( K[4, 2] + K[2, 4]))

            if sigma == 1
                acc += T(-12 * Sj2Si * K[3, 3]) * X[chain_idx, me, 1]
                acc += CT(4 * Sj2Si * (K[1, 3] + im*K[5, 3])) * X[chain_idx, me, 2]
                acc += CT(Q * phase) * X[chain_idx, j, 1]
                acc += CT(conj(P) * phase) * X[chain_idx, j, 2]
            else
                acc += T(-12 * Sj2Si * K[3, 3]) * X[chain_idx, me, 2]
                acc += CT(4 * Sj2Si * (K[1, 3] - im*K[5, 3])) * X[chain_idx, me, 1]
                acc += CT(conj(Q) * phase) * X[chain_idx, j, 2]
                acc += CT(P * phase) * X[chain_idx, j, 1]
            end
        end
    end

    (; incoming_pairs, incoming_indices) = incoming
    for idx in incoming_indices[me]:incoming_indices[me + 1] - 1
        coupling = incoming_pairs[idx]
        bond = coupling.bond
        nbr = bond.i

        phase = exp(T(2π)*im * dot(q, bond.n))

        si_local  = sqrtS[nbr]^2
        sj_local  = sqrtS[me]^2
        sij_local = sqrtS[nbr] * sqrtS[me]

        if !iszero(coupling.bilin)
            J = coupling.bilin
            P = T(0.5) * sij_local * (J[1, 1] - J[2, 2] - im*J[1, 2] - im*J[2, 1])
            Q = T(0.5) * sij_local * (J[1, 1] + J[2, 2] - im*J[1, 2] + im*J[2, 1])

            if sigma == 1
                acc += CT(conj(P) * conj(phase)) * X[chain_idx, nbr, 2]
                acc += CT(conj(Q) * conj(phase)) * X[chain_idx, nbr, 1]
                acc -= T(si_local * J[3, 3]) * X[chain_idx, me, 1]
            else
                acc += CT(Q * conj(phase)) * X[chain_idx, nbr, 2]
                acc += CT(P * conj(phase)) * X[chain_idx, nbr, 1]
                acc -= T(si_local * J[3, 3]) * X[chain_idx, me, 2]
            end
        end

        if !iszero(coupling.biquad)
            K = coupling.biquad
            Si2Sj = si_local^2 * sj_local
            Q = T(0.5) * sij_local^3 * ( K[4, 4] + K[2, 2] - im*(-K[4, 2] + K[2, 4]))
            P = T(0.5) * sij_local^3 * (-K[4, 4] + K[2, 2] - im*( K[4, 2] + K[2, 4]))

            if sigma == 1
                acc += T(-12 * Si2Sj * K[3, 3]) * X[chain_idx, me, 1]
                acc += CT(4 * Si2Sj * (K[3, 1] + im*K[3, 5])) * X[chain_idx, me, 2]
                acc += CT(conj(Q) * conj(phase)) * X[chain_idx, nbr, 1]
                acc += CT(conj(P) * conj(phase)) * X[chain_idx, nbr, 2]
            else
                acc += T(-12 * Si2Sj * K[3, 3]) * X[chain_idx, me, 2]
                acc += CT(4 * Si2Sj * (K[3, 1] - im*K[3, 5])) * X[chain_idx, me, 1]
                acc += CT(Q * conj(phase)) * X[chain_idx, nbr, 2]
                acc += CT(P * conj(phase)) * X[chain_idx, nbr, 1]
            end
        end
    end

    acc += T(regularization) * X[chain_idx, me, sigma]

    Y[chain_idx, me, sigma] = acc
end

function multiply_by_hamiltonian_dipole_batched_ka!(
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    swt::SpinWaveTheoryDeviceKA,
    incoming::IncomingBondCSRKA,
    q_reshaped_chain::AbstractVector,
)
    L = Sunny.natoms(swt.sys.crystal)
    N_chains = size(Y, 1)
    @assert size(Y) == size(X) == (N_chains, 2L) "Bogoliubov vector size mismatch"
    @assert length(q_reshaped_chain) == N_chains "q_reshaped_chain length must match N_chains"
    @assert eltype(Y) == eltype(X) "Y and X must have the same element type"

    CT = eltype(Y)
    T  = real(CT)
    fill!(Y, CT(0))

    Y3 = reshape(Y, (N_chains, L, 2))
    X3 = reshape(X, (N_chains, L, 2))

    total = N_chains * L * 2
    backend = KernelAbstractions.get_backend(Y)
    workgroup = min(total, 256)
    workgroup = workgroup < 1 ? 1 : workgroup
    kernel! = multiply_by_hamiltonian_dipole_batched_kernel_ka!(backend, workgroup)

    kernel!(Y3, X3, swt.sys, incoming, swt.data,
            q_reshaped_chain, T(swt.regularization), L; ndrange=total)

    return Y
end

function mul_dynamical_matrix_batched_ka!(
    swt::SpinWaveTheoryDeviceKA,
    incoming::IncomingBondCSRKA,
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    q_reshaped_chain::AbstractVector,
)
    @assert swt.sys.mode in (dipole, dipole_uncorrected) "GPU batched matvec only supports dipole / dipole_uncorrected mode"
    multiply_by_hamiltonian_dipole_batched_ka!(Y, X, swt, incoming, q_reshaped_chain)
    return Y
end

@kernel inbounds=true function precompute_bond_phases_kernel_ka!(
    phases,
    q_per_iq,
    bond_arr,
    Nq,
)
    tid = @index(Global, Linear)
    T  = @uniform real(eltype(phases))

    iq       = ((tid - 1) % Nq) + 1
    bond_idx = ((tid - 1) ÷ Nq) + 1

    q = q_per_iq[iq]
    n = bond_arr[bond_idx].bond.n
    phases[iq, bond_idx] = exp(T(2π) * im * dot(q, n))
end

function precompute_bond_phases_ka!(
    phases,
    q_per_iq,
    bond_arr,
    Nq::Int,
)
    total_bonds = size(phases, 2)
    (total_bonds == 0 || Nq == 0) && return phases
    total     = Nq * total_bonds
    backend   = KernelAbstractions.get_backend(phases)
    workgroup = min(total, 256)
    workgroup = max(workgroup, 1)
    kernel!   = precompute_bond_phases_kernel_ka!(backend, workgroup)
    kernel!(phases, q_per_iq, bond_arr, Nq; ndrange=total)
    return phases
end

@kernel inbounds=true function multiply_by_hamiltonian_dipole_batched_phased_precomp_kernel_ka!(
    Y, X, indices, incoming_indices,
    isculled_out, bond_j_out,
    bond_i_in,
    onsite_BA1, onsite_A2,
    bilin_P_out, bilin_Q_out, bilin_Jzz_sj_out,
    biquad_P_out, biquad_Q_out, biquad_Kzz_Sj2Si_out, biquad_K13_K53_Sj2Si_out,
    bilin_P_in, bilin_Q_in, bilin_Jzz_si_in,
    biquad_P_in, biquad_Q_in, biquad_Kzz_Si2Sj_in, biquad_K31_K35_Si2Sj_in,
    phases_out, phases_in,
    regularization, L, Nobs,
)
    tid = @index(Global, Linear)
    N_chains = size(Y, 1)

    chain_idx = ((tid - 1) % N_chains) + 1
    rest      =  (tid - 1) ÷ N_chains
    me        =  (rest % L) + 1
    sigma     =  (rest ÷ L) + 1

    iq = ((chain_idx - 1) ÷ Nobs) + 1

    T = real(eltype(Y))
    CT = Complex{T}

    acc = CT(0)

    BA1 = onsite_BA1[me]
    A2  = onsite_A2[me]
    if sigma == 1
        acc += BA1 * X[chain_idx, me, 1]
        acc += A2  * X[chain_idx, me, 2]
    else
        acc += BA1 * X[chain_idx, me, 2]
        acc += conj(A2) * X[chain_idx, me, 1]
    end

    for idx in indices[me]:indices[me + 1] - 1
        isculled_out[idx] && break
        j = Int(bond_j_out[idx])
        phase = phases_out[iq, idx]

        P  = bilin_P_out[idx]
        Q  = bilin_Q_out[idx]
        Jz = bilin_Jzz_sj_out[idx]

        bqP  = biquad_P_out[idx]
        bqQ  = biquad_Q_out[idx]
        bqKd = biquad_Kzz_Sj2Si_out[idx]
        bqKc = biquad_K13_K53_Sj2Si_out[idx]

        if sigma == 1
            acc += Q * phase * X[chain_idx, j, 1]
            acc += conj(P) * phase * X[chain_idx, j, 2]
            acc -= Jz * X[chain_idx, me, 1]

            acc += bqKd * X[chain_idx, me, 1]
            acc += bqKc * X[chain_idx, me, 2]
            acc += bqQ * phase * X[chain_idx, j, 1]
            acc += conj(bqP) * phase * X[chain_idx, j, 2]
        else
            acc += conj(Q) * phase * X[chain_idx, j, 2]
            acc += P * phase * X[chain_idx, j, 1]
            acc -= Jz * X[chain_idx, me, 2]

            acc += bqKd * X[chain_idx, me, 2]
            acc += conj(bqKc) * X[chain_idx, me, 1]
            acc += conj(bqQ) * phase * X[chain_idx, j, 2]
            acc += bqP * phase * X[chain_idx, j, 1]
        end
    end

    for idx in incoming_indices[me]:incoming_indices[me + 1] - 1
        nbr = Int(bond_i_in[idx])
        phase = phases_in[iq, idx]

        P  = bilin_P_in[idx]
        Q  = bilin_Q_in[idx]
        Jz = bilin_Jzz_si_in[idx]
        bqP  = biquad_P_in[idx]
        bqQ  = biquad_Q_in[idx]
        bqKd = biquad_Kzz_Si2Sj_in[idx]
        bqKc = biquad_K31_K35_Si2Sj_in[idx]

        if sigma == 1
            acc += conj(P) * conj(phase) * X[chain_idx, nbr, 2]
            acc += conj(Q) * conj(phase) * X[chain_idx, nbr, 1]
            acc -= Jz * X[chain_idx, me, 1]

            acc += bqKd * X[chain_idx, me, 1]
            acc += bqKc * X[chain_idx, me, 2]
            acc += conj(bqQ) * conj(phase) * X[chain_idx, nbr, 1]
            acc += conj(bqP) * conj(phase) * X[chain_idx, nbr, 2]
        else
            acc += Q * conj(phase) * X[chain_idx, nbr, 2]
            acc += P * conj(phase) * X[chain_idx, nbr, 1]
            acc -= Jz * X[chain_idx, me, 2]

            acc += bqKd * X[chain_idx, me, 2]
            acc += conj(bqKc) * X[chain_idx, me, 1]
            acc += bqQ * conj(phase) * X[chain_idx, nbr, 2]
            acc += bqP * conj(phase) * X[chain_idx, nbr, 1]
        end
    end

    acc += T(regularization) * X[chain_idx, me, sigma]

    Y[chain_idx, me, sigma] = acc
end

function multiply_by_hamiltonian_dipole_batched_phased_precomp_ka!(
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    swt::SpinWaveTheoryDeviceKA,
    incoming::IncomingBondCSRKA,
    precomp,
    Nobs::Int,
)
    L = Sunny.natoms(swt.sys.crystal)
    N_chains = size(Y, 1)
    @assert size(Y) == size(X) == (N_chains, 2L) "Bogoliubov vector size mismatch"
    @assert eltype(Y) == eltype(X) "Y and X must have the same element type"

    CT = eltype(Y)
    T  = real(CT)
    fill!(Y, CT(0))

    Y3 = reshape(Y, (N_chains, L, 2))
    X3 = reshape(X, (N_chains, L, 2))

    total     = N_chains * L * 2
    backend   = KernelAbstractions.get_backend(Y)
    workgroup = max(min(total, 256), 1)
    kernel!   = multiply_by_hamiltonian_dipole_batched_phased_precomp_kernel_ka!(backend, workgroup)
    kernel!(Y3, X3, swt.sys.indices, incoming.incoming_indices,
            precomp.isculled_out, precomp.bond_j_out,
            precomp.bond_i_in,
            precomp.onsite_BA1, precomp.onsite_A2,
            precomp.bilin_P_out, precomp.bilin_Q_out, precomp.bilin_Jzz_sj_out,
            precomp.biquad_P_out, precomp.biquad_Q_out,
            precomp.biquad_Kzz_Sj2Si_out, precomp.biquad_K13_K53_Sj2Si_out,
            precomp.bilin_P_in, precomp.bilin_Q_in, precomp.bilin_Jzz_si_in,
            precomp.biquad_P_in, precomp.biquad_Q_in,
            precomp.biquad_Kzz_Si2Sj_in, precomp.biquad_K31_K35_Si2Sj_in,
            precomp.phases_out, precomp.phases_in,
            T(swt.regularization), L, Nobs;
            ndrange=total)
    return Y
end

function mul_dynamical_matrix_batched_phased_precomp_ka!(
    swt::SpinWaveTheoryDeviceKA,
    incoming::IncomingBondCSRKA,
    Y::AbstractMatrix{<:Complex},
    X::AbstractMatrix{<:Complex},
    precomp,
    Nobs::Int,
)
    @assert swt.sys.mode in (dipole, dipole_uncorrected) "GPU batched matvec only supports dipole / dipole_uncorrected mode"
    multiply_by_hamiltonian_dipole_batched_phased_precomp_ka!(Y, X, swt, incoming, precomp, Nobs)
    return Y
end
