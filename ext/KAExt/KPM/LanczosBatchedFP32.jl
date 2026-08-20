@kernel inbounds=true function batched_dot_chains_fp32_kernel_ka!(result, a, b, twoL)
    chain = @index(Group, Linear)
    lid   = @index(Local, Linear)
    groupsize = @uniform @groupsize()[1]

    lmem_real = @localmem Float32 (BATCHED_DOT_WG,)
    lmem_imag = @localmem Float32 (BATCHED_DOT_WG,)

    acc = ComplexF32(0)
    k = lid
    while k <= twoL
        acc += conj(ComplexF32(a[chain, k])) * ComplexF32(b[chain, k])
        k += groupsize
    end

    lmem_real[lid] = real(acc)
    lmem_imag[lid] = imag(acc)
    @synchronize

    if lid <= 32
        lmem_real[lid] += lmem_real[lid + 32]
        lmem_imag[lid] += lmem_imag[lid + 32]
    end
    @synchronize
    if lid <= 16
        lmem_real[lid] += lmem_real[lid + 16]
        lmem_imag[lid] += lmem_imag[lid + 16]
    end
    @synchronize
    if lid <= 8
        lmem_real[lid] += lmem_real[lid + 8]
        lmem_imag[lid] += lmem_imag[lid + 8]
    end
    @synchronize
    if lid <= 4
        lmem_real[lid] += lmem_real[lid + 4]
        lmem_imag[lid] += lmem_imag[lid + 4]
    end
    @synchronize
    if lid <= 2
        lmem_real[lid] += lmem_real[lid + 2]
        lmem_imag[lid] += lmem_imag[lid + 2]
    end
    @synchronize
    if lid == 1
        result[chain] = ComplexF32(lmem_real[1] + lmem_real[2], lmem_imag[1] + lmem_imag[2])
    end
end

function batched_dot_chains_fp32_ka!(
    result::AbstractVector{ComplexF32},
    a::AbstractMatrix{<:Complex},
    b::AbstractMatrix{<:Complex},
)
    N_chains = size(a, 1)
    twoL = size(a, 2)
    @assert size(a) == size(b) "batched_dot_chains_fp32_ka! a and b must have the same shape"
    @assert length(result) == N_chains "batched_dot_chains_fp32_ka! result length must match N_chains"
    @assert eltype(result) == ComplexF32 "batched_dot_chains_fp32_ka! result must be ComplexF32"

    backend = KernelAbstractions.get_backend(a)
    kernel! = batched_dot_chains_fp32_kernel_ka!(backend, BATCHED_DOT_WG)
    kernel!(result, a, b, twoL; ndrange = N_chains * BATCHED_DOT_WG)
    return result
end

function build_fp32_hv_precomp(swt::Sunny.SpinWaveTheory, backend,
                               L::Int, total_out::Int, total_in::Int,
                               phases_out_dev, phases_in_dev)
    host_sys  = swt.sys
    host_data = swt.data::Sunny.SWTDataDipole

    sqrtS_h     = host_data.sqrtS
    sc_h        = host_data.stevens_coefs
    lr_h        = host_data.local_rotations
    gs_h        = host_sys.gs
    ext_h       = host_sys.extfield
    interactions_h = host_sys.interactions_union

    @assert length(interactions_h) == L "interactions_union length $(length(interactions_h)) != L=$L"

    onsite_BA1_host = Vector{Float32}(undef, L)
    onsite_A2_host  = Vector{ComplexF32}(undef, L)
    @inbounds for me in 1:L
        s_me = Float64(sqrtS_h[me])^2
        c2 = sc_h[me].c2
        c4 = sc_h[me].c4
        c6 = sc_h[me].c6
        A1 = -6 * s_me * Float64(c2[3]) - 80 * s_me^3 * Float64(c4[5]) - 336 * s_me^5 * Float64(c6[7])
        A2 =  2 * s_me * (Float64(c2[1]) + im*Float64(c2[5])) +
             12 * s_me^3 * (Float64(c4[3]) + im*Float64(c4[7])) +
             32 * s_me^5 * (Float64(c6[5]) + im*Float64(c6[9]))
        B_f64      = Vector{Float64}(gs_h[1,1,1,me]' * ext_h[1,1,1,me])
        lr3_f64    = Vector{Float64}(lr_h[me][:, 3])
        Bprime     = -dot(B_f64, lr3_f64)
        onsite_BA1_host[me] = Float32(Bprime + A1)
        onsite_A2_host[me]  = ComplexF32(A2)
    end

    @inline get_bilin(J, i, j) = J isa Real ? Float64(J) : Float64(J[i, j])
    @inline get_biquad(K, i, j) = K isa Real ? Float64(K) : Float64(K[i, j])

    isculled_out_host              = Vector{Bool}(undef, total_out)
    bond_j_out_host                = Vector{Int32}(undef, total_out)
    bilin_P_out_host               = Vector{ComplexF32}(undef, total_out)
    bilin_Q_out_host               = Vector{ComplexF32}(undef, total_out)
    bilin_Jzz_sj_out_host          = Vector{Float32}(undef, total_out)
    biquad_P_out_host              = Vector{ComplexF32}(undef, total_out)
    biquad_Q_out_host              = Vector{ComplexF32}(undef, total_out)
    biquad_Kzz_Sj2Si_out_host      = Vector{Float32}(undef, total_out)
    biquad_K13_K53_Sj2Si_out_host  = Vector{ComplexF32}(undef, total_out)

    idx_flat = 0
    @inbounds for me in 1:L
        sqrtS_me_f64 = Float64(sqrtS_h[me])
        for pair in interactions_h[me].pair
            idx_flat += 1
            j = pair.bond.j
            sqrtS_j_f64 = Float64(sqrtS_h[j])
            sij_f64 = sqrtS_me_f64 * sqrtS_j_f64
            si_me   = sqrtS_me_f64^2
            sj      = sqrtS_j_f64^2

            J = pair.bilin
            J11 = get_bilin(J, 1, 1); J22 = get_bilin(J, 2, 2)
            J12 = get_bilin(J, 1, 2); J21 = get_bilin(J, 2, 1)
            J33 = get_bilin(J, 3, 3)
            P = 0.5 * sij_f64 * (J11 - J22 - im*J12 - im*J21)
            Q = 0.5 * sij_f64 * (J11 + J22 - im*J12 + im*J21)

            K = pair.biquad
            K44 = get_biquad(K, 4, 4); K22 = get_biquad(K, 2, 2)
            K42 = get_biquad(K, 4, 2); K24 = get_biquad(K, 2, 4)
            K33 = get_biquad(K, 3, 3); K13 = get_biquad(K, 1, 3); K53 = get_biquad(K, 5, 3)
            sij3 = sij_f64^3
            bqQ    = 0.5 * sij3 * (K44 + K22 - im*(-K42 + K24))
            bqP    = 0.5 * sij3 * (-K44 + K22 - im*(K42 + K24))
            Sj2Si  = sj^2 * si_me
            bqKzz  = -12 * Sj2Si * K33
            bqK13  =   4 * Sj2Si * (K13 + im*K53)

            isculled_out_host[idx_flat]              = pair.isculled
            bond_j_out_host[idx_flat]                = Int32(j)
            bilin_P_out_host[idx_flat]               = ComplexF32(P)
            bilin_Q_out_host[idx_flat]               = ComplexF32(Q)
            bilin_Jzz_sj_out_host[idx_flat]          = Float32(sj * J33)
            biquad_P_out_host[idx_flat]              = ComplexF32(bqP)
            biquad_Q_out_host[idx_flat]              = ComplexF32(bqQ)
            biquad_Kzz_Sj2Si_out_host[idx_flat]      = Float32(bqKzz)
            biquad_K13_K53_Sj2Si_out_host[idx_flat]  = ComplexF32(bqK13)
        end
    end
    @assert idx_flat == total_out "outgoing bond flat count mismatch: $idx_flat vs $total_out"

    incoming_counts = zeros(Int, L)
    for int in interactions_h
        for pair in int.pair
            pair.isculled && continue
            incoming_counts[pair.bond.j] += 1
        end
    end
    incoming_idx_offsets = Vector{Int}(undef, L + 1)
    incoming_idx_offsets[1] = 1
    for jj in 1:L
        incoming_idx_offsets[jj+1] = incoming_idx_offsets[jj] + incoming_counts[jj]
    end
    @assert incoming_idx_offsets[L+1] - 1 == total_in "incoming total mismatch: $(incoming_idx_offsets[L+1] - 1) vs $total_in"

    bond_i_in_host                 = Vector{Int32}(undef, total_in)
    bilin_P_in_host                = Vector{ComplexF32}(undef, total_in)
    bilin_Q_in_host                = Vector{ComplexF32}(undef, total_in)
    bilin_Jzz_si_in_host           = Vector{Float32}(undef, total_in)
    biquad_P_in_host               = Vector{ComplexF32}(undef, total_in)
    biquad_Q_in_host               = Vector{ComplexF32}(undef, total_in)
    biquad_Kzz_Si2Sj_in_host       = Vector{Float32}(undef, total_in)
    biquad_K31_K35_Si2Sj_in_host   = Vector{ComplexF32}(undef, total_in)

    write_pos = copy(@view incoming_idx_offsets[1:L])
    @inbounds for me in 1:L
        sqrtS_me_f64 = Float64(sqrtS_h[me])
        for pair in interactions_h[me].pair
            pair.isculled && continue
            j = pair.bond.j
            sqrtS_j_f64 = Float64(sqrtS_h[j])
            sij_f64 = sqrtS_me_f64 * sqrtS_j_f64

            si_sender   = sqrtS_me_f64^2
            sj_receiver = sqrtS_j_f64^2

            J = pair.bilin
            J11 = get_bilin(J, 1, 1); J22 = get_bilin(J, 2, 2)
            J12 = get_bilin(J, 1, 2); J21 = get_bilin(J, 2, 1)
            J33 = get_bilin(J, 3, 3)
            P = 0.5 * sij_f64 * (J11 - J22 - im*J12 - im*J21)
            Q = 0.5 * sij_f64 * (J11 + J22 - im*J12 + im*J21)

            K = pair.biquad
            K44 = get_biquad(K, 4, 4); K22 = get_biquad(K, 2, 2)
            K42 = get_biquad(K, 4, 2); K24 = get_biquad(K, 2, 4)
            K33 = get_biquad(K, 3, 3); K31 = get_biquad(K, 3, 1); K35 = get_biquad(K, 3, 5)
            sij3 = sij_f64^3
            bqQ    = 0.5 * sij3 * (K44 + K22 - im*(-K42 + K24))
            bqP    = 0.5 * sij3 * (-K44 + K22 - im*(K42 + K24))
            Si2Sj  = si_sender^2 * sj_receiver
            bqKzz  = -12 * Si2Sj * K33
            bqK31  =   4 * Si2Sj * (K31 + im*K35)

            pos = write_pos[j]
            bond_i_in_host[pos]                 = Int32(me)
            bilin_P_in_host[pos]                = ComplexF32(P)
            bilin_Q_in_host[pos]                = ComplexF32(Q)
            bilin_Jzz_si_in_host[pos]           = Float32(si_sender * J33)
            biquad_P_in_host[pos]               = ComplexF32(bqP)
            biquad_Q_in_host[pos]               = ComplexF32(bqQ)
            biquad_Kzz_Si2Sj_in_host[pos]       = Float32(bqKzz)
            biquad_K31_K35_Si2Sj_in_host[pos]   = ComplexF32(bqK31)
            write_pos[j] += 1
        end
    end
    @assert all(write_pos[jj] == incoming_idx_offsets[jj+1] for jj in 1:L) "incoming write_pos sanity check failed"

    return (
        isculled_out             = _backend_array(backend, isculled_out_host),
        bond_j_out               = _backend_array(backend, bond_j_out_host),
        bond_i_in                = _backend_array(backend, bond_i_in_host),
        onsite_BA1               = _backend_array(backend, onsite_BA1_host),
        onsite_A2                = _backend_array(backend, onsite_A2_host),
        bilin_P_out              = _backend_array(backend, bilin_P_out_host),
        bilin_Q_out              = _backend_array(backend, bilin_Q_out_host),
        bilin_Jzz_sj_out         = _backend_array(backend, bilin_Jzz_sj_out_host),
        biquad_P_out             = _backend_array(backend, biquad_P_out_host),
        biquad_Q_out             = _backend_array(backend, biquad_Q_out_host),
        biquad_Kzz_Sj2Si_out     = _backend_array(backend, biquad_Kzz_Sj2Si_out_host),
        biquad_K13_K53_Sj2Si_out = _backend_array(backend, biquad_K13_K53_Sj2Si_out_host),
        bilin_P_in               = _backend_array(backend, bilin_P_in_host),
        bilin_Q_in               = _backend_array(backend, bilin_Q_in_host),
        bilin_Jzz_si_in          = _backend_array(backend, bilin_Jzz_si_in_host),
        biquad_P_in              = _backend_array(backend, biquad_P_in_host),
        biquad_Q_in              = _backend_array(backend, biquad_Q_in_host),
        biquad_Kzz_Si2Sj_in      = _backend_array(backend, biquad_Kzz_Si2Sj_in_host),
        biquad_K31_K35_Si2Sj_in  = _backend_array(backend, biquad_K31_K35_Si2Sj_in_host),
        phases_out               = phases_out_dev,
        phases_in                = phases_in_dev,
    )
end

@kernel inbounds=true function accumulate_lhs_adj_Q_kernel_ka!(
    lhs_adj_Q_dev, lhs, v_batch, Nobs, twoL, iter,
)
    tid = @index(Global, Linear)
    CT = @uniform eltype(lhs_adj_Q_dev)

    μ = ((tid - 1) % Nobs) + 1
    c = ((tid - 1) ÷ Nobs) + 1

    iq = ((c - 1) ÷ Nobs) + 1

    acc = CT(0)
    for k in 1:twoL
        acc += conj(CT(lhs[k, μ, iq])) * CT(v_batch[c, k])
    end
    lhs_adj_Q_dev[μ, iter, c] = acc
end

function accumulate_lhs_adj_Q_ka!(
    lhs_adj_Q_dev::AbstractArray{<:Complex, 3},
    lhs::AbstractArray{<:Complex, 3},
    v_batch::AbstractMatrix{<:Complex},
    Nobs::Int, iter::Int,
)
    N_chains = size(v_batch, 1)
    twoL = size(v_batch, 2)
    total = Nobs * N_chains
    backend = KernelAbstractions.get_backend(lhs_adj_Q_dev)
    workgroup = max(min(total, 256), 1)
    kernel! = accumulate_lhs_adj_Q_kernel_ka!(backend, workgroup)
    kernel!(lhs_adj_Q_dev, lhs, v_batch, Nobs, twoL, iter; ndrange=total)
    return lhs_adj_Q_dev
end

@kernel inbounds=true function fused_vp_svp_rescale_kernel_ka!(vp, svp, w, sw, β_chain)
    tid = @index(Global, Linear)
    N_chains = size(vp, 1)
    CT = @uniform eltype(vp)

    c = ((tid - 1) % N_chains) + 1
    k = ((tid - 1) ÷ N_chains) + 1

    inv_β = one(eltype(β_chain)) / β_chain[c]
    vp[c, k]  = CT(w[c, k]  * inv_β)
    svp[c, k] = CT(sw[c, k] * inv_β)
end

function fused_vp_svp_rescale_ka!(
    vp::AbstractMatrix{<:Complex},
    svp::AbstractMatrix{<:Complex},
    w::AbstractMatrix{<:Complex},
    sw::AbstractMatrix{<:Complex},
    β_chain::AbstractVector{<:Real},
)
    total = length(vp)
    @assert size(vp) == size(svp) == size(w) == size(sw) "vp/svp/w/sw must all be (N_chains, 2L)"
    backend = KernelAbstractions.get_backend(vp)
    workgroup = max(min(total, 256), 1)
    kernel! = fused_vp_svp_rescale_kernel_ka!(backend, workgroup)
    kernel!(vp, svp, w, sw, β_chain; ndrange=total)
    return vp
end

function lanczos_device_batched_fp32_ka(
    mulA!,
    mulS!,
    v_batch::AbstractMatrix{<:Complex},
    vp_batch::AbstractMatrix{<:Complex},
    Sv_batch::AbstractMatrix{<:Complex},
    Svp_batch::AbstractMatrix{<:Complex},
    w_batch::AbstractMatrix{<:Complex},
    Sw_batch::AbstractMatrix{<:Complex},
    lhs_dev::AbstractArray{<:Complex, 3},
    lhs_adj_Q_dev::AbstractArray{ComplexF32, 3},
    α_chains_dev::AbstractVector{ComplexF32},
    β²_chains_dev::AbstractVector{ComplexF32},
    β_chains_dev::AbstractVector{<:Real},
    α_chains_host::Vector{ComplexF32},
    β²_chains_host::Vector{ComplexF32};
    max_iters::Int,
    min_iters::Int,
    resolution::Float64,
    verbose::Bool = false,
)::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Int}, Vector{Float64}}

    N_chains = size(v_batch, 1)
    twoL = size(v_batch, 2)
    Nobs = size(lhs_dev, 2)
    Nq   = size(lhs_dev, 3)

    @assert size(v_batch) == size(vp_batch) == size(Sv_batch) == size(Svp_batch) == size(w_batch) == size(Sw_batch) "Recurrence buffers must all be (N_chains, 2L)"
    @assert size(lhs_dev, 1) == twoL "lhs_dev first axis must equal 2L"
    @assert Nobs * Nq == N_chains "lhs_dev's Nobs*Nq must equal N_chains"
    @assert size(lhs_adj_Q_dev) == (Nobs, max_iters, N_chains) "lhs_adj_Q_dev must be (Nobs, max_iters, N_chains)"
    @assert length(α_chains_dev) == N_chains
    @assert length(β²_chains_dev) == N_chains
    @assert length(β_chains_dev) == N_chains
    @assert length(α_chains_host) == N_chains
    @assert length(β²_chains_host) == N_chains
    @assert max_iters >= 1 "max_iters must be >= 1"

    αs_matrix = Matrix{Float64}(undef, max_iters, N_chains)
    βs_matrix = Matrix{Float64}(undef, max_iters - 1, N_chains)
    n_iters_done_per_chain = zeros(Int, N_chains)
    chain_done = falses(N_chains)
    c_norm_per_chain = Vector{Float64}(undef, N_chains)
    niters_eff = fill(max_iters, N_chains)

    mulS!(Sv_batch, v_batch)

    batched_dot_chains_fp32_ka!(α_chains_dev, v_batch, Sv_batch)

    copyto!(α_chains_host, α_chains_dev)

    T_real = real(eltype(v_batch))
    hermitian_atol = max(twoL, 1) * max(1f-5, 100 * eps(T_real))
    near_zero_atol = max(twoL, 1) * max(1f-6, 100 * eps(T_real))
    patched_alpha = false
    @inbounds for c in 1:N_chains
        vSv_c = α_chains_host[c]
        if abs(imag(vSv_c)) > hermitian_atol
            error("Chain $c: S not Hermitian (imag(v†Sv) = $(imag(vSv_c)), |imag| = $(abs(imag(vSv_c))) > $hermitian_atol)")
        end
        r = Float64(real(vSv_c))
        if r < -near_zero_atol
            error("Chain $c: v†Sv strongly negative ($vSv_c); S not positive definite")
        elseif r <= near_zero_atol
            chain_done[c] = true
            c_norm_per_chain[c] = 0.0
            α_chains_host[c] = one(eltype(α_chains_host))
            patched_alpha = true
        else
            c_norm_per_chain[c] = sqrt(r)
        end
    end
    if patched_alpha
        copyto!(α_chains_dev, α_chains_host)
    end

    v_batch  ./= sqrt.(real.(reshape(α_chains_dev, N_chains, 1)))
    Sv_batch ./= sqrt.(real.(reshape(α_chains_dev, N_chains, 1)))

    mulA!(w_batch, Sv_batch)

    batched_dot_chains_fp32_ka!(α_chains_dev, w_batch, Sv_batch)

    copyto!(α_chains_host, α_chains_dev)

    @inbounds for c in 1:N_chains
        αs_matrix[1, c] = chain_done[c] ? 0.0 : Float64(real(α_chains_host[c]))
    end

    w_batch .-= reshape(α_chains_dev, N_chains, 1) .* v_batch

    mulS!(Sw_batch, w_batch)

    accumulate_lhs_adj_Q_ka!(lhs_adj_Q_dev, lhs_dev, v_batch, Nobs, 1)

    @inbounds for c in 1:N_chains
        n_iters_done_per_chain[c] = 1
    end

    @inbounds for i in 1:max_iters - 1
        batched_dot_chains_fp32_ka!(β²_chains_dev, w_batch, Sw_batch)

        copyto!(β²_chains_host, β²_chains_dev)

        for c in 1:N_chains
            if !chain_done[c]
                β²_real = real(β²_chains_host[c])
                if !(β²_real > 0.0f0)
                    chain_done[c] = true
                end
            end
        end

        if i == min_iters && isfinite(resolution)
            for c in 1:N_chains
                if !chain_done[c]
                    n_so_far = n_iters_done_per_chain[c]
                    if n_so_far >= 2
                        T_partial = SymTridiagonal(
                            αs_matrix[1:n_so_far, c],
                            βs_matrix[1:n_so_far-1, c],
                        )
                        Δϵ_c = eigmax(T_partial) - eigmin(T_partial)
                        if isfinite(Δϵ_c) && Δϵ_c < Float64(max_iters) * resolution
                            niters_c = max(min_iters, fld(Δϵ_c, resolution))
                            niters_c += mod(niters_c, 2)
                            niters_eff[c] = min(Int(niters_c), max_iters)
                        end
                    end
                end
            end
            if verbose
                println("Per-chain niters_eff at min_iters=$min_iters: " *
                        "min=$(minimum(niters_eff)), max=$(maximum(niters_eff)), " *
                        "median=$(sort(niters_eff)[N_chains ÷ 2])")
            end
        end

        for c in 1:N_chains
            if !chain_done[c] && n_iters_done_per_chain[c] >= niters_eff[c]
                chain_done[c] = true
            end
        end

        all(chain_done) && break

        β_chains_dev .= ifelse.(real.(β²_chains_dev) .> zero(T_real),
                                sqrt.(max.(real.(β²_chains_dev), zero(T_real))),
                                one(T_real))

        fused_vp_svp_rescale_ka!(vp_batch, Svp_batch, w_batch, Sw_batch, β_chains_dev)

        mulA!(w_batch, Svp_batch)

        batched_dot_chains_fp32_ka!(α_chains_dev, w_batch, Svp_batch)

        copyto!(α_chains_host, α_chains_dev)

        @inbounds for c in 1:N_chains
            if !chain_done[c]
                αs_matrix[i + 1, c] = Float64(real(α_chains_host[c]))
                βs_matrix[i, c] = sqrt(Float64(real(β²_chains_host[c])))
                n_iters_done_per_chain[c] = i + 1
            end
        end

        batched_w_update_ka!(w_batch, α_chains_dev, vp_batch, β_chains_dev, v_batch)

        mulS!(Sw_batch, w_batch)

        v_batch, vp_batch = vp_batch, v_batch

        accumulate_lhs_adj_Q_ka!(lhs_adj_Q_dev, lhs_dev, v_batch, Nobs, i + 1)

        if verbose
            n_done = count(chain_done)
            println("lanczos_device_batched_fp32_ka iter $i: $n_done/$N_chains chains done")
        end
    end

    return αs_matrix, βs_matrix, n_iters_done_per_chain, c_norm_per_chain
end

function intensities_lanczos_device_batched_fp32_ka!(data, swt_kry_d::SpinWaveTheoryKPMBatchedDeviceKA{T}, qpts;
                                                     energies, kernel, kT, verbose) where T
    swt = swt_kry_d.swt_host
    sys = swt.sys
    measure = swt.measure
    cryst = Sunny.orig_crystal(sys)

    isnothing(kernel.fwhm) && error("Cannot determine the kernel fwhm")

    @assert eltype(data) == eltype(measure)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    fill!(data, zero(eltype(data)))

    Na      = Sunny.nsites(sys)
    Ncells  = Na / Sunny.natoms(cryst)
    Nf      = Sunny.nflavors(swt)
    L       = Nf * Na
    twoL    = 2 * L

    Nobs    = size(measure.observables, 1)
    Ncorr   = length(measure.corr_pairs)
    Nq      = length(qpts.qs)
    N_chains = Nq * Nobs

    tol           = swt_kry_d.tol
    niters_user   = swt_kry_d.niters
    niters_bounds = swt_kry_d.niters_bounds

    CT = Complex{T}

    swt_d_sys = swt_kry_d.swt_d.sys
    backend = if swt_d_sys isa SystemDeviceSUNF32KA
        KernelAbstractions.get_backend(swt_d_sys.W_onsite)
    else
        KernelAbstractions.get_backend(swt_d_sys.extfield)
    end

    resolution = Inf
    if niters_user > 0
        @assert tol == 1

        max_iters_global = niters_user
        max_iters_global += mod(max_iters_global, 2)
        max_iters_global = min(max_iters_global, twoL)
    else
        @assert 0.0 < tol <= 1
        q_idx_rep = max(1, length(qpts.qs) ÷ 2)
        q_rep = qpts.qs[q_idx_rep]
        q_rep_reshaped = Sunny.to_reshaped_rlu(sys, q_rep)
        lo, hi = Sunny.eigbounds(swt, q_rep_reshaped, niters_bounds)
        Δϵ = hi - lo
        resolution = (kernel.fwhm/2) / (-log10(tol))
        max_iters_unbounded = max(niters_bounds, fld(Δϵ, resolution))
        max_iters_with_safety = ceil(Int, 4.0 * max_iters_unbounded)
        max_iters_with_safety += mod(max_iters_with_safety, 2)
        max_iters_global = min(max_iters_with_safety, twoL)
        if verbose
            println("Eigbounds at q_idx=$q_idx_rep: Δϵ=$Δϵ")
            println("max_iters_global=$max_iters_global (4× safety on single-q estimate, buffer upper bound only)")
            println("Per-chain niters_eff determined adaptively at iteration $niters_bounds")
        end
    end
    @assert max_iters_global >= 1 "max_iters_global must be ≥ 1"

    q_reshaped_host = Sunny.Vec3[
        Sunny.to_reshaped_rlu(sys, q) for q in vec(qpts.qs)
    ]

    v_batch   = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    vp_batch  = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    Sv_batch  = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    Svp_batch = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    w_batch   = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    Sw_batch  = KernelAbstractions.zeros(backend, CT, N_chains, twoL)

    lhs_dev       = KernelAbstractions.zeros(backend, CT,         twoL, Nobs, Nq)
    lhs_adj_Q_dev = KernelAbstractions.zeros(backend, ComplexF32, Nobs, max_iters_global, N_chains)

    α_chains_dev   = KernelAbstractions.zeros(backend, ComplexF32, N_chains)
    β²_chains_dev  = KernelAbstractions.zeros(backend, ComplexF32, N_chains)
    β_chains_dev   = KernelAbstractions.zeros(backend, T,          N_chains)
    α_chains_host  = Vector{ComplexF32}(undef, N_chains)
    β²_chains_host = Vector{ComplexF32}(undef, N_chains)

    u_host_per_chain = zeros(ComplexF64, twoL, N_chains)
    Avec_pref = zeros(ComplexF64, Na)

    if sys.mode == :SUN

        let
            (; observables_localized) = swt.data::Sunny.SWTDataSUN
            N_hilbert = sys.Ns[1]
            for iq in 1:Nq
                q        = qpts.qs[iq]
                q_r      = q_reshaped_host[iq]
                q_global = cryst.recipvecs * q
                for i in 1:Na
                    r  = sys.crystal.positions[i]
                    ff = Sunny.get_swt_formfactor(measure, 1, i)
                    Avec_pref[i]  = exp(2π*im * dot(q_r, r))
                    Avec_pref[i] *= Sunny.compute_form_factor(ff, Sunny.norm2(q_global))
                end
                for ξ in 1:Nobs
                    c = (iq - 1) * Nobs + ξ
                    for i in 1:Na
                        O = observables_localized[ξ, i]
                        for f in 1:Nf
                            l = f + (i - 1) * Nf
                            u_host_per_chain[l,   c] = Avec_pref[i] * O[f, N_hilbert]
                            u_host_per_chain[l+L, c] = Avec_pref[i] * O[N_hilbert, f]
                        end
                    end
                end
            end
        end
    else

        let
            (; sqrtS, observables_localized) = swt.data::Sunny.SWTDataDipole
            for iq in 1:Nq
                q = qpts.qs[iq]
                q_reshaped = q_reshaped_host[iq]
                q_global = cryst.recipvecs * q

                for i in 1:Na
                    r = sys.crystal.positions[i]
                    ff = Sunny.get_swt_formfactor(measure, 1, i)
                    Avec_pref[i] = exp(2π*im * dot(q_reshaped, r))
                    Avec_pref[i] *= Sunny.compute_form_factor(ff, Sunny.norm2(q_global))
                end

                for ξ in 1:Nobs
                    c = (iq - 1) * Nobs + ξ
                    for i in 1:Na
                        O = observables_localized[ξ, i]
                        u_host_per_chain[i,   c] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] + im*O[2])
                        u_host_per_chain[i+L, c] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] - im*O[2])
                    end
                end
            end
        end
    end

    chain_skip = falses(N_chains)
    @inbounds for c in 1:N_chains
        if iszero(view(u_host_per_chain, :, c))
            chain_skip[c] = true
        end
    end

    copyto!(reshape(lhs_dev, twoL, N_chains), convert(Array{CT}, u_host_per_chain))

    u_host_dummy = copy(u_host_per_chain)
    @inbounds for c in 1:N_chains
        if chain_skip[c]
            u_host_dummy[1, c] = ComplexF64(1, 0)
        end
    end

    v_batch_host = permutedims(u_host_dummy, (2, 1))
    @views v_batch_host[:, L+1:twoL] .*= -1
    copyto!(v_batch, convert(Array{CT}, v_batch_host))

    mulA! = let L = L, twoL = twoL
        (out, in_) -> begin
            @views @. out[:, 1:L]      = +in_[:, 1:L]
            @views @. out[:, L+1:twoL] = -in_[:, L+1:twoL]
            return out
        end
    end

    mulS! = if sys.mode == :SUN

        q_chain_sun_dev = let
            q_chain_h = Vector{SVector{3, T}}(undef, N_chains)
            for iq in 1:Nq
                for ξ in 1:Nobs
                    c = (iq - 1) * Nobs + ξ
                    q_chain_h[c] = SVector{3, T}(q_reshaped_host[iq])
                end
            end
            _backend_array(backend, q_chain_h)
        end
        let swt_d = swt_kry_d.swt_d,
            incoming = swt_kry_d.incoming,
            q_chain = q_chain_sun_dev
            (out, in_) -> begin
                mul_dynamical_matrix_batched_ka!(swt_d, incoming, out, in_, q_chain)
                return out
            end
        end
    else

        total_out      = length(swt_kry_d.swt_d.sys.pairs)
        total_in       = length(swt_kry_d.incoming.incoming_pairs)
        q_per_iq_dev   = _backend_array(backend,
                             SVector{3,T}[SVector{3,T}(q_reshaped_host[iq]) for iq in 1:Nq])
        phases_out_dev = KernelAbstractions.zeros(backend, CT, Nq, total_out)
        phases_in_dev  = KernelAbstractions.zeros(backend, CT, Nq, total_in)
        precompute_bond_phases_ka!(phases_out_dev, q_per_iq_dev,
                                   swt_kry_d.swt_d.sys.pairs, Nq)
        precompute_bond_phases_ka!(phases_in_dev,  q_per_iq_dev,
                                   swt_kry_d.incoming.incoming_pairs, Nq)
        precomp = build_fp32_hv_precomp(swt, backend, L, total_out, total_in,
                                        phases_out_dev, phases_in_dev)
        let swt_d = swt_kry_d.swt_d,
            incoming = swt_kry_d.incoming,
            precomp_ = precomp,
            Nobs_ = Nobs
            (out, in_) -> begin
                mul_dynamical_matrix_batched_phased_precomp_ka!(swt_d, incoming, out, in_,
                                                                precomp_, Nobs_)
                return out
            end
        end
    end

    αs_matrix, βs_matrix, n_iters_done_per_chain, c_norm_per_chain = try
        lanczos_device_batched_fp32_ka(
            mulA!, mulS!,
            v_batch, vp_batch, Sv_batch, Svp_batch, w_batch, Sw_batch,
            lhs_dev, lhs_adj_Q_dev,
            α_chains_dev, β²_chains_dev, β_chains_dev,
            α_chains_host, β²_chains_host;
            max_iters = max_iters_global,
            min_iters = niters_bounds,
            resolution = resolution,
            verbose = verbose,
        )
    catch e
        if e isa ErrorException && occursin("not a positive definite measure", e.msg)
            rethrow(ErrorException("GPU batched FP32 Lanczos: not an energy-minimum; some q wavevector unstable. " *
                                   "Original error: $(e.msg)"))
        else
            rethrow()
        end
    end

    lhs_adj_Q_host = Array{ComplexF32}(undef, Nobs, max_iters_global, N_chains)
    copyto!(lhs_adj_Q_host, lhs_adj_Q_dev)

    Threads.@threads for iq_linear in 1:Nq
        iq = CartesianIndices(qpts.qs)[iq_linear]
        q = qpts.qs[iq]
        q_global = cryst.recipvecs * q

        corrbuf = zeros(ComplexF64, Ncorr)

        for ξ in 1:Nobs
            c = (iq_linear - 1) * Nobs + ξ
            chain_skip[c] && continue
            c_norm_per_chain[c] == 0 && continue

            n = n_iters_done_per_chain[c]

            αs_chain = view(αs_matrix, 1:n, c)
            βs_chain = n > 1 ? view(βs_matrix, 1:n-1, c) : Float64[]

            tridiag = SymTridiagonal(collect(αs_chain), collect(βs_chain))

            (; values, vectors) = try
                eigen(tridiag)
            catch
                eigen(Hermitian(collect(tridiag)))
            end

            c_norm = c_norm_per_chain[c]

            lhs_adj_Q_chain = ComplexF64.(view(lhs_adj_Q_host, :, 1:n, c))

            for (iω, ω) in enumerate(energies)
                f(x) = kernel(x, ω) * Sunny.thermal_prefactor(x; kT)

                corr_ξ = c_norm * lhs_adj_Q_chain * vectors * Diagonal(f.(values)) * (vectors'[:, 1])

                corrbuf .= 0
                for (i, (μ, ν)) in enumerate(measure.corr_pairs)
                    ξ == ν && (corrbuf[i] += (1/2) *     (corr_ξ[μ] / Ncells))
                    ξ == μ && (corrbuf[i] += (1/2) * conj(corr_ξ[ν] / Ncells))
                end

                data[iω, iq] += measure.combiner(q_global, corrbuf)
            end
        end
    end

    return Sunny.Intensities(cryst, qpts, collect(energies), data)
end
