using LinearAlgebra: SymTridiagonal, eigmax, eigmin

function lanczos_device_batched_ka(
    mulA!,
    mulS!,
    v_batch::AbstractMatrix{<:Complex},
    vp_batch::AbstractMatrix{<:Complex},
    Sv_batch::AbstractMatrix{<:Complex},
    Svp_batch::AbstractMatrix{<:Complex},
    w_batch::AbstractMatrix{<:Complex},
    Sw_batch::AbstractMatrix{<:Complex},
    lhs_dev::AbstractArray{<:Complex, 3},
    corr_dev::AbstractArray{ComplexF64, 3},
    corr_host::Array{ComplexF64, 3},
    α_chains_dev::AbstractVector{ComplexF64},
    β²_chains_dev::AbstractVector{ComplexF64},
    β_chains_dev::AbstractVector{<:Real},
    α_chains_host::Vector{ComplexF64},
    β²_chains_host::Vector{ComplexF64},
    lhs_adj_Q_buf::Array{ComplexF64, 3};
    max_iters::Int,
    min_iters::Int,
    resolution::Float64,
    verbose::Bool = false,
)::Tuple{Matrix{Float64}, Matrix{Float64}, Array{ComplexF64, 3}, Vector{Int}, Vector{Float64}}

    N_chains = size(v_batch, 1)
    twoL = size(v_batch, 2)
    Nobs = size(lhs_dev, 2)
    Nq   = size(lhs_dev, 3)

    @assert size(v_batch) == size(vp_batch) == size(Sv_batch) == size(Svp_batch) == size(w_batch) == size(Sw_batch) "Recurrence buffers must all be (N_chains, 2L)"
    @assert size(lhs_dev, 1) == twoL "lhs_dev first axis must equal 2L"
    @assert Nobs * Nq == N_chains "lhs_dev's Nobs*Nq must equal N_chains"
    @assert size(corr_dev) == (Nobs, Nobs, Nq) "corr_dev must be (Nobs, Nobs, Nq)"
    @assert size(corr_host) == (Nobs, Nobs, Nq) "corr_host must be (Nobs, Nobs, Nq)"
    @assert length(α_chains_dev) == N_chains
    @assert length(β²_chains_dev) == N_chains
    @assert length(β_chains_dev) == N_chains
    @assert length(α_chains_host) == N_chains
    @assert length(β²_chains_host) == N_chains
    @assert size(lhs_adj_Q_buf) == (Nobs, max_iters, N_chains) "lhs_adj_Q_buf must be (Nobs, max_iters, N_chains)"
    @assert max_iters >= 1 "max_iters must be >= 1"

    αs_matrix = Matrix{Float64}(undef, max_iters, N_chains)
    βs_matrix = Matrix{Float64}(undef, max_iters - 1, N_chains)
    n_iters_done_per_chain = zeros(Int, N_chains)
    chain_done = falses(N_chains)
    c_norm_per_chain = Vector{Float64}(undef, N_chains)
    niters_eff = fill(max_iters, N_chains)

    mulS!(Sv_batch, v_batch)

    batched_dot_chains_ka!(α_chains_dev, v_batch, Sv_batch)

    copyto!(α_chains_host, α_chains_dev)

    T_real = real(eltype(v_batch))
    hermitian_atol = max(twoL, 1) * max(1e-10, 100 * eps(T_real))
    near_zero_atol = max(twoL, 1) * max(1e-12, 100 * eps(T_real))
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

    batched_dot_chains_ka!(α_chains_dev, w_batch, Sv_batch)

    copyto!(α_chains_host, α_chains_dev)

    @inbounds for c in 1:N_chains
        αs_matrix[1, c] = chain_done[c] ? 0.0 : Float64(real(α_chains_host[c]))
    end

    w_batch .-= reshape(real.(α_chains_dev), N_chains, 1) .* v_batch

    mulS!(Sw_batch, w_batch)

    batched_projection_ka!(corr_dev, lhs_dev, v_batch)

    copyto!(corr_host, corr_dev)

    @inbounds for c in 1:N_chains
        iq = ((c - 1) ÷ Nobs) + 1
        ξ = ((c - 1) % Nobs) + 1
        for μ in 1:Nobs
            lhs_adj_Q_buf[μ, 1, c] = ComplexF64(corr_host[μ, ξ, iq])
        end
        n_iters_done_per_chain[c] = 1
    end

    @inbounds for i in 1:max_iters - 1
        batched_dot_chains_ka!(β²_chains_dev, w_batch, Sw_batch)

        copyto!(β²_chains_host, β²_chains_dev)

        for c in 1:N_chains
            if !chain_done[c]
                β²_real = real(β²_chains_host[c])
                if !(β²_real > 0.0)
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

        vp_batch  .= w_batch  ./ reshape(β_chains_dev, N_chains, 1)
        Svp_batch .= Sw_batch ./ reshape(β_chains_dev, N_chains, 1)

        mulA!(w_batch, Svp_batch)

        batched_dot_chains_ka!(α_chains_dev, w_batch, Svp_batch)

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

        copyto!(v_batch, vp_batch)

        batched_projection_ka!(corr_dev, lhs_dev, v_batch)

        copyto!(corr_host, corr_dev)

        @inbounds for c in 1:N_chains
            if !chain_done[c]
                iq = ((c - 1) ÷ Nobs) + 1
                ξ = ((c - 1) % Nobs) + 1
                for μ in 1:Nobs
                    lhs_adj_Q_buf[μ, i + 1, c] = ComplexF64(corr_host[μ, ξ, iq])
                end
            end
        end

        if verbose
            n_done = count(chain_done)
            println("lanczos_device_batched_ka iter $i: $n_done/$N_chains chains done")
        end
    end

    return αs_matrix, βs_matrix, lhs_adj_Q_buf, n_iters_done_per_chain, c_norm_per_chain
end
