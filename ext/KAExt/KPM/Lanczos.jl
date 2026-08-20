using LinearAlgebra: SymTridiagonal, dot, mul!, eigmax, eigmin

function lanczos_device_ka(
    mulA!,
    mulS!,
    v::AbstractVector{ComplexF64},
    vp::AbstractVector{ComplexF64},
    Sv::AbstractVector{ComplexF64},
    Svp::AbstractVector{ComplexF64},
    w::AbstractVector{ComplexF64},
    Sw::AbstractVector{ComplexF64},
    lhs::AbstractMatrix{ComplexF64},
    corr_dev::AbstractVector{ComplexF64},
    corr_host::Vector{ComplexF64},
    lhs_adj_Q_buf::Matrix{ComplexF64};
    min_iters::Int,
    max_iters::Int,
    resolution::Float64 = Inf,
    verbose::Bool = false,
)::Tuple{SymTridiagonal{Float64, Vector{Float64}}, Matrix{ComplexF64}}

    Ncorr = size(lhs, 2)
    @assert size(lhs_adj_Q_buf, 1) == Ncorr "lhs_adj_Q_buf row count must match size(lhs, 2)"
    @assert size(lhs_adj_Q_buf, 2) >= max_iters "lhs_adj_Q_buf column count too small"
    @assert length(corr_dev) == Ncorr "corr_dev length must match size(lhs, 2)"
    @assert length(corr_host) == Ncorr "corr_host length must match size(lhs, 2)"
    @assert length(v) == length(vp) == length(Sv) == length(Svp) == length(w) == length(Sw)
    @assert size(lhs, 1) == length(v)

    αs = Float64[]
    βs = Float64[]
    n_iters_done = 0

    mulS!(Sv, v)

    vSv = dot(v, Sv)

    hermitian_atol = max(length(v), 1) * 1e-10
    abs(imag(vSv)) <= hermitian_atol || error("S not Hermitian (imag(v†Sv) = $(imag(vSv)), |imag| > $hermitian_atol)")
    abs(real(vSv) - 1) <= hermitian_atol || error("Initial v not normalized (real(v†Sv) = $(real(vSv)))")

    norm_factor = sqrt(real(vSv))
    @. v  /= norm_factor
    @. Sv /= norm_factor

    mulA!(w, Sv)

    α = real(dot(w, Sv))

    @. w = w - α * v
    mulS!(Sw, w)

    push!(αs, α)
    n_iters_done += 1

    mul!(corr_dev, lhs', v)
    copyto!(corr_host, corr_dev)
    lhs_adj_Q_buf[:, n_iters_done] .= corr_host

    niters_eff = max_iters
    @inbounds for i in 1:length(v)-1
        if i == min_iters
            T_partial = SymTridiagonal(αs, βs)
            Δϵ = eigmax(T_partial) - eigmin(T_partial)
            niters_from_resolution = max(min_iters, fld(Δϵ, resolution))
            niters_from_resolution += mod(niters_from_resolution, 2)
            niters_eff = min(max_iters, Int(niters_from_resolution))
            if verbose
                println("Δϵ=$Δϵ, niters_eff=$niters_eff (capped by max_iters=$max_iters)")
            end
        end

        i >= niters_eff && break

        β² = real(dot(w, Sw))
        iszero(β²) && break
        β² < 0 && error("S is not a positive definite measure (β² = $β² < 0)")

        β = sqrt(β²)

        @. vp  = w  / β
        @. Svp = Sw / β

        mulA!(w, Svp)
        α = real(dot(w, Svp))
        @. w = w - α * vp - β * v
        mulS!(Sw, w)

        @. v = vp

        push!(αs, α)
        push!(βs, β)
        n_iters_done += 1

        mul!(corr_dev, lhs', v)
        copyto!(corr_host, corr_dev)
        lhs_adj_Q_buf[:, n_iters_done] .= corr_host
    end

    T = SymTridiagonal(αs, βs)

    lhs_adj_Q = lhs_adj_Q_buf[:, 1:n_iters_done]

    return T, lhs_adj_Q
end
