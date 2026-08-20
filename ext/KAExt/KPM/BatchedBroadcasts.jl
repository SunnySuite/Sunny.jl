@kernel inbounds=true function batched_w_update_kernel_ka!(w, α, vp, β, v)
    tid = @index(Global, Linear)
    N_chains = size(w, 1)
    CT = @uniform eltype(w)

    c = ((tid - 1) % N_chains) + 1
    k = ((tid - 1) ÷ N_chains) + 1
    w[c, k] = CT(w[c, k] - real(α[c]) * vp[c, k] - β[c] * v[c, k])
end

function batched_w_update_ka!(
    w::AbstractMatrix{<:Complex},
    α::AbstractVector{<:Complex},
    vp::AbstractMatrix{<:Complex},
    β::AbstractVector,
    v::AbstractMatrix{<:Complex},
)
    total = length(w)
    backend = KernelAbstractions.get_backend(w)
    workgroup = min(total, 256)
    workgroup = workgroup < 1 ? 1 : workgroup
    kernel! = batched_w_update_kernel_ka!(backend, workgroup)
    kernel!(w, α, vp, β, v; ndrange=total)
    return w
end
