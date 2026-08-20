@kernel inbounds=true function batched_projection_kernel_ka!(
    corr, lhs, v_batch, Nobs, Nq, twoL,
)
    tid = @index(Global, Linear)
    CT = @uniform eltype(corr)

    μ    = ((tid - 1) % Nobs) + 1
    rest =  (tid - 1) ÷ Nobs
    ξ    =  (rest % Nobs) + 1
    iq   =  (rest ÷ Nobs) + 1

    c = (iq - 1) * Nobs + ξ

    acc = CT(0)
    for k in 1:twoL
        acc += conj(CT(lhs[k, μ, iq])) * CT(v_batch[c, k])
    end
    corr[μ, ξ, iq] = acc
end

function batched_projection_ka!(
    corr::AbstractArray{<:Complex, 3},
    lhs::AbstractArray{<:Complex, 3},
    v_batch::AbstractMatrix{<:Complex},
)
    Nobs = size(lhs, 2)
    Nq   = size(lhs, 3)
    twoL = size(lhs, 1)
    total = Nobs * Nobs * Nq

    backend = KernelAbstractions.get_backend(corr)
    workgroup = min(total, 256)
    workgroup = workgroup < 1 ? 1 : workgroup
    kernel! = batched_projection_kernel_ka!(backend, workgroup)
    kernel!(corr, lhs, v_batch, Nobs, Nq, twoL; ndrange=total)
    return corr
end
