const BATCHED_DOT_WG = 64

@kernel inbounds=true function batched_dot_chains_kernel_ka!(result, a, b, twoL)
    chain = @index(Group, Linear)
    lid   = @index(Local, Linear)
    groupsize = @uniform @groupsize()[1]

    lmem_real = @localmem Float64 (BATCHED_DOT_WG,)
    lmem_imag = @localmem Float64 (BATCHED_DOT_WG,)

    acc = ComplexF64(0)
    k = lid
    while k <= twoL
        acc += conj(ComplexF64(a[chain, k])) * ComplexF64(b[chain, k])
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
        result[chain] = ComplexF64(lmem_real[1] + lmem_real[2], lmem_imag[1] + lmem_imag[2])
    end
end

function batched_dot_chains_ka!(
    result::AbstractVector{ComplexF64},
    a::AbstractMatrix{<:Complex},
    b::AbstractMatrix{<:Complex},
)
    N_chains = size(a, 1)
    twoL = size(a, 2)
    @assert size(a) == size(b) "batched_dot_chains_ka! a and b must have the same shape"
    @assert length(result) == N_chains "batched_dot_chains_ka! result length must match N_chains"
    @assert eltype(result) == ComplexF64 "batched_dot_chains_ka! result must be ComplexF64"

    backend = KernelAbstractions.get_backend(a)
    kernel! = batched_dot_chains_kernel_ka!(backend, BATCHED_DOT_WG)
    kernel!(result, a, b, twoL; ndrange = N_chains * BATCHED_DOT_WG)
    return result
end
