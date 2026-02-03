diffnorm2(a::Number, b::Number) = abs2(a - b)
function diffnorm2(a, b, i, L, tol)
    buf = CuDynamicSharedArray(Float64, (blockDim().x, blockDim().y))
    bufq = view(buf, :, threadIdx().y)

    acc = 0.0
    for j in 1:size(a,1)
        acc += diffnorm2(a[i, j], b[i ,j])
    end
    for j in 1:size(a,1)
        acc += diffnorm2(a[i+L,j], b[i+L,j])
    end
    bufq[threadIdx().x] = acc
    CUDA.sync_threads()

    tid = threadIdx().x
    # Perform reduction within the block
    stride = blockDim().x ÷ 2
    while stride >= 1
        if tid <= stride
            bufq[tid] += bufq[tid + stride]
        end
        CUDA.sync_threads() # Synchronize threads at each reduction step
        stride ÷= 2
    end

    if tid == 1
        @assert bufq[1] < tol
    end
    return
end

function matrix_cleanup(H, swt, L)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if i > L
        return
    end

    iq = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    if iq > size(H, 3)
        return
    end

    Hq = view(H,:,:,iq)

    @inbounds begin
        # H must be hermitian up to round-off errors
        diffnorm2(Hq, Hq', i, L, 1e-12)

        # Make H exactly hermitian
        for j in 1:i-1
            Hq[i, j] = val = 0.5 * (Hq[i, j] + adjoint(Hq[j, i]))
            Hq[j, i] = adjoint(val)
        end
        Hq[i,i] = real(Hq[i,i])

        ip = 2L - i + 1
        for jp in 1:ip -i
            Hq[ip, jp] = val = 0.5 * (Hq[ip, jp] + adjoint(Hq[jp, ip]))
            Hq[jp, ip] = adjoint(val)
        end
        Hq[ip,ip] = real(Hq[ip,ip])

        # Add small constant shift for positive-definiteness
        Hq[i,i] += swt.regularization
        Hq[i+L,i+L] += swt.regularization
    end
    return
end

function get_shmem_matrix_cleanup(threads)
    if length(threads) == 2
        return threads[1] * threads[2] * sizeof(Float64)
    else
        return threads * sizeof(Float64)
    end
end
