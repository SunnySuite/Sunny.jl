using LinearAlgebra: BlasFloat, checksquare

using CUDA: unsafe_free!, with_workspaces, findfirst
using CUDA.CUBLAS: unsafe_strided_batch, handle, cublasZtrsmBatched_64, cublasCtrsmBatched_64
using CUDA.CUSOLVER: dense_handle, CuSolverParameters, cusolverDnZpotrfBatched, cusolverDnCpotrfBatched, cusolverDnXsyevBatched_bufferSize, cusolverDnXsyevBatched
## (TR) triangular triangular matrix solution batched
for (fname, elty) in ((:cublasZtrsmBatched_64, :ComplexF64),
                      (:cublasCtrsmBatched_64, :ComplexF32))
    @eval begin
        function trsm_batched!(side::Char,
                               uplo::Char,
                               transa::Char,
                               diag::Char,
                               alpha,
                               n,
                               lda,
                               batch_size,
                               A::CuArray{CuPtr{$elty}, 1},
                               B::CuArray{CuPtr{$elty}, 1})
            $fname(handle(), side, uplo, transa, diag, n, n, alpha, A, lda, B, lda, batch_size)
        end
    end
end

for (fname, elty) in ((:cusolverDnCpotrfBatched, :ComplexF32),
                      (:cusolverDnZpotrfBatched, :ComplexF64))
    @eval begin
        function potrfBatched!(dh, uplo::Char, n, lda, batch_size, A::CuArray{CuPtr{$elty}, 1})
            # Run the solver
            $fname(dh, uplo, n, A, lda, dh.info, batch_size)

            @assert isnothing(findfirst(!iszero, dh.info))
            return A
        end
    end
end

# XsyevBatched
function XsyevBatched!(dh, jobz::Char, uplo::Char, n, lda, batch_size, A::StridedCuArray{T, 3}) where {T <: BlasFloat}
    minimum_version = v"11.7.1"
    CUSOLVER.version() < minimum_version && throw(ErrorException("This operation requires cuSOLVER
        $(minimum_version) or later. Current cuSOLVER version: $(CUSOLVER.version())."))
    R = real(T)
    W = CuMatrix{R}(undef, n, batch_size)
    params = CuSolverParameters()

    function bufferSize()
        out_cpu = Ref{Csize_t}(0)
        out_gpu = Ref{Csize_t}(0)
        cusolverDnXsyevBatched_bufferSize(
            dh, params, jobz, uplo, n,
            T, A, lda, R, W, T, out_gpu, out_cpu, batch_size
        )
        return out_gpu[], out_cpu[]
    end

    with_workspaces(dh.workspace_gpu, dh.workspace_cpu, bufferSize()...) do buffer_gpu, buffer_cpu
        cusolverDnXsyevBatched(
            dh, params, jobz, uplo, n, T, A,
            lda, R, W, T, buffer_gpu, sizeof(buffer_gpu),
            buffer_cpu, sizeof(buffer_cpu), dh.info, batch_size
        )
    end

    @assert isnothing(findfirst(!iszero, dh.info))

    return W
end

function hegvd_batched!(H_d, I_d)
    n = checksquare(H_d)
    batch_size = size(H_d, 3)
    @assert size(H_d) == size(I_d)

    lda = max(1,stride(H_d, 2))
    @assert max(1,stride(I_d, 2)) == lda

    H_dp = unsafe_strided_batch(H_d)
    I_dp = unsafe_strided_batch(I_d)

    dh = dense_handle()
    resize!(dh.info, batch_size)

    potrfBatched!(dh, 'L', n, lda, batch_size, H_dp)
    trsm_batched!('R', 'L', 'C', 'N', ComplexF64(1.), n, lda, batch_size, H_dp, I_dp)
    trsm_batched!('L', 'L', 'N', 'N', ComplexF64(1.), n, lda, batch_size, H_dp, I_dp)
    evalues_d = XsyevBatched!(dh, 'V', 'L', n, lda, batch_size, I_d)
    trsm_batched!('L', 'L', 'C', 'N', ComplexF64(1.), n, lda, batch_size, H_dp, I_dp)
    unsafe_free!(H_dp)
    unsafe_free!(I_dp)
    return evalues_d
end
