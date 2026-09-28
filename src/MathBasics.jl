# Frequently used static types
const Vec3 = SVector{3, Float64}
const Vec5 = SVector{5, Float64}
const Mat3 = SMatrix{3, 3, Float64, 9}
const Mat5 = SMatrix{5, 5, Float64, 25}
const CMat3 = SMatrix{3, 3, ComplexF64, 9}
const CVec{N} = SVector{N, ComplexF64}
const HermitianC64 = Hermitian{ComplexF64, Matrix{ComplexF64}}

# Convenience for Kronecker-δ syntax. Note that boolean (false, true) result
# acts as multiplicative (0, 1).
@inline δ(x, y) = (x==y)

# Square matrix times SVector without allocating. Unrolled via @generated so
# that the compiler has freedom to rearrange order of summation. Beats
# SMatrix{N, N}(H) * v in both run-time and compile-time benchmarks.
@generated function mul_svec(H::AbstractMatrix, v::SVector{N}) where N
    rows = [:( +($([:(H[$i,$j] * v[$j]) for j in 1:N]...)) ) for i in 1:N]
    return :(SVector{N}($(rows...)))
end

# Calculates norm(a)^2 without allocating
norm2(a::Number) = abs2(a)
function norm2(a)
    acc = 0.0
    for i in eachindex(a)
        acc += norm2(a[i])
    end
    return acc
end

# Calculates norm(a - b)^2 without allocating
diffnorm2(a::Number, b::Number) = abs2(a - b)
function diffnorm2(a, b)
    @assert size(a) == size(b) "Non-matching dimensions"
    acc = 0.0
    for i in eachindex(a)
        acc += diffnorm2(a[i], b[i])
    end
    return acc
end

# Calculates norm(a - b, Inf) without allocating
function maxdiff(a, b)
    @assert size(a) == size(b) "Non-matching dimensions"
    ret = 0.0
    for i in eachindex(a)
        ret = max(abs(a[i] - b[i]), ret)
    end
    return ret
end

function is_integer(x; tol)
    return abs(x - round(x)) < tol
end

function all_integer(xs; tol)
    return all(is_integer(x; tol) for x in xs)
end

# Periodic variant of Base.isapprox. When comparing lattice quantities like
# positions or bonds, prefer is_periodic_copy because it works element-wise.
function isapprox_mod1(x::AbstractArray, y::AbstractArray; opts...)
    @assert size(x) == size(y) "Non-matching dimensions"
    Δ = @. mod(x - y + 0.5, 1) - 0.5
    return isapprox(Δ, zero(Δ); opts...)
end

# Project `v` onto space perpendicular to `n`
@inline proj(v, n) = v - n * ((n' * v) / norm2(n))

# Avoid linter false positives per
# https://github.com/julia-vscode/julia-vscode/issues/1497
kron(a...) = Base.kron(a...)

function tracelesspart(A)
    @assert allequal(size(A))
    return A - tr(A) * I / size(A,1)
end

# https://github.com/JuliaLang/julia/issues/44996
function findfirstval(f, a)
    i = findfirst(f, a)
    return isnothing(i) ? nothing : a[i]
end

# Returns the QL decomposition `Q, L = ql(A)` satisfying `Q * L ≈ A` with Q
# orthogonal and L lower-triangular. 
#
# Let (Q, R) be the usual QR decomposition of A. Let F be the matrix with ones
# on the antidiagonal. Then AF is the matrix A with columns reversed and FRF is
# the matrix R with all elements reversed. With this notation, the return value
# (; Q=QF, L=FRF) gives the desired QL decomposition of A.
function ql(A)
    AF = reverse!(Matrix(A); dims=2)
    (; Q, R) = qr!(AF)
    QF = reverse!(Matrix(Q); dims=2)
    FRF = reverse!(R)
    return (; Q=QF, L=FRF)
end

# Perform the SVD decomposition A = U Σ Vᵀ and return an iterator over the
# principal triples (σₖ, U[:,k], V[:,k]), dropping components where σₖ < atol.
function svd_iterator(A; atol)
    F = svd(A)
    return ((σ, u, v) for (σ, u, v) in zip(F.S, eachcol(F.U), eachcol(F.V)) if σ > atol)
end

# flatten_to_vec([1, ([2, 3], 4, [5, 6])]) == [1, 2, 3, 4, 5, 6]
flatten_to_vec(x::Number) = [x]
flatten_to_vec(x::AbstractArray{<: Number}) = vec(x)
flatten_to_vec(xs) = reduce(vcat, (flatten_to_vec(x) for x in xs))

same_shape(x::Number, y::Number) = true
same_shape(x::AbstractArray{<:Number}, y::AbstractArray{<:Number}) = size(x) == size(y)
function same_shape(xs, ys)
    size(xs) == size(ys) || return false
    all(same_shape(x, y) for (x, y) in zip(xs, ys))
end


# Rescale v such that sum(v) = 1
fractionalize(v) = iszero(v) ? one.(v) / length(v) : v ./ sum(v)

# Student's t-distribution, but normalized to 1 at x=0. Converges to exp(-x²/2)
# when ν → ∞.
function studentt_kernel(x::Real, ν::Real)
    ν > 0 || error("ν must be positive")
    if isinf(ν)
        return exp(-x^2/2)
    else
        return exp(-((ν+1)/2) * log1p(x^2/ν))
    end
end

"""
    softplus(x; β=1) = log(1 + exp(β x)) / β

Smooth approximation to `max(x, 0)`, exact in the limit `β = Inf`.
"""
function softplus(x; β=1)
    β > 0 || error("β must be positive")
    t = β*x
    if t > 40
        return x                 # log(1+exp(t)) ~ t
    elseif t < -40
        return exp(t) / β        # log(1+exp(t)) ~ exp(t)
    else
        return log1p(exp(t)) / β
    end
end

"""
    softcap(x, cap; β=1) = cap - softplus(cap - x; β)

Smooth approximation to `min(x, cap)`, exact in the limit `β = Inf`.
"""
softcap(x, cap; β=1) = cap - softplus(cap - x; β)

# Recommended alternative to OpenBLAS, or `nothing` if none for this platform
fast_blas_backend() = Sys.isapple() ? :AppleAccelerate : Sys.ARCH === :x86_64 ? :MKL : nothing

"""
    load_fast_blas()

Loads a BLAS backend for the whole Julia process: either AppleAccelerate or MKL,
depending on the platform. This significantly accelerates spin-wave
[`intensities`](@ref) calculations with the `threaded=true` option.

Julia's default backend is OpenBLAS. It performs well serially, but suffers from
global lock contention when parallelizing over many small matrix calculations.
"""
function load_fast_blas()
    backend = fast_blas_backend()
    isnothing(backend) &&
        error("Cannot recommend an alternative to OpenBLAS for this platform: ", Sys.MACHINE)

    isnothing(Base.find_package(String(backend))) &&
        error("Backend $backend is recommended; install it in the Julia package manager.")

    Base.eval(Main, :(using $backend))
    using_openblas() && error("Loaded $backend, but OpenBLAS is still the BLAS backend.")

    println("Loaded $backend as the BLAS backend.")
    return nothing
end

# Whether OpenBLAS is the library servicing Julia's ILP64 BLAS calls
function using_openblas()
    lib = BLAS.lbt_find_backing_library("zgemm_", :ilp64)
    return !isnothing(lib) && occursin("openblas", lowercase(basename(lib.libname)))
end

# Calls `f(buf, i)` for each `i` in `indices`, where `buf = newbuf()` is a
# buffer owned by the calling task. If `threaded`, spawns one task per thread.
# These dynamically claim blocks of indices from a shared counter, which
# balances the load when indices vary in cost or cores vary in speed (e.g.,
# performance vs. efficiency cores). About 16 blocks per task keeps the load
# balanced, while amortizing the cost of claiming a block. This scheme is
# equivalent to OhMyThreads.jl's `GreedyScheduler(; chunking=true)` with buffers
# held in a `TaskLocalValue`, and performed similarly on benchmarks.
# Set `warn_blas` if `f` calls BLAS, to flag the poor thread scaling of OpenBLAS.
function foreach_chunked(f, newbuf, indices; threaded, warn_blas=false)
    if threaded && Threads.nthreads() == 1
        @warn "Option `threaded=true` has no effect here. Restart with `julia --threads=auto`." maxlog=1
    end
    if warn_blas && threaded && using_openblas()
        backend = fast_blas_backend()
        if isnothing(backend)
            @warn "OpenBLAS scales poorly with `threaded=true` (but cannot recommend alternative for $(Sys.MACHINE))" maxlog=1
        else
            @warn "OpenBLAS scales poorly with `threaded=true` (consider loading $backend)" maxlog=1
        end
    end
    if threaded
        n = length(indices)
        ntasks = min(n, Threads.nthreads())
        blocksize = max(1, n ÷ 16ntasks)
        next = Threads.Atomic{Int}(1)
        try
            @sync for _ in 1:ntasks
                Threads.@spawn let buf = newbuf()
                    while (start = Threads.atomic_add!(next, blocksize)) <= n
                        for j in start:min(start + blocksize - 1, n)
                            f(buf, indices[j])
                        end
                    end
                end
            end
        catch err
            # Unwrap task failure to preserve its type, e.g., `InstabilityError`
            throw(err isa CompositeException ? first(err).task.exception : err)
        end
    else
        buf = newbuf()
        foreach(i -> f(buf, i), indices)
    end
end
