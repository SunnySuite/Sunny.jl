# Recommended alternative to OpenBLAS for threaded workloads, or `nothing` if none
# for this platform. Serial workloads may run faster with OpenBLAS.
threading_blas_backend() = Sys.isapple() ? :AppleAccelerate : Sys.ARCH === :x86_64 ? :MKL : nothing

"""
    load_blas_for_threading()

Loads a BLAS backend that scales better under multi-threaded workloads.
Specifically, this function replaces Julia's OpenBLAS default with either MKL or
AppleAccelerate, depending on the platform.

Changing the BLAS backend is especially recommended for [`SpinWaveTheory`](@ref)
calculations whenever `threaded=true` is enabled. Although OpenBLAS is fast for
serial workloads, it employs a global thread lock that can strongly bottleneck
parallelized operations on small matrices.
"""
function load_blas_for_threading()
    backend = threading_blas_backend()
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
# held in a `TaskLocalValue`, and performed similarly on benchmarks. Set
# `warn_blas` if `f` calls BLAS with small matrices, to flag the poor thread
# scaling of OpenBLAS. A progress bar labeled `desc` advances per block, unless
# `desc` is nothing.
function foreach_chunked(f, newbuf, indices; threaded, warn_blas=false, desc=nothing)
    if threaded && Threads.nthreads() == 1
        @warn "Option `threaded=true` has no effect here. Restart with `julia --threads=auto`." maxlog=1
    end
    if warn_blas && threaded && using_openblas()
        backend = threading_blas_backend()
        if isnothing(backend)
            @warn "OpenBLAS scales poorly with `threaded=true` (but cannot recommend alternative for $(Sys.MACHINE))" maxlog=1
        else
            @warn "OpenBLAS scales poorly with `threaded=true` (consider loading $backend)" maxlog=1
        end
    end

    n = length(indices)
    ntasks = threaded ? clamp(n, 1, Threads.nthreads()) : 1
    # A single task has no contention to amortize, and advances the progress
    # bar per index
    blocksize = threaded ? max(1, n ÷ 16ntasks) : 1
    next = Threads.Atomic{Int}(1)
    # To animate, the bar must rewrite the line with `\r`, which requires a TTY
    enabled = !isnothing(desc) && stdout isa Base.TTY
    meter = ProgressMeter.Progress(n; desc=something(desc, ""), output=stdout, enabled)
    function work()
        buf = newbuf()
        while (start = Threads.atomic_add!(next, blocksize)) <= n
            stop = min(start + blocksize - 1, n)
            for j in start:stop
                f(buf, indices[j])
            end
            ProgressMeter.next!(meter; step=stop-start+1)
        end
    end

    if threaded
        try
            @sync for _ in 1:ntasks
                Threads.@spawn work()
            end
        catch err
            # Unwrap task failure to preserve its type, e.g., `InstabilityError`
            throw(err isa CompositeException ? first(err).task.exception : err)
        end
    else
        work()
    end
end
