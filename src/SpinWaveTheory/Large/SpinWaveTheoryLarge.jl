"""
    SpinWaveTheoryLarge(sys::System; measure, tol=1e-2, regularization=1e-8, maxiters=nothing)

A variant of [`SpinWaveTheory`](@ref) that estimates [`intensities`](@ref) using
iterated matrix-vector products. By avoiding direct matrix diagonalization, this
method reduces computational cost from cubic to linear-scaling in the system
size ``N``. Large system sizes can arise, e.g., for models of quenched disorder
or models with nearly incommensurate ordering wavevectors.

At each wavevector, iteration continues until the intensities agree with their
previous estimate to within the relative tolerance `tol`. The number of
iterations scales like ``(Δϵ/w) \\log(1/\\mathrm{tol})``, where ``Δϵ``
is the spectral bandwidth of excitations and ``w`` is the energy resolution of
interest, and the cost of each scales like ``𝒪(N)``. An optional `maxiters`
caps the iterations at each wavevector. Long-range dipole-dipole interactions
are not yet supported.

!!! tip "Consider `SampledCorrelations` when calculating powder averages"

    Spin wave theory requires an independent calculation for each ``𝐪`` point of
    interest. Consequently, it can be very slow to sample a 3D volume of
    ``𝐪``-space, e.g., as required for a [`powder_average`](@ref). A compelling
    alternative may be [`SampledCorrelations`](@ref). It uses real-time spin
    dynamics to calculate structure factor data over the entire 3D grid of
    commensurate ``𝐪``-vectors in one shot. This may provide a considerable speedup
    at the cost of: limited ``𝐪``-space resolution and stochastic error due to
    statistical sampling.

The calculation uses the block Lanczos method [1] to approximate the resolvent
of the dynamical matrix, which is near-optimal at a fixed number of iterations
[2, 3]. It supersedes the Kernel Polynomial Method of Ref. [4], and the former
`SpinWaveTheoryKPM`.

## References

1. [G. H. Golub and R. Underwood, _The block Lanczos method for computing
   eigenvalues_, in Mathematical Software III, 361-377
   (1977)](https://doi.org/10.1016/B978-0-12-587260-7.50018-2).
2. [T. Chen, _The Lanczos algorithm for matrix functions: a handbook for
   scientists_ (2024) [arXiv:2410.11090]](https://arxiv.org/abs/2410.11090).
3. [N. Amsel, T. Chen, A. Greenbaum, C. Musco, C. Musco, _Near-Optimal
   Approximation of Matrix Functions by the Lanczos Method_ (2023)
   [arXiv:2303.03358]](https://arxiv.org/abs/2303.03358).
4. [H. Lane et al., _Kernel Polynomial Method for Linear Spin Wave Theory_
   (2023) [arXiv:2312.08349]](https://arxiv.org/abs/2312.08349).
"""
struct SpinWaveTheoryLarge
    swt :: SpinWaveTheory
    H :: SparseDynamicalMatrix
    tol :: Float64
    maxiters :: Int

    function SpinWaveTheoryLarge(sys::System; measure::Union{Nothing, MeasureSpec}, regularization=1e-8,
                                 tol=1e-2, maxiters=nothing)
        0 < tol < 1 || error("Require 0 < tol < 1")
        maxiters = @something maxiters typemax(Int)
        maxiters > 0 || error("Require maxiters > 0")
        swt = SpinWaveTheory(sys; measure, regularization)
        return new(swt, SparseDynamicalMatrix(swt), tol, maxiters)
    end
end

function Base.show(io::IO, ::MIME"text/plain", swt_large::SpinWaveTheoryLarge)
    (; swt) = swt_large
    printstyled(io, "SpinWaveTheoryLarge ", mode_to_str(swt.sys), "\n"; bold=true, color=:underline)
    println(io, "  ", nsites(swt.sys), " sites")
end

"""
    intensities(swt::SpinWaveTheoryLarge, qpts; energies, η, kernel=nothing, kT=0, threaded=false)

Like [`intensities`](@ref) for [`SpinWaveTheory`](@ref), with the intrinsic
broadening `η` required. The iterative method evaluates the susceptibility at
``ω + iη`` directly, and cannot resolve energies finer than `η`. The optional
`kernel` must be analytic, and its width may vary smoothly with energy.

The accuracy is set by `tol` in [`SpinWaveTheoryLarge`](@ref), and the cost by
the energy resolution, the larger of `η` and the kernel width, except near soft
modes, where it is `η`. Set `threaded=true` to parallelize over wavevectors. As
for `SpinWaveTheory`, an instability of the magnetic structure is an error, but
here only one that couples to the measured observables is detected.
"""
function intensities(swt_large::SpinWaveTheoryLarge, qpts; energies, η=nothing, kernel=nothing, kT=0.0, threaded=false)
    qpts = convert(AbstractQPoints, qpts)
    data = zeros(eltype(swt_large.swt.measure), length(energies), size(qpts.qs)...)
    return intensities!(data, swt_large, qpts; energies, η, kernel, kT, threaded)
end

function intensities!(data, swt_large::SpinWaveTheoryLarge, qpts; energies, η=nothing, kernel=nothing, kT=0.0, threaded=false)
    qpts = convert(AbstractQPoints, qpts)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    isnothing(η) && error("Keyword `η` is required, e.g., an intrinsic broadening small compared to the kernel width.")
    η > 0 || error("Require η > 0")
    intensities_lanczos!(data, swt_large, qpts; energies, η, kernel, kT, threaded)
    return Intensities(orig_crystal(swt_large.swt.sys), qpts, collect(Float64, energies), data)
end

# At each wavevector, block Lanczos over the observables builds the block
# tridiagonal matrix T that represents ĨH in the H inner product, starting from
# the block V = Ĩu of observables u. The susceptibility at complex energy z
# follows from a partial-fraction identity,
#
#     χ(z) = Σₖ wₖwₖ† / (λₖ (z - λₖ)) = (G + R₀† (z - T)⁻¹₁₁ R₀) / z,
#
# where the sum runs over eigenpairs of ĨH with weights wₖ, and the commutator G
# = Σₖ wₖwₖ†/λₖ = V† Ĩ V is known exactly. Lanczos approximates only the
# resolvent, which is analytic within a distance Im z of the real axis, and is
# evaluated by a block continued fraction. Pairs of exact zero modes at ±λ,
# e.g., of a finite magnetic cell with continuous symmetry, cancel inside χ'' as
# λ → 0.
#
# After m steps, Lanczos resolves the spectrum to about b log(1/tol) / 4m, where
# b is the bandwidth [1]. The resolvent is evaluated at Im z = η. Iteration
# starts at the steps that resolve the larger of η and half the kernel width.
# Every 50% of extra steps, the result after thermal and instrumental broadening
# is checked, until it changes by less than `tol` relative to its maximum.
# Lanczos converges at the coarser scale set by the kernel, except near a soft
# mode, where the thermal factor exposes the structure of χ'' at ω → 0, and the
# scale is η.
#
# 1. T. Chen, The Lanczos algorithm for matrix functions: a handbook for
#    scientists (2024), Sec. 5.
function intensities_lanczos!(data, swt_large, qpts; energies, η, kernel, kT, threaded)
    (; swt, H, tol, maxiters) = swt_large
    (; sys, measure) = swt
    cryst = orig_crystal(sys)

    @assert eltype(data) == eltype(measure)
    fill!(data, zero(eltype(data)))

    Ncells = nsites(sys) / natoms(cryst)
    n = size(H, 1)
    Nobs = num_observables(measure)
    rb = RetardedBroadening(energies, η, kernel, kT)
    unconverged = Threads.Atomic{Bool}(false)

    newbuf() = (
        u    = zeros(ComplexF64, n, Nobs),
        Hq   = similar(H.vals),
        corr = zeros(ComplexF64, length(measure.corr_pairs)),
    )

    function calc_iq!(buf, iq)
        (; u, Hq, corr) = buf
        q = qpts.qs[iq]
        q_global = cryst.recipvecs * q
        observable_vectors!(u, swt, q_global)
        values_at!(Hq, H, to_supercell_rlu(sys, q))

        # Starting block Ĩ u, as rows, for the Lanczos iteration in the measure
        # S = H. Observables that vanish at this wavevector are dropped.
        ξs = [ξ for ξ in 1:Nobs if !iszero(view(u, :, ξ))]
        isempty(ξs) && return
        V = zeros(ComplexF64, length(ξs), n)
        mul_Ĩ_rows!(V, transpose(view(u, :, ξs)))
        G = gram_rows(V, mul_Ĩ_rows!(similar(V), V))
        mulS!(W, X) = mul_rows!(W, H, X, Hq)

        C = try
            lanczos_correlations(mulS!, V, G, rb; η, tol, maxiters, unconverged)
        catch e
            qstr = vec3_to_string(q)
            e isa InstabilityError && throw(InstabilityError("Not an energy-minimum; wavevector q = $qstr unstable."))
            rethrow()
        end

        # Indices of the pair correlations among the observables kept
        pairs = map(measure.corr_pairs) do (μ, ν)
            (findfirst(==(μ), ξs), findfirst(==(ν), ξs))
        end
        for iω in eachindex(C)
            for (i, (a, b)) in enumerate(pairs)
                corr[i] = (isnothing(a) || isnothing(b)) ? 0.0 : C[iω][a, b] / Ncells
            end
            data[iω, iq] = measure.combiner(q_global, corr)
        end
    end

    foreach_chunked(calc_iq!, newbuf, CartesianIndices(qpts.qs); threaded)

    if unconverged[]
        @warn "Not converged within maxiters = $maxiters at some wavevectors. Increase `maxiters` or `η`."
    end

    return data
end

# The pair correlations C[iω][a, b] of the observables whose starting block is V,
# at each of the `energies`, by block Lanczos iteration to convergence
function lanczos_correlations(mulS!, V, G, rb; η, tol, maxiters, unconverged)
    lz = BlockLanczos(mulS!, V)

    # The bandwidth b, from a short run, sets the steps that resolve an energy
    # scale w, about b log(1/tol) / 4w
    extend!(lz, min(maxiters, 10))
    θ = eigvals(lanczos_matrix(lz))
    b = θ[end] - θ[1]
    estimate() = broadened_correlations(lz, G, rb)

    # Start optimistically, at the steps that resolve half the kernel width, or
    # half the steps that resolve η, and continue until converged. Beyond twice
    # the steps that resolve η, further steps change nothing.
    steps(w) = ceil(Int, b * log(1/tol) / 4w)
    m₀ = steps(max(rb.width/2, 2η))
    m_max = min(maxiters, 2 * steps(η))
    extend!(lz, clamp(m₀ - nsteps(lz), 0, m_max - nsteps(lz)))
    C = estimate()
    while !lz.exhausted && nsteps(lz) < m_max
        extend!(lz, min(max(2, nsteps(lz) ÷ 2), m_max - nsteps(lz)))
        C′ = estimate()
        converged = maximum(maximum(abs, x) for x in C′ - C) <= tol * maximum(maximum(abs, x) for x in C′)
        C = C′
        converged && return C
    end
    nsteps(lz) >= maxiters && !lz.exhausted && (unconverged[] = true)
    return C
end

# Correlations from the structure factor (1 + n_B(ω)) χ''(ω + iη), where χ'' =
# -(χ - χ†)/2πi, broadened on the precomputed `RetardedBroadening`
function broadened_correlations(lz::BlockLanczos, G, rb)
    # The block size is known only at runtime, so a function barrier specializes
    # the rest on the static block type
    p = size(G, 1)
    G = SMatrix{p, p, ComplexF64}(G)
    R₀ = SMatrix{size(lz.R₀)..., ComplexF64}(lz.R₀)
    sb = StaticBlocks(lz)
    return broaden_retarded(z -> susceptibility(sb, R₀, G, z), rb)
end

# The p×p matrix χ(z) = (G + R₀† X₁(z) R₀)/z
function susceptibility(sb::StaticBlocks{r}, R₀::SMatrix{r, p}, G::SMatrix{p, p}, z) where {r, p}
    X₁ = resolvent_leading_block(sb, z)
    return (G + R₀' * X₁ * R₀) / z
end
