"""
    SpinWaveTheoryLarge(sys::System; measure, tol=1e-2, regularization=1e-8, maxiters=1000)

A variant of [`SpinWaveTheory`](@ref) optimized for calculations on large
magnetic cells. This method avoids direct matrix diagonalization, and thereby
reduces computational cost from cubic- to linear-scaling in the system size.
Large supercells can be useful for modeling, e.g., quenched chemical disorder or
magnetic states with long-wavelength modulations.

The [`intensities`](@ref) are estimated independently per ``𝐪`` using the
Lanczos method [1, 2]. Lanczos iterations continue until a heuristic error
estimate is below `tol` relative to the max intensity at ``𝐪``. The number of
iterations scales like ``(Δϵ/w) \\log(1/\\mathrm{tol})``, where ``Δϵ`` is the
spectral bandwidth of excitations and ``w ≥ η`` is the required energy
resolution. Small ``w`` may arise from nearly degenerate soft modes; here,
consider using the intrinsic broadening ``η`` as a regulator. Exceeding
`maxiters` iterations will throw an error. Pass `verbose=true` to `intensities`
to print a progress bar and the iteration statistics.

!!! tip "Consider `SampledCorrelations` when calculating powder averages"

    Spin wave theory performs an independent calculation for each ``𝐪`` point of
    interest. Consequently, it can be very slow to sample a 3D volume of
    ``𝐪``-space, e.g., as required for a [`powder_average`](@ref). Here a
    compelling alternative may be [`SampledCorrelations`](@ref). It uses real-time
    spin dynamics to calculate structure factor data over the entire 3D grid of
    commensurate ``𝐪``-vectors in one shot. This may provide a considerable speedup
    at the cost of limited ``𝐪``-space resolution and statistical sampling error.

Lanczos is near optimal at a fixed number of iterations [2]. A related scheme is
the Kernel Polynomial Method [3], which Sunny had previously implemented as
`SpinWaveTheoryKPM`.

## References

1. [T. Chen, _The Lanczos algorithm for matrix functions: a handbook for
   scientists_ (2024) [arXiv:2410.11090]](https://arxiv.org/abs/2410.11090).
2. [N. Amsel, T. Chen, A. Greenbaum, C. Musco, C. Musco, _Near-Optimal
   Approximation of Matrix Functions by the Lanczos Method_ (2023)
   [arXiv:2303.03358]](https://arxiv.org/abs/2303.03358).
3. [H. Lane et al., _Kernel Polynomial Method for Linear Spin Wave Theory_
   (2023) [arXiv:2312.08349]](https://arxiv.org/abs/2312.08349).
"""
struct SpinWaveTheoryLarge
    swt :: SpinWaveTheory
    H :: SparseDynamicalMatrix
    tol :: Float64
    maxiters :: Int

    function SpinWaveTheoryLarge(sys::System; measure::Union{Nothing, MeasureSpec},
                                 regularization=1e-8, tol=1e-2, maxiters=1000)
        0 < tol < 1 || error("Require 0 < tol < 1")
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

function intensities(swt_large::SpinWaveTheoryLarge, qpts; energies, η=0.0, kernel=nothing, kT=0.0, threaded=false, verbose=false)
    qpts = convert(AbstractQPoints, qpts)
    data = zeros(eltype(swt_large.swt.measure), length(energies), size(qpts.qs)...)
    return intensities!(data, swt_large, qpts; energies, η, kernel, kT, threaded, verbose)
end

function intensities!(data, swt_large::SpinWaveTheoryLarge, qpts; energies, η=0.0, kernel=nothing, kT=0.0, threaded=false, verbose=false)
    qpts = convert(AbstractQPoints, qpts)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    η >= 0 || error("Require η ≥ 0")
    if iszero(η)
        # The limit η → 0⁺, approximated by η = √ε × kernel width, for machine
        # precision ε. A Lorentzian of half width η biases intensities by about
        # η/width relative to the peak, here √ε, below any attainable `tol`. The
        # floor is needed because, at kT > 0, the pole of 1 + n_B at ω = 0 is
        # evaluated at z = iη. There, χ(z) = (G + β₀² x₁)/z loses about ε |G|/η
        # to cancellation, i.e., √ε |G|/width, which is negligible.
        isnothing(kernel) && error("Specify a `kernel`, a broadening `η`, or both.")
        # A measure that is not real-valued needs cross correlations between
        # observables, which lose accuracy as z → 0, from cancellation in the
        # probe sums of `LanczosResolvent`
        if !(eltype(swt_large.swt.measure) <: Real) && kT > 0
            error("Specify a broadening η > 0 for a measure that is not real-valued at nonzero kT")
        end
        η = sqrt(eps()) * minimum(ω -> kernel_fwhm(kernel, ω), energies)
    end
    intensities_lanczos!(data, swt_large, qpts; energies, η, kernel, kT, threaded, verbose)
    return Intensities(orig_crystal(swt_large.swt.sys), qpts, collect(Float64, energies), data)
end

# At each wavevector, Lanczos builds the tridiagonal matrix T that represents ĨH
# in the H inner product, starting from v = Ĩu, for a combination u of the
# observable vectors. The susceptibility of u at complex energy z follows from a
# partial-fraction identity,
#
#     χ(z) = Σₖ |wₖ|² / (λₖ (z - λₖ)) = (G + β₀² [(z - T)⁻¹]₁₁) / z,
#
# where the sum runs over eigenpairs of ĨH with weights wₖ, β₀² = v† H v, and
# the commutator G = Σₖ |wₖ|²/λₖ = v† Ĩ v is known exactly. Lanczos
# approximates only the resolvent, which is analytic within a distance Im z of
# the real axis. Pairs of exact zero modes at ±λ, e.g., of a finite magnetic
# cell with continuous symmetry, cancel inside χ'' as λ → 0.
#
# A run from v_b also gives the cross terms χ_ab(z) = (G_ab + β₀ dₐᵀ (z - T)⁻¹ e₁)
# / z, with the overlaps dₐⱼ = vₐ† H q_j of probes vₐ, and G_ab = vₐ† Ĩ v_b. A
# general combiner needs all correlations C, from one run per observable, with
# all observables as probes. A real-valued `combiner` is real-linear, so that it
# is tr(M C) for a Hermitian matrix M, found by probing. With the eigenpairs M =
# Σₖ μₖ eₖ eₖ†, the intensity is Σₖ μₖ eₖ† C eₖ, which needs only one run per
# nonzero eigenvalue, from u = Σₐ eₖₐ uₐ, and no probes, e.g., 2 runs for
# `ssf_perp`.
#
# 1. T. Chen, The Lanczos algorithm for matrix functions: a handbook for
#    scientists (2024), Sec. 5.
function intensities_lanczos!(data, swt_large, qpts; energies, η, kernel, kT, threaded, verbose)
    (; swt, H, tol, maxiters) = swt_large
    (; sys, measure) = swt
    cryst = orig_crystal(sys)

    @assert eltype(data) == eltype(measure)
    fill!(data, zero(eltype(data)))

    Ncells = nsites(sys) / natoms(cryst)
    n = size(H, 1)
    Nobs = num_observables(measure)
    rb = RetardedBroadening(energies; η, kernel, kT)
    resolution = max(rb.width, 2η)
    # Lanczos steps per wavevector, summed over runs, for the `verbose` report
    iters = zeros(Int, size(qpts.qs))

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
        mulS!(y, x) = mul_sparse!(y, H, x, Hq)

        V = mul_Ĩ!(similar(u), u)

        # The starting vectors Ĩu of the runs, and whether they are also probes
        if eltype(measure) <: Real
            (μs, es) = eigen(combiner_matrix(measure, q_global))
            ks = findall(μ -> abs(μ) > 1e-12 * maximum(abs, μs), μs)
            μs = μs[ks]
            # Σₐ eₖₐ vₐ by broadcast, which avoids BLAS
            starts = [sum(es[a, k] .* view(V, :, a) for a in 1:Nobs) for k in ks]
            cross = false
        else
            starts = collect(eachcol(V))
            cross = true
        end

        (As, m) = try
            lanczos_correlations(mulS!, starts, rb; cross, resolution, tol, maxiters)
        catch err
            err isa InstabilityError || rethrow()
            throw(InstabilityError("Not an energy-minimum; wavevector q = $(vec3_to_string(q)) unstable."))
        end
        if isnothing(As)
            error("Exceeded maxiters = $maxiters at wavevector q = $(vec3_to_string(q)); \
                   consider increasing η or maxiters")
        end
        iters[iq] = m

        # The structure factor -(χ - χ†)/2πi per chemical cell. For a real
        # combiner, from the start of each run, and otherwise from the probes,
        # where run b gives χ[:, b].
        for iω in axes(data, 1)
            data[iω, iq] = if eltype(measure) <: Real
                sum(-μs[k] * imag(As[k][iω, 1]) / π for k in eachindex(μs)) / Ncells
            else
                χ = [As[b][iω, 1+a] for a in 1:Nobs, b in 1:Nobs]
                C = -(χ - χ') / (2π*im*Ncells)
                map!(((a, b),) -> C[a, b], corr, measure.corr_pairs)
                measure.combiner(q_global, corr)
            end
        end
    end

    verbose && println("SpinWaveTheoryLarge intensities\n  dimension       $n, $Nobs observables")
    desc = verbose ? "  wavevectors     " : nothing
    elapsed = @elapsed foreach_chunked(calc_iq!, newbuf, CartesianIndices(qpts.qs); threaded, desc)

    if verbose
        nthreads = threaded ? min(Threads.nthreads(), length(qpts.qs)) : 1
        per_q = number_to_simple_string(1000 * elapsed / length(qpts.qs); digits=2)
        println("  elapsed         $(number_to_simple_string(elapsed; digits=2)) s on $nthreads \
                 thread$(nthreads == 1 ? "" : "s"), $per_q ms per 𝐪")
        if !isempty(iters)
            mean = number_to_simple_string(Statistics.mean(iters); digits=3)
            median = number_to_simple_string(Statistics.median(iters); digits=3)
            println("  Lanczos steps   mean $mean, median $median, max $(maximum(iters))")
        end
    end
    return data
end

# The broadened χ(z), as Σⱼ Wᵢⱼ χ(zⱼ) at the energies of `rb`, from Lanczos runs
# for each of the starting vectors `vs`. For run k, column 1 of As[k] is χ of
# the start. With `cross`, the starts are also the probes of every run, and
# column 1 + a is the cross term with start a. Starts that are negligible
# relative to the largest, e.g., that vanish at this wavevector, give zero,
# rather than iterating on roundoff. Returns also the number of Lanczos steps
# taken, summed over runs, or `nothing` for As if not converged within
# `maxiters` steps.
#
# Iteration stops when both:
# 1. Every run has enough rows to resolve the energy scale `resolution` = w, the
#    larger of the kernel width and the line width 2η, unless its Krylov space
#    is exhausted. Like the nodes of Gauss quadrature, the m Ritz values of T
#    are spaced by about πΛ/m near the center of a spectrum in [-Λ, Λ], i.e.,
#    near ω = 0, where soft modes are. Before m ≥ πΛ/w, the rows so far say
#    little about modes a distance w apart, and the error estimate, which
#    assumes that the further rows resemble those so far, cannot be trusted.
#    Here Λ is estimated by 2 max β_j, which bounds the Ritz values when α_j =
#    0, by the Gershgorin theorem.
# 2. Every error estimate of `LanczosResolvent`, for the start and the probes,
#    broadened by the magnitudes of the same weights, is below `tol/2` relative
#    to the largest intensity -Im χ / π of the starts, a factor 2 margin for the
#    estimates.
#
# The runs advance together, because this largest intensity is itself an
# estimate, which can fall by orders of magnitude while some run is unconverged,
# e.g., at a soft mode as η → 0.
#
# Checks are at even numbers of rows of T. For a real combination of Hermitian
# observables, the weights at ±λ are equal, so that α_j = 0, and the spectrum of
# T is symmetric under θ → -θ. An odd number of rows then forces a spurious Ritz
# value at zero. The steps between checks grow by 25%, so that the overhead of
# checks is a fixed fraction.
function lanczos_correlations(mulS!, vs, rb; cross, resolution, tol, maxiters)
    probes = cross ? reduce(hcat, vs) : zeros(ComplexF64, length(vs[1]), 0)
    As = [zeros(ComplexF64, size(rb.W, 1), 1 + size(probes, 2)) for _ in vs]
    norms = [sqrt(real(dot_noblas(v, v))) for v in vs]
    ks = findall(>(sqrt(eps()) * maximum(norms)), norms)
    rs = [LanczosResolvent(mulS!, vs[k], probes, rb.zs) for k in ks]
    isempty(rs) && return (As, 0)
    # The largest H-norm of a probe, for the bound on cross terms
    pnorm = cross ? maximum(r -> r.β₀, rs) : 0.0
    errs = zeros(length(rs))
    m = 10
    while true
        for (k, r) in enumerate(rs)
            extend!(r, m)
            errs[k] = broaden!(As[ks[k]], r, pnorm, rb)
        end
        scale = maximum(A -> maximum(x -> abs(imag(x)), view(A, :, 1)), As) / π
        steps = sum(r -> r.nrows, rs)
        resolved = all(r -> exhausted(r) || r.nrows >= 2π * r.βmax / resolution, rs)
        resolved && 2 * maximum(errs) <= tol * scale && return (As, steps)
        m >= maxiters && return (nothing, steps)
        m += max(2, m ÷ 4)
        m += isodd(m)
    end
end

# Writes into the columns of A the broadened χ = (G₁ + β₀² x₁) / z of the
# start, and χₐ = (G₁₊ₐ + β₀ sₐ) / z of each probe a. Returns the largest
# broadened error estimate, which for probes of H-norm up to `pnorm` includes
# the bound on cross terms.
function broaden!(A, r::LanczosResolvent, pnorm, rb::RetardedBroadening)
    (; zs, W, ranges) = rb
    (; β₀, G, x₁, s) = r
    fill!(A, 0)
    err = zeros(size(A, 1))
    @inbounds for j in eachindex(zs)
        (δx, δe) = error_estimates(r, j)
        δχ = max(β₀^2 * δx, β₀ * pnorm * δe) / abs(zs[j])
        χ = (G[1] + β₀^2 * x₁[j]) / zs[j]
        for i in ranges[j]
            A[i, 1] += W[i, j] * χ
            err[i] += abs(W[i, j]) * δχ
        end
        for a in axes(s, 1)
            χ = (G[1+a] + β₀ * s[a, j]) / zs[j]
            for i in ranges[j]
                A[i, 1+a] += W[i, j] * χ
            end
        end
    end
    return maximum(err; init=0.0) / π
end
