"""
    SpinWaveTheoryLarge(sys::System; measure, tol=1e-3, regularization=1e-8, maxiters=1000)

A variant of [`SpinWaveTheory`](@ref) optimized for calculations on large
magnetic cells. By avoiding direct matrix diagonalization, the computational
cost is reduced from cubic- to linear-scaling in the system size. Large
supercells can be useful for modeling quenched chemical disorder or magnetic
states with long-wavelength modulations.

At each ``𝐪``, the [`intensities`](@ref) are estimated using iterated
matrix-vector multiplications. Iterations continue until the estimated error is
below `tol` relative to the max intensity for the given `energies` grid. The
number of iterations scales like ``(Δϵ/w) \\log(1/\\mathrm{tol})``, where ``Δϵ``
is the spectral bandwidth of excitations and ``w`` is the required energy
resolution. Exceeding `maxiters` iterations will throw an error.

Pass `verbose=true` to `intensities` to print a progress bar and calculation
statistics.

The calculation is fastest for a `measure` that returns a real number, e.g.
[`ssf_perp`](@ref). A measure that returns a complex number or an array is
evaluated as one real-valued measure per real component of its value, at
proportionally higher cost.

!!! tip "Numerical challenges with near-gapless spectra"

    If the spectrum is near-gapless, and if the requested `energies` grid extends
    towards zero, then the soft mode energies may set a vanishing scale for the
    required energy resolution ``w``. To reduce computational cost, one may consider
    lifting the `energies` window or introducing an intrinsic broadening ``η > 0``
    to damp soft-mode intensity. Otherwise, increasing the `maxiters` budget may be
    necessary. Soft modes (or instabilities of the magnetic order) will be detected
    on a best-effort basis. In particular, spectral features that are significantly
    below the `kernel` width may be missed.

!!! tip "Consider `SampledCorrelations` when calculating powder averages"

    Spin wave theory performs an independent calculation for each ``𝐪`` point of
    interest. Consequently, it can be very slow to sample a 3D volume of
    ``𝐪``-space, e.g., as required for a [`powder_average`](@ref). Here a
    compelling alternative may be [`SampledCorrelations`](@ref). It uses real-time
    spin dynamics to calculate structure factor data over the entire 3D grid of
    commensurate ``𝐪``-vectors in one shot. This may provide a considerable speedup
    at the cost of limited ``𝐪``-space resolution and statistical sampling error.

This calculator uses the Lanczos method [1], which exhibits a near-optimal
convergence rate in the number of iterations. Lanczos may be viewed as a
refinement of the kernel polynomial method, which was previously implemented as
`SpinWaveTheoryKPM` [2].

## References

1. [T. Chen, _The Lanczos algorithm for matrix functions: a handbook for
   scientists_ (2024) [arXiv:2410.11090]](https://arxiv.org/abs/2410.11090).
2. [H. Lane et al., _Kernel Polynomial Method for Linear Spin Wave Theory_
   (2023) [arXiv:2312.08349]](https://arxiv.org/abs/2312.08349).
"""
struct SpinWaveTheoryLarge
    swt :: SpinWaveTheory
    tol :: Float64
    maxiters :: Int

    SpinWaveTheoryLarge(swt::SpinWaveTheory, tol, maxiters) = new(swt, tol, maxiters)
end

function SpinWaveTheoryLarge(sys::System; measure::Union{Nothing, MeasureSpec}, regularization=1e-8, tol=1e-3, maxiters=1000)
    0 < tol < 1 || error("Require 0 < tol < 1")
    maxiters > 0 || error("Require maxiters > 0")
    return SpinWaveTheoryLarge(SpinWaveTheory(sys; measure, regularization), tol, maxiters)
end

function Base.show(io::IO, ::MIME"text/plain", swt_large::SpinWaveTheoryLarge)
    (; swt) = swt_large
    printstyled(io, "SpinWaveTheoryLarge ", mode_to_str(swt.sys), "\n"; bold=true, color=:underline)
    println(io, "  ", nsites(swt.sys), " sites")
end

function intensities(swt_large::SpinWaveTheoryLarge, qpts; energies, η=0.0, kernel=nothing, kT=0.0, threaded=false, verbose=false)
    qpts = convert(AbstractQPoints, qpts)
    η >= 0 || error("Require η ≥ 0")
    if iszero(η)
        # The limit η → 0⁺, approximated by η = √ε × kernel width, for machine
        # precision ε. A Lorentzian of half width η biases intensities by about
        # η/width relative to the peak, here √ε, below any attainable `tol`. The
        # floor is needed because, at kT > 0, the pole of 1 + n_B at ω = 0 is
        # evaluated at z = iη. There, χ(z) = -(c₀ + β₀² x₁)/z loses about ε
        # |c₀|/η to cancellation, i.e., √ε |c₀|/width, which is negligible.
        isnothing(kernel) && error("Specify a `kernel`, a broadening `η`, or both.")
        η = sqrt(eps()) * minimum(ω -> kernel_fwhm(kernel, ω), energies)
    end

    # Lanczos serves a real-valued measure. Any other measure is evaluated as
    # one real-valued measure per real component of its value. If this becomes a
    # bottleneck, runs could be shared instead. For r independent observables,
    # the matrices M of all components are Hermitian r × r, so r² runs from
    # fixed directions (eₐ, eₐ + e_b, and eₐ + i e_b) determine every component
    # by polarization identities, still without probes.
    (; swt, tol, maxiters) = swt_large
    parts = map(real_component_measures(swt.measure)) do measure
        swt_m = SpinWaveTheory(swt.sys, measure, swt.regularization, swt.vacuum, swt.classical_energy,
                               swt.quadratic_table, swt.bare_table, swt.ewald, swt.monomials, swt.observable_coefs)
        intensities_lanczos(SpinWaveTheoryLarge(swt_m, tol, maxiters), qpts; energies, η, kernel, kT, threaded, verbose)
    end
    T = eltype(swt.measure)
    data = map((xs...) -> from_real_components(T, xs), parts...)
    return Intensities(orig_crystal(swt.sys), qpts, collect(Float64, energies), data)
end

# At each wavevector, Lanczos builds the tridiagonal matrix T that represents ĨH
# in the H inner product, starting from v = Ĩu, for a combination u of the
# observable vectors. The susceptibility of u at complex energy z follows from a
# partial-fraction identity,
#
#     χ(z) = Σₖ |wₖ|² / (λₖ (λₖ - z)) = -(c₀ + β₀² [(z - T)⁻¹]₁₁) / z,
#
# where the sum runs over eigenpairs of ĨH with weights wₖ, β₀² = v† H v, and
# the commutator c₀ = Σₖ |wₖ|²/λₖ = v† Ĩ v is known exactly. Lanczos
# approximates only the resolvent, which is analytic within a distance Im z of
# the real axis. Pairs of exact zero modes at ±λ, e.g., of a finite magnetic
# cell with continuous symmetry, cancel inside χ″ as λ → 0.
#
# A real-valued intensity is tr(M C) for one Hermitian matrix M acting on the
# correlations C of the observables; see `measure_forms`. With the eigenpairs
# M = Σⱼ μⱼ eⱼ eⱼ†, it is Σⱼ μⱼ eⱼ† C eⱼ, a sum of the autocorrelations of the
# combinations u = Σₐ eⱼₐ uₐ, so that one run per nonzero eigenvalue suffices,
# e.g., 2 runs for `ssf_perp`.
#
# 1. T. Chen, The Lanczos algorithm for matrix functions: a handbook for
#    scientists (2024), Sec. 5.
function intensities_lanczos(swt_large, qpts; energies, η, kernel, kT, threaded, verbose)
    (; swt, tol, maxiters) = swt_large
    (; sys) = swt
    cryst = orig_crystal(sys)
    @assert eltype(swt.measure) == Float64

    Ncells = nsites(sys) / natoms(cryst)
    n = 2nbands(swt)
    Nobs = num_observables(swt.measure)
    rb = RetardedBroadening(energies; η, kernel, kT)
    Λ = spectral_scale(swt)
    # Lanczos steps of each run per wavevector, for the `verbose` report
    iters = [Int[] for _ in qpts.qs]
    data = zeros(length(energies), size(qpts.qs)...)

    # The closures read the abstract field `swt.measure` rather than capture
    # the measure, so that they compile once for every combiner
    newbuf() = (
        u = zeros(ComplexF64, n, Nobs),
        H = HamiltonianOperator(swt),
    )

    function calc_iq!(buf, iq)
        (; u, H) = buf
        q = qpts.qs[iq]
        q_global = cryst.recipvecs * q
        observable_vectors!(u, swt, q_global)
        set_wavevector!(H, to_supercell_rlu(sys, q))
        mulS!(y, x) = mul!(y, H, x)

        # A run from Ĩ u e for each eigenvector e of M with nonzero eigenvalue,
        # written by broadcast, which avoids BLAS
        (μs, es) = eigen(only(measure_forms(swt.measure, q_global)))
        ks = findall(μ -> abs(μ) > 1e-12 * maximum(abs, μs), μs)
        (μs, es) = (μs[ks], es[:, ks])
        V = mul_Ĩ!(similar(u), u)
        starts = [sum(e[a] .* view(V, :, a) for a in 1:Nobs) for e in eachcol(es)]

        (As, nrows) = try
            lanczos_correlations(mulS!, starts, rb; Λ, tol, maxiters)
        catch err
            err isa InstabilityError || rethrow()
            throw(InstabilityError("Not an energy-minimum; wavevector q = $(vec3_to_string(q)) unstable."))
        end
        if isnothing(As)
            error("Exceeded maxiters = $maxiters at wavevector q = $(vec3_to_string(q)); \
                   consider increasing η or maxiters")
        end
        iters[iq] = nrows

        # The structure factor per chemical cell, Σⱼ μⱼ Im χⱼ / π
        for iω in axes(data, 1)
            data[iω, iq] = sum(μs[j] * imag(As[j][iω]) for j in eachindex(As); init=0.0) / (π*Ncells)
        end
    end

    verbose && println("SpinWaveTheoryLarge intensities\n  Supercell       $(supercell_to_str(sys.dims, sys.crystal)), $(n ÷ 2) bosons")
    desc = verbose ? "  Wavevectors     " : nothing
    elapsed = @elapsed foreach_chunked(calc_iq!, newbuf, CartesianIndices(qpts.qs); threaded, desc)

    if verbose
        nthreads = threaded ? min(Threads.nthreads(), length(qpts.qs)) : 1
        per_q = number_to_simple_string(1000 * elapsed / length(qpts.qs); digits=2)
        println("  Elapsed         $(number_to_simple_string(elapsed; digits=2)) s on $nthreads \
                 thread$(nthreads == 1 ? "" : "s"), $per_q ms per 𝐪")
        (lo, hi) = extrema(length, iters)
        println("  Lanczos runs    $(lo == hi ? lo : "$lo to $hi") per 𝐪")
        steps = reduce(vcat, iters)
        if !isempty(steps)
            mean = number_to_simple_string(Statistics.mean(steps); digits=3)
            median = number_to_simple_string(Statistics.median(steps); digits=3)
            println("  Lanczos steps   mean $mean, median $median, max $(maximum(steps)) per run")
        end
    end
    return data
end

# The broadened χ(z), as Σⱼ Wᵢⱼ χ(zⱼ) at the energies of `rb`, from a Lanczos
# run for each of the starting vectors `vs`, as As[k] for run k. Starts that are
# negligible relative to the largest, e.g., that vanish at this wavevector, give
# zero, rather than iterating on roundoff. Returns also the number of Lanczos
# steps taken by each run, or `nothing` for As if not converged within
# `maxiters` steps.
#
# Iteration stops when every error estimate of `resolvent_and_error!` is below
# `tol/2` relative to the largest intensity Im χ / π of the runs, a factor 2
# margin for the estimates. The runs advance together, because this largest
# intensity is itself an estimate, which can fall by orders of magnitude while
# some run is unconverged, e.g., at a soft mode as η → 0.
#
# Checks are at even numbers of rows of T. For a real combination of Hermitian
# observables, the weights at ±λ are equal, so that α_j = 0, and the spectrum of
# T is symmetric under θ → -θ. An odd number of rows then forces a spurious Ritz
# value at zero. The steps between checks grow by 25%, so that the overhead of
# checks is a fixed fraction.
function lanczos_correlations(mulS!, vs, rb; Λ, tol, maxiters)
    As = [zeros(ComplexF64, size(rb.W, 1)) for _ in vs]
    norms = [sqrt(real(dot_noblas(v, v))) for v in vs]
    ks = findall(>(sqrt(eps()) * maximum(norms; init=0.0)), norms)
    rs = [LanczosResolvent(mulS!, vs[k], rb.zs, Λ) for k in ks]
    isempty(rs) && return (As, Int[])
    errs = zeros(length(rs))
    m = 10
    while true
        for (k, r) in enumerate(rs)
            extend!(r, m)
            errs[k] = resolvent_and_error!(As[ks[k]], r, rb)
        end
        scale = maximum(A -> maximum(x -> abs(imag(x)), A), As) / π
        nrows = [r.nrows for r in rs]
        2 * maximum(errs) <= tol * scale && return (As, nrows)
        m >= maxiters && return (nothing, nrows)
        m += max(2, m ÷ 4)
        m += isodd(m)
    end
end

# Writes into A the broadened χ = -(c₀ + β₀² x₁) / z of the start. Returns the
# largest broadened error estimate.
#
# The estimate is the change of the intensity Im A / π by two corrections to χ
# of the start, broadened with their phases:
# 1. The correction δx₁ of `error_estimate`, for the rows of T not yet computed.
#    The corrections at nearby nodes largely cancel under the kernel, and their
#    magnitudes would overestimate by orders of magnitude.
# 2. Spectral weight near ω = 0 that T has not resolved. The residue of a mode
#    in χ″ is its weight divided by its energy λ, and the thermal factor 1 + n_B
#    adds about kT/λ, so that a soft mode can be bright although it has
#    negligible weight in the moments that T matches. Here the largest such
#    weight ν = β₀² / Σⱼ pⱼ(0)², by the Christoffel function (see
#    `LanczosResolvent`), is a pair at ±δ, for the smallest energy δ that T
#    resolves, about πβ/m for m rows and Ritz values in [-2β, 2β], with χ = ν /
#    (δ² - z²). Its broadening by the same weights includes the thermal factor.
# These are a best effort. The intensity of a mode far below δ, with weight much
# less than ν, may be missed. A broadening η suppresses the intensity of modes
# below η, which changes the target itself, but reproducibly.
function resolvent_and_error!(A, r::LanczosResolvent, rb::RetardedBroadening)
    (; zs, W, ranges) = rb
    (; β₀, c₀, x₁) = r
    fill!(A, 0)
    err = zeros(ComplexF64, length(A))
    ν = exhausted(r) ? 0.0 : β₀^2 / r.p₀sum
    δ = π * r.βmax / max(r.nrows, 1)
    @inbounds for j in eachindex(zs)
        χ = -(c₀ + β₀^2 * x₁[j]) / zs[j]
        δχ = -(β₀^2 * error_estimate(r, j) + ν * zs[j] / (zs[j]^2 - δ^2)) / zs[j]
        for i in ranges[j]
            A[i] += W[i, j] * χ
            err[i] += W[i, j] * δχ
        end
    end
    return maximum(x -> abs(imag(x)), err; init=0.0) / π
end
