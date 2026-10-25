# Dynamical structure factor at O(1/s). The default `:dyson` is the exact
# resolvent χ of the auxiliary quadratic model in Corrections.jl, with  
# S = (χ' - χ)/2πi at z = ω + iη.
#
# The alternative `:particle` swaps in the propagator of Eq. (12) of
# arXiv:1306.1231: the Dyson equation projected onto the particle block, with
# the source channel frozen on shell. It then assembles the same χ, interference
# included, minus the O(1/s²) term K_dm G K_md. It is not Goldstone-protected,
# and its O(1/s) spectral weights are inexact. The published figures omit the
# interference, so reproducing them means summing `transverse + direct`.

"""
    corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                          mean_field_maxevals=100_000, threaded=false, verbose=false)

Dynamical spin structure factor at temperature ``T = 0``, including corrections
of relative order ``1/s`` to linear spin wave theory. Magnon energies are
shifted by the mean fields and the cubic self-energy, magnons that can decay
into two magnons are broadened, and the two-magnon continuum is included, both
the part fed by decay and the part the observable creates directly.

The regulator `η` is required. The Green function is evaluated at ``ω + iη``,
which is the same as convolving the spectrum with `lorentzian(fwhm=2η)`. Choose
it small compared to the linewidths of interest, but no smaller: the wavevector
grid of the loop integrals grows as `1/η` in each dispersing dimension.
Instrumental resolution should be applied to the result afterwards.

The accuracy target `tol` sets the loop grid (unless `loop_grid` is given) and
the adaptive cubature of the static mean fields, which stops after
`mean_field_maxevals` evaluations.

The result is exact at order ``1/s``, preserves Goldstone modes, and satisfies
the frequency sum rule. As a retarded response it includes the tails of the
negative-frequency poles. A warning is issued at wavevectors where the ``1/s``
correction exceeds the harmonic energy, so that the expansion is not controlled.
There the spectrum may show spurious low-energy bands.

Set `threaded=true` to parallelize over `qpts`, and `verbose=true` to print
parameters, progress and linewidths.
"""
function corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, kwargs...)
    (; cryst, qpts, energies, transverse, cross, direct) = corrected_channels(swt, qpts; energies, η, kwargs...)
    data = transverse + cross + direct
    return Intensities(cryst, qpts, energies, reshape(data, length(energies), size(qpts.qs)...))
end

# Workhorse of `corrected_intensities`, returning the three terms of S
# separately as (energy × wavevector) matrices: `transverse` from w'Gw, `cross`
# from the terms linear in K_md, and `direct` from the rest. Also returns
# `disp`, the harmonic energies; `stability`, described below; and if
# `spectral=true` then `specfunc`, the particle block of the magnon spectral
# matrix (G' - G)/2πi.
function corrected_channels(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                            mean_field_maxevals=100_000, resummation=:dyson, threaded=false,
                            verbose=false, spectral=false)
    check_corrections_supported(swt)
    η > 0 || error("Regulator `η` must be positive.")
    resummation in (:dyson, :particle) || error("Unknown resummation `$resummation`; use :dyson or :particle.")

    (; sys, measure) = swt
    cryst = orig_crystal(sys)
    L = nbands(swt)
    Nobs = num_observables(measure)
    Ncells = nsites(sys) / natoms(cryst)

    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")
    qpts = convert(AbstractQPoints, qpts)

    loop_grid = @something loop_grid auto_loop_grid(swt, η, tol)
    # Binning error is O((Δ/η)²), so Δ/η = √tol contributes of order `tol`
    bin_width = η * min(1/2, √tol)

    if length(energies) > 1
        dω = (energies[end] - energies[begin]) / (length(energies) - 1)
        if dω > (η/2) * (1 + 1e-8)
            @warn """Requested `energies` are spaced by $(round(dω, sigdigits=2)) on \
                     average, which will not resolve features of width η = $η. A spacing \
                     of η/2 or less adds little cost."""
        end
    end

    # A negligible cubic vertex is dropped outright, e.g. collinear order in
    # dipole mode, so that everything it generates is exactly zero
    terms3 = cubic_monomials(swt)
    cubic = !cubic_vertex_vanishes(swt, terms3)
    cubic || empty!(terms3)

    tad = cubic ? tadpole_correction(swt; tol, maxevals=mean_field_maxevals) : nothing
    terms2 = [hartree_fock_correction(swt; maxiters=1, tol, maxevals=mean_field_maxevals).terms2
              isnothing(tad) ? BosonMonomial{2}[] : tad.terms2
              anisotropy_correction(swt).terms2]
    δc = observable_corrections(swt; v = isnothing(tad) ? nothing : tad.v,
                                tol, maxevals=mean_field_maxevals)

    chans = (; transverse = zeros(eltype(measure), length(energies), length(qpts.qs)),
               cross = zeros(eltype(measure), length(energies), length(qpts.qs)),
               direct = zeros(eltype(measure), length(energies), length(qpts.qs)))
    specfunc = spectral ? zeros(ComplexF64, L, L, length(energies), length(qpts.qs)) : nothing
    disp = zeros(L, length(qpts.qs))
    # On-shell decay rate of each band, for the `verbose` report
    linewidths = fill(NaN, L, length(qpts.qs))
    # Least eigenvalue of |ε| + Σstat + K_mm(0), relative to |ε| + Σstat. It is
    # negative where the auxiliary model is unstable, i.e., the 1/s expansion
    # breaks down. The normalization matters near a Goldstone mode, where Σstat
    # and the loop each grow like 1/ε and cancel.
    stability = zeros(length(qpts.qs))

    # Propagated Nambu indices, and the rows of K for the direct amplitudes
    p = resummation == :dyson ? (1:2L) : (1:L)
    d = 2L .+ (1:Nobs)
    Ĩ = Diagonal([ones(L); -ones(L)][p])

    function calc_iq!(iq)
        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Nobs)
        δH = zeros(ComplexF64, 2L, 2L)
        corr = zeros(ComplexF64, num_correlations(measure))

        q = qpts.qs[iq]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        ε = excitations!(T, H, swt, q)
        view(disp, :, iq) .= view(ε, 1:L)
        E = Diagonal(abs.(ε[p]))

        accum_quadratic!(δH, terms2, q_reshaped)
        Σstat = (T' * δH * T)[p, p]

        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
        accum_observable_corrections!(u, swt, q_reshaped, q_global, δc)
        w = (T' * u)[p, :]

        grid = loop_wavevectors(loop_grid, q_reshaped)
        words2 = observable_pair_words(swt, q_reshaped, q_global)
        (; decay, source) = pair_measures(swt, terms3, q_reshaped, grid; bin_width, words2)
        # The final frequency, just above ω = 0, serves the stability check
        zs = [energies .+ im*η; im*η]
        if resummation == :dyson
            K = cauchy_transform((decay, source), zs)
        else
            # Source channel frozen at the mean on-shell energy of the two legs
            K = cauchy_transform((decay,), zs)
            Σstat += [sum(((j, v),) -> bin_entry(v, m, m′) / ((ε[m] + ε[m′])/2 - j * bin_width), source.bins; init=0im)
                      for m in p, m′ in p]
        end

        Σ0 = E + Σstat + K[p, p, end]
        stability[iq] = eigmin(Hermitian((Σ0 + Σ0') / 2)) / opnorm(E + Σstat)

        # Contracts an Nobs×Nobs χ through the measure, taking S = (χ' - χ)/2πi
        function accum_channel!(accum, iω, χ)
            map!(((μ, ν),) -> (conj(χ[ν, μ]) - χ[μ, ν]) / (2π*im*Ncells), corr, measure.corr_pairs)
            accum[iω, iq] = measure.combiner(q_global, corr)
        end

        for iω in eachindex(energies)
            Kz = view(K, :, :, iω)
            (Kmm, Kmd, Kdm, Kdd) = (Kz[p, p], Kz[p, d], Kz[d, p], Kz[d, d])
            G = inv(zs[iω]*Ĩ - E - Σstat - Kmm)
            accum_channel!(chans.transverse, iω, w' * G * w)
            accum_channel!(chans.cross, iω, Kdm * G * w + w' * G * Kmd)
            # Pair → magnon → pair is of order 1/s², and published schemes omit it
            accum_channel!(chans.direct, iω, resummation == :dyson ? Kdd + Kdm * G * Kmd : Kdd)
            isnothing(specfunc) || (view(specfunc, :, :, iω, iq) .= ((G' - G) / (2π*im))[1:L, 1:L])
        end

        for n in 1:L
            if energies[begin] ≤ ε[n] ≤ energies[end]
                iω = argmin(iω -> abs(energies[iω] - ε[n]), eachindex(energies))
                linewidths[n, iq] = -imag(K[n, n, iω])
            end
        end
    end

    if verbose
        println("""
            corrected_intensities with tol = $tol
              loop grid       $(join(loop_grid, "×")) = $(prod(loop_grid)) points""")
    end

    t0 = time()
    foreach_maybe_threaded(calc_iq!, threaded, eachindex(qpts.qs);
                           desc = verbose ? "  wavevectors     " : nothing)
    elapsed = time() - t0

    unstable = count(<(-tol), stability)
    if resummation == :dyson && unstable > 0
        @warn """The 1/s expansion is not controlled at $unstable of $(length(qpts.qs)) \
                 wavevectors, where the correction exceeds the harmonic energy (relative \
                 least eigenvalue of ε + Σ(0) is $(round(minimum(stability), sigdigits=2))). \
                 Spurious low-energy bands may appear there."""
    end

    if verbose
        r2 = x -> round(x; sigdigits=2)
        nthreads = threaded ? min(Threads.nthreads(), length(qpts.qs)) : 1
        per_q = 1000 * elapsed / length(qpts.qs)
        Γs = sort!(filter(isfinite, vec(linewidths)))
        report = isempty(Γs) ? "none of the LSWT poles lie within `energies`" :
            "median $(r2(Γs[cld(end, 2)])), 90th pct $(r2(Γs[ceil(Int, 0.9end)])), against η = $(r2(η))"
        println("  elapsed         $(r2(elapsed)) s on $nthreads \
                 thread$(nthreads == 1 ? "" : "s"), $(r2(per_q)) ms per 𝐪")
        println("  on-shell -Im Σ  $report")
    end

    return (; cryst, qpts, energies, chans..., specfunc, disp, stability)
end
