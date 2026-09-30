# Dynamical structure factor at O(1/s). Conventions are in Corrections.jl. The
# default resummation, `:dyson`, is the exact resolvent of the auxiliary
# quadratic model described there,
#
#     X(z) = (w + K_dm)' G (w + K_md) + K_dd,   G = (zĨ - |ε| - Σstat - K_mm)⁻¹,
#
# over the full 2L Nambu space, evaluated at z = ω + iη. The structure factor is
# its anti-Hermitian part, S = (X' - X)/2πi. Both alternatives are for
# comparison only: `:perturbative` expands G to first order and drops the O(1/s²)
# term K_dm G K_md; `:particle` projects onto the particle block and freezes the
# source channel on shell, as in Eq. (12) of arXiv:1306.1231.

"""
    corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                          mean_field_maxevals=100_000, resummation=:dyson, threaded=false,
                          verbose=false)

Dynamical spin structure factor at temperature ``T = 0``, including corrections
of relative order ``1/s`` to linear spin wave theory. Magnon energies are shifted
by the mean fields and the cubic self-energy, magnons that can decay into two
magnons are broadened, and the two-magnon continuum is included, both the part
fed by decay and the part the observable creates directly.

The regulator `η` is required. The Green function is evaluated at ``ω + iη``,
which is the same as convolving the spectrum with `lorentzian(fwhm=2η)`. Choose
it small compared to the linewidths of interest, but no smaller: the wavevector
grid of the loop integrals grows as `1/η` in each dispersing dimension.
Instrumental resolution should be applied to the result afterwards.

The accuracy target `tol` sets the loop grid (unless `loop_grid` is given) and
the adaptive cubature of the static mean fields, which stops after
`mean_field_maxevals` evaluations.

The default `resummation=:dyson` solves the Dyson equation exactly in the full
Nambu space. It is exact at order ``1/s``, preserves Goldstone modes, and is
positive at ``ω > 0`` up to Lorentzian tails of the negative-frequency mirror
poles. A warning is issued where the resummed quadratic model is unstable, which
signals a breakdown of the ``1/s`` expansion. The alternatives `:perturbative`
and `:particle` are for testing and comparison with the literature.

Set `threaded=true` to parallelize over `qpts`, and `verbose=true` to print
parameters, progress and linewidths.
"""
function corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, kwargs...)
    (; cryst, qpts, energies, transverse, cross, direct) = corrected_channels(swt, qpts; energies, η, kwargs...)
    data = transverse + cross + direct
    return Intensities(cryst, qpts, energies, reshape(data, length(energies), size(qpts.qs)...))
end

# Workhorse of `corrected_intensities`, returning the three terms of S separately
# as (energy × wavevector) matrices: `transverse` from w' G w, `direct` from K_dd +
# K_dm G K_md, and `cross` from the terms linear in w. Published 1/s calculations
# omit `cross`, which cancels in a zone sum of a trace measure. Also returns
# `disp`, the harmonic energies, and if `spectral=true` then `specfunc`, the
# particle block of the magnon spectral matrix (G' - G)/2πi.
function corrected_channels(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                            mean_field_maxevals=100_000, resummation=:dyson, threaded=false,
                            verbose=false, spectral=false)
    check_corrections_supported(swt)
    η > 0 || error("Regulator `η` must be positive.")
    resummation in (:dyson, :perturbative, :particle) ||
        error("Unknown resummation `$resummation`; use :dyson, :perturbative or :particle.")

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
    # Least eigenvalue of |ε| + Σ(0), negative where the auxiliary model of :dyson
    # is unstable. It is relative to the norm of |ε| + Σstat, because near a
    # Goldstone mode both Σstat and the loop grow like 1/ε and cancel.
    stability = zeros(length(qpts.qs))

    # Nambu indices that are propagated, and the remaining rows of K
    p = resummation == :particle ? (1:L) : (1:2L)
    d = 2L .+ (1:Nobs)

    function calc_iq!(iq)
        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Nobs)
        corr = zeros(ComplexF64, num_correlations(measure))

        q = qpts.qs[iq]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        ε = excitations!(T, H, swt, q)
        view(disp, :, iq) .= view(ε, 1:L)

        δH = zeros(ComplexF64, 2L, 2L)
        accum_quadratic!(δH, terms2, q_reshaped)
        Σstat = (T' * δH * T)[p, p]

        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
        accum_observable_corrections!(u, swt, q_reshaped, q_global, δc)
        w = (T' * u)[p, :]

        grid = loop_wavevectors(loop_grid, q_reshaped)
        words2 = observable_pair_words(swt, q_reshaped, q_global)
        (; decay, source) = pair_measures(swt, terms3, q_reshaped, grid; bin_width, words2)
        # The last frequency, just above ω = 0, is for the stability check
        zs = [energies .+ im*η; im*η]
        if resummation == :particle
            # Source channel frozen at the mean on-shell energy of the two legs
            K = cauchy_transform((decay,), zs)
            Σstat += [sum(((j, v),) -> bin_entry(v, m, m′) / ((ε[m] + ε[m′])/2 - j * bin_width), source.bins; init=0im)
                      for m in p, m′ in p]
        else
            K = cauchy_transform((decay, source), zs)
        end

        Ĩ = Diagonal([ones(L); -ones(L)][p])
        E = Diagonal(abs.(ε[p]))
        Σ0 = E + Σstat + K[p, p, end]
        stability[iq] = eigmin(Hermitian((Σ0 + Σ0') / 2)) / opnorm(E + Σstat)

        # Static sizes make the ω loop allocation-free; `Val` dispatch is the barrier
        Np = length(p)
        Kpd = K[[p; d], [p; d], eachindex(energies)]
        spec = isnothing(specfunc) ? nothing : view(specfunc, :, :, :, iq)
        dyson_loop!(chans, iq, spec, Val{Np}(), Val{Nobs}(), Kpd, zs, Ĩ, E, Σstat, w, resummation,
                    measure, q_global, corr, Ncells, L)

        for n in 1:L
            if energies[begin] ≤ ε[n] ≤ energies[end]
                iω = argmin(iω -> abs(energies[iω] - ε[n]), eachindex(energies))
                linewidths[n, iq] = -imag(K[n, n, iω])
            end
        end
    end

    if verbose
        println("""
            corrected_intensities with tol = $tol, resummation = $resummation
              loop grid       $(join(loop_grid, "×")) = $(prod(loop_grid)) points""")
    end

    t0 = time()
    foreach_maybe_threaded(calc_iq!, threaded, eachindex(qpts.qs);
                           desc = verbose ? "  wavevectors     " : nothing)
    elapsed = time() - t0

    unstable = count(<(-tol), stability)
    if resummation == :dyson && unstable > 0
        @warn """Resummed spin wave theory is unstable at $unstable of $(length(qpts.qs)) \
                 wavevectors (relative least eigenvalue of ε + Σ(0) is \
                 $(round(minimum(stability), sigdigits=2))). The 1/s expansion is not \
                 controlled there."""
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

# Frequency loop of `corrected_channels`, over the static blocks of K. Writes the
# three channels S = (X' - X)/2πi at wavevector index `iq`.
function dyson_loop!(chans, iq, spec, ::Val{Np}, ::Val{Nobs}, Kpd, zs, Ĩ, E, Σstat, w, resummation,
                     measure, q_global, corr, Ncells, L) where {Np, Nobs}
    Nk = Np + Nobs
    (Ĩ, E, Σstat, w) = (SMatrix{Np, Np}(Ĩ), SMatrix{Np, Np}(E), SMatrix{Np, Np, ComplexF64}(Σstat), SMatrix{Np, Nobs}(w))
    (ip, id) = (SVector{Np}(1:Np), SVector{Nobs}(Np+1:Nk))
    function accum!(accum, iω, X)
        map!(((μ, ν),) -> (conj(X[ν, μ]) - X[μ, ν]) / (2π*im*Ncells), corr, measure.corr_pairs)
        accum[iω, iq] = measure.combiner(q_global, corr)
    end
    for iω in axes(Kpd, 3)
        z = zs[iω]
        Kz = SMatrix{Nk, Nk}(view(Kpd, :, :, iω))
        (Kmm, Kmd, Kdm, Kdd) = (Kz[ip, ip], Kz[ip, id], Kz[id, ip], Kz[id, id])
        if resummation == :perturbative
            G0 = inv(z*Ĩ - E)
            G = G0 + G0 * (Σstat + Kmm) * G0
            Gc = G0
        else
            G = Gc = inv(z*Ĩ - E - Σstat - Kmm)
        end
        # Pair → magnon → pair, of order 1/s², is kept only by the exact resolvent
        Gd = resummation == :dyson ? Kdm * G * Kmd : zero(Kdd)
        accum!(chans.transverse, iω, w' * G * w)
        accum!(chans.cross, iω, Kdm * Gc * w + w' * Gc * Kmd)
        accum!(chans.direct, iω, Kdd + Gd)
        isnothing(spec) || (view(spec, :, :, iω) .= ((G' - G) / (2π*im))[1:L, 1:L])
    end
end
