# Dynamical structure factor at O(1/s), S = (χ' - χ)/2πi at z = ω + iη. Each
# `dyson` scheme propagates the magnons with the 2L×2L Nambu resolvent G of a
# quadratic model of magnons coupled to the binned two-magnon bath of
# Corrections.jl. The schemes differ in that model, and in how much of its
# response they assemble.
#
# `:nambu` is the full auxiliary model, exact at O(1/s). It assembles the whole
# Schur complement of the model's resolvent onto the observables,
#
#     χ = w'Gw + K_dm G w + w'G K_md + K_dm G K_md + K_dd,
#
# which inherits the frequency sum rule and is positive at ω > 0, up to the O(η)
# tails of the mirror poles, whenever the model is stable.
#
# The reproduction schemes keep the published assembly, the magnon term w'Gw
# plus the bare continuum K_dd, each positive up to the same tails while the
# model is stable.
#
# `:particle` is the rotating-wave truncation of the auxiliary model, the
# canonical treatment of Chernyshev and Zhitomirsky, PRB 79, 144416 (2009) and
# RMP 85, 219 (2013). Each bath couples only to the legs it resonates with,
# decay pairs to particles and source pairs to holes, and the anomalous blocks
# of the self-energy are dropped, which the papers justify for the poles at
# O(1/s) (PRB Sec. IV). Each block keeps its non-resonant channel frozen on
# shell, as in PRB 88, 094407 (2013), which removes a spurious branch pushed up
# from negative frequency near ±Q. The particle block is then exactly G₁₁ of
# RMP Eq. (34), PRB Eq. (86), and the hole block its mirror G₁₁(-𝐪, -ω). The
# dropped counter-rotating terms are suppressed by |Σ_anom|/(ω + ε), so they are
# O(1/s) in the intensities and not small near a soft mode whose rotation
# generator has a quadratic part (e.g. at ±K of the triangular lattice), where
# the 1/ε divergences cancel only in the full Nambu inverse. There the scheme
# overestimates the intensity at ω ≫ ε.
#
# `:on_shell` freezes the resonant channel too, at ε + iη, so that each band
# becomes one complex pole, the first-order result of RMP Eq. (35), PRB Eq.
# (55). Each element of the frozen matrix is evaluated at the mean energy of
# its two legs, reducing to the published diagonal form when bands are well
# separated. The pole has logarithmic singularities where a band crosses a
# saddle point of the continuum (PRB Sec. V). Near a Goldstone mode, ε̃ and Γ
# are a residue of cancelling 1/ε terms and are sensitive to `η` and to the loop
# grid once ε falls below a few η.
#
# No scheme dresses the internal lines: RMP Sec. IV.C notes that self-consistent
# dressing gaps the Goldstone modes.

"""
    corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, dyson=:nambu,
                          mark_uncontrolled=true, threaded=false, verbose=false)

Dynamical spin structure factor at temperature ``T = 0``, including one-loop
spin-wave corrections, i.e., order ``1/s`` in dipole mode. Magnon energies are
shifted by the mean fields and the cubic self-energy, magnons that can decay
into two magnons are broadened, and the two-magnon continuum is included, both
the part fed by decay and the part the observable creates directly.

The regulator `η` is required. The Green function is evaluated at ``ω + iη``,
which is the same as convolving the spectrum with `lorentzian(fwhm=2η)`. Choose
it small compared to the linewidths of interest, but no smaller: the wavevector
grid of the loop integrals grows as `1/η` in each dispersing dimension.
Instrumental resolution should be applied to the result afterwards.

The accuracy target `tol` sets the wavevector grid of the loop integrals and the
adaptive cubature of the static mean fields.

The `dyson` option selects how the one-loop self-energy enters the Dyson
equation for the magnon propagator.

- `:nambu` (default) solves the Dyson equation in the full particle/hole (Nambu)
  space. It includes all terms at order ``1/s``, preserves Goldstone modes at
  all frequencies, and satisfies the frequency sum rule.
- `:particle` drops the anomalous (particle/hole mixing) self-energy, following
  Eq. (86) of [Chernyshev and Zhitomirsky, PRB **79**, 144416
  (2009)](https://doi.org/10.1103/PhysRevB.79.144416). The magnon poles are
  correct at order ``1/s``, but near some Goldstone modes it produces spurious
  intensity well above the magnon energy.
- `:on_shell` evaluates the self-energy at the harmonic magnon energies, so that
  each magnon becomes a Lorentzian of shifted energy and finite width, following
  Eq. (55) of the same reference. Peak energies have logarithmic singularities
  where a magnon crosses a saddle point of the two-magnon continuum, and are
  sensitive to `η` near a Goldstone mode.

In the latter two schemes the intensity is the magnon term plus the bare
two-magnon continuum. As a retarded response, the result of each scheme includes
the tails of the negative-frequency poles, which cancel the Lorentzian tail of a
Goldstone mode at high energy.

In each scheme the renormalized magnon poles may move anywhere, including to
imaginary frequency. If any pole does so, or falls below half of the lowest
harmonic magnon energy at its wavevector ``𝐪``, then Sunny interprets this as
an _uncontrolled_ application of the one-loop correction for that ``𝐪``. This
breakdown is marked by setting all intensity within `±η` of the uncontrolled
pole to `NaN`. Set `mark_uncontrolled=false` to retain all raw intensity data.

Set `threaded=true` to parallelize over `qpts`, and `verbose=true` to print a
progress bar and other diagnostics.
"""
function corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                               mean_field_maxevals=100_000, mark_uncontrolled=true, dyson=:nambu,
                               threaded=false, verbose=false)
    (; cryst, qpts, energies, transverse, cross, direct, artifacts) =
        corrected_channels(swt, qpts; energies, η, tol, loop_grid, mean_field_maxevals, dyson,
                           threaded, verbose)
    data = transverse + cross + direct
    mark_uncontrolled && (data[artifacts] .= NaN)
    return Intensities(cryst, qpts, energies, reshape(data, length(energies), size(qpts.qs)...))
end

"""
    corrected_intensities_bands(swt::SpinWaveTheory, qpts; η, tol=0.01, threaded=false,
                                verbose=false)

Magnon bands at temperature ``T = 0`` with one-loop corrections, i.e., order
``1/s`` in dipole mode, in the on-shell approximation. This is the
`dyson=:on_shell` option of [`corrected_intensities`](@ref), keeping only
the magnon poles. The two-magnon continuum is omitted.

Each band carries a shifted energy, its intensity, and a half width at half
maximum in the field `widths`, arising from decay into the two-magnon continuum.
Here `η` regularizes the loop integrals only, and does not broaden the result.
"""
function corrected_intensities_bands(swt::SpinWaveTheory, qpts; η, tol=0.01, loop_grid=nothing,
                                     mean_field_maxevals=100_000, threaded=false, verbose=false)
    (; cryst, qpts, bands) = corrected_channels(swt, qpts; energies=Float64[], η, tol, loop_grid,
                                                mean_field_maxevals, dyson=:on_shell, threaded, verbose)
    sz = (size(bands.disp, 1), size(qpts.qs)...)
    return BandIntensities(cryst, qpts, reshape(bands.disp, sz), reshape(bands.data, sz), reshape(bands.widths, sz))
end

# Workhorse of `corrected_intensities`, returning the three terms of S
# separately as (energy × wavevector) matrices: `transverse` from w'Gw, `cross`
# from the terms linear in K_md, and `direct` from the rest. Also returns
# `disp`, the harmonic energies; `artifacts`, a mask over (energy, wavevector)
# of the uncontrolled poles described above; for `:on_shell`, `bands`, the
# energies, half widths and intensities of the poles; and if `spectral=true`
# then `specfunc`, the particle block of the magnon spectral matrix (G' - G)/2πi.
function corrected_channels(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                            mean_field_maxevals=100_000, dyson=:nambu, threaded=false,
                            verbose=false, spectral=false)
    check_corrections_supported(swt)
    η > 0 || error("Regulator `η` must be positive.")
    dyson in (:nambu, :particle, :on_shell) ||
        error("Unknown `dyson=:$dyson`; use :nambu, :particle or :on_shell.")

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
    # Poles of `:on_shell`: energies, half widths and intensities
    bands = dyson != :on_shell ? nothing :
        (; disp = zeros(L, length(qpts.qs)), widths = zeros(L, length(qpts.qs)),
           data = zeros(eltype(measure), L, length(qpts.qs)))
    artifacts = falses(length(energies), length(qpts.qs))

    # Positive real mode frequencies of Ĩ·Re M, the linearized problem at fixed
    # M, and whether any mode is unstable, i.e. has a complex frequency. The
    # resummed poles are where a mode frequency of M(ω) = |ε| + Σstat + K_mm(ω)
    # crosses ω.
    function mode_frequencies(M)
        zs = eigvals(Ĩ * (M + M') / 2)
        real_zs = filter(z -> abs(imag(z)) ≤ 1e-8 * opnorm(M), zs)
        return (sort!([real(z) for z in real_zs if real(z) > 0]), length(real_zs) < length(zs))
    end

    # Nambu indices of the magnon legs, and the rows of K for the direct amplitudes
    p = 1:2L
    d = 2L .+ (1:Nobs)
    Ĩ = Diagonal([ones(L); -ones(L)])
    # Particle and hole blocks of a 2L×2L matrix
    Pp = [m ≤ L && m′ ≤ L for m in p, m′ in p]
    Ph = [m > L && m′ > L for m in p, m′ in p]

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
        E = Diagonal(abs.(ε))

        accum_quadratic!(δH, terms2, q_reshaped)
        Σstat = T' * δH * T

        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
        accum_observable_corrections!(u, swt, q_reshaped, q_global, δc)
        w = T' * u

        grid = loop_wavevectors(loop_grid, q_reshaped)
        words2 = observable_pair_words(swt, q_reshaped, q_global)
        (; decay, source) = pair_measures(swt, terms3, q_reshaped, grid; bin_width, words2)
        # After `energies` come a frequency just above ω = 0, for the
        # instability check, and the harmonic band energies, for the on-shell
        # linewidths that `verbose` reports
        nω = length(energies)
        zs = [energies; 0; ε[1:L]] .+ im*η
        Kdec = cauchy_transform((decay,), zs)
        Ksrc = cauchy_transform((source,), zs)

        # Channel ρ frozen at the mean on-shell energy of the two legs, which
        # is negative in the hole block
        frozen(ρ, δ) = [sum(((j, v),) -> bin_entry(v, m, m′) / ((ε[m] + ε[m′])/2 + im*δ - j * bin_width), ρ.bins; init=0im)
                        for m in p, m′ in p]
        if dyson != :nambu
            # The rotating-wave truncation of the auxiliary model. Each bath
            # couples only to the legs it resonates with, decay pairs to
            # particles and source pairs to holes, the anomalous blocks are
            # dropped, and each block keeps its non-resonant channel frozen on
            # shell, which is Hermitian.
            Kdec[L+1:2L, :, :] .= 0
            Kdec[:, L+1:2L, :] .= 0
            Ksrc[1:L, :, :] .= 0
            Ksrc[:, 1:L, :] .= 0
            Σstat = Σstat .* (Pp .| Ph) + frozen(source, 0) .* Pp + frozen(decay, 0) .* Ph
            if dyson == :on_shell
                # The resonant channel frozen too, at ε + iη
                Σstat += frozen(decay, η) .* Pp + frozen(source, η) .* Ph
                Kdec[p, p, :] .= 0
                Ksrc[p, p, :] .= 0
            end
        end
        K = Kdec + Ksrc

        # Contracts an Nobs×Nobs χ through the measure, taking S = (χ' - χ)/2πi
        function accum_channel!(accum, iω, χ)
            map!(((μ, ν),) -> (conj(χ[ν, μ]) - χ[μ, ν]) / (2π*im*Ncells), corr, measure.corr_pairs)
            accum[iω, iq] = measure.combiner(q_global, corr)
        end

        for iω in eachindex(energies)
            Kz = view(K, :, :, iω)
            (Kmd, Kdm, Kdd) = (Kz[p, d], Kz[d, p], Kz[d, d])
            G = inv(zs[iω]*Ĩ - E - Σstat - Kz[p, p])
            accum_channel!(chans.transverse, iω, w' * G * w)
            # The routes through the bath complete the resolvent of the full
            # model. The reproduction schemes keep the published assembly
            # without them.
            if dyson == :nambu
                accum_channel!(chans.cross, iω, Kdm * G * w + w' * G * Kmd)
                accum_channel!(chans.direct, iω, Kdd + Kdm * G * Kmd)
            else
                accum_channel!(chans.direct, iω, Kdd)
            end
            isnothing(specfunc) || (view(specfunc, :, :, iω, iq) .= ((G' - G) / (2π*im))[1:L, 1:L])
        end

        # Mark the uncontrolled poles, those below half the lowest harmonic
        # energy. An unstable mode is a pole at ω = 0. Skip wavevectors whose
        # harmonic energy is below resolution, e.g. at a Goldstone mode, where
        # the static terms and the loop each grow like 1/ε and their residue
        # after cancellation is quadrature noise.
        εmin = minimum(abs, view(ε, 1:L))
        if εmin ≥ η
            mark!(ω) = (view(artifacts, :, iq) .|= abs.(energies .- ω) .≤ η)
            last(mode_frequencies(E + Σstat + K[p, p, nω+1])) && mark!(0.0)
            low = findall(<(εmin/2), energies)
            fs = [first(mode_frequencies(E + Σstat + K[p, p, iω])) .- energies[iω] for iω in low]
            for i in 1:length(low)-1
                for n in 1:min(length(fs[i]), length(fs[i+1]))
                    fs[i][n] * fs[i+1][n] < 0 && mark!((energies[low[i]] + energies[low[i+1]]) / 2)
                end
            end
        end

        if dyson == :on_shell
            # Poles ε̃ - iΓ of the particle block. Diagonalize its Hermitian
            # part; each pole takes its width from the anti-Hermitian part in
            # that basis.
            Σ = Σstat[1:L, 1:L]
            (ε̃, U) = eigen(Hermitian(E[1:L, 1:L] + (Σ + Σ') / 2))
            Γ = max.(0, imag.(diag(U' * (Σ' - Σ) * U)) / 2)
            amps = U' * w[1:L, :]
            for n in 1:L
                map!(((μ, ν),) -> conj(amps[n, μ]) * amps[n, ν] / Ncells, corr, measure.corr_pairs)
                (bands.disp[n, iq], bands.widths[n, iq]) = (ε̃[n], Γ[n])
                bands.data[n, iq] = measure.combiner(q_global, corr)
            end
        end

        for n in 1:L
            linewidths[n, iq] = -imag(Kdec[n, n, nω + 1 + n])
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

    # A wavevector with harmonic energy below resolution was skipped above. It
    # is marked unstable if any of its near neighbours in `qpts`, those within
    # twice the nearest distance, are uncontrolled.
    let
        εmins = vec(minimum(abs, disp; dims=1))
        resolved = findall(≥(η), εmins)
        ks = [cryst.recipvecs * q for q in qpts.qs]
        for iq in findall(<(η), εmins)
            isempty(resolved) && break
            ds = [norm(ks[j] - ks[iq]) for j in resolved]
            near = resolved[ds .≤ 2 * minimum(ds)]
            if any(j -> any(view(artifacts, :, j)), near)
                view(artifacts, :, iq) .|= abs.(energies) .≤ η
            end
        end
    end

    if verbose
        r2 = x -> round(x; sigdigits=2)
        nthreads = threaded ? min(Threads.nthreads(), length(qpts.qs)) : 1
        per_q = 1000 * elapsed / length(qpts.qs)
        Γs = sort!(filter(isfinite, vec(linewidths)))
        report = isempty(Γs) ? "none" :
            "median $(r2(Γs[cld(end, 2)])), 90th pct $(r2(Γs[ceil(Int, 0.9end)])), against η = $(r2(η))"
        println("  elapsed         $(r2(elapsed)) s on $nthreads \
                 thread$(nthreads == 1 ? "" : "s"), $(r2(per_q)) ms per 𝐪")
        println("  on-shell -Im Σ  $report")
        nbad = count(any, eachcol(artifacts))
        println("  uncontrolled    $nbad of $(length(qpts.qs)) wavevectors (NaN near poles below ε_min/2)")
    end

    return (; cryst, qpts, energies, chans..., specfunc, disp, artifacts, bands)
end
