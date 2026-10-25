# Dynamical structure factor at O(1/s), S = (χ' - χ)/2πi at z = ω + iη.
#
# The docstring of `corrected_intensities` states what each `scheme` gives and
# where it comes from; references [n] are to its list. These notes describe the
# construction and justify the claims made there.
#
# Each `scheme` propagates the magnons with the 2L×2L Nambu resolvent G of a
# quadratic model of magnons coupled to the binned two-magnon bath of
# Corrections.jl. The schemes make three independent choices: which blocks of
# the self-energy to keep, whether to invert the resolvent or take its
# first-order poles, and how much of the model's response to assemble.
#
# `:one_loop` keeps the whole model and assembles the Schur complement of its
# resolvent onto the observables,
#
#     χ = w'Gw + K_dm G w + w'G K_md + K_dm G K_md + K_dd,
#
# which inherits the frequency sum rule and is positive at ω > 0, up to the O(η)
# tails of the mirror poles, whenever the model is stable. `:two_magnon_ladder`
# is the same model with the pair interaction of H₄ added to the bath (see
# `ladder_transform`), which is solved exactly and so keeps all of this. Its
# static counterterm, which pins the Goldstone modes, removes a physical shift
# in a gapped magnet, so there it departs slightly from the exact ladder. Near a
# Goldstone mode in three dimensions its results are sensitive to the loop grid;
# in two dimensions it binds soft pairs below ω = 0 (see SelfEnergy.jl). The
# other two schemes assemble only the magnon term w'Gw and the bare continuum
# K_dd.
#
# `:one_loop_projected` is the rotating-wave truncation. Each bath couples only
# to the legs it resonates with, decay pairs to particles and source pairs to
# holes, and the anomalous blocks are dropped. Each block takes its non-resonant
# channel on shell, as in Ref. [6], which removes a spurious branch pushed up
# from negative frequency near ±Q. Off the diagonal, it takes the van Vleck form
# ½[K(ε_m) + K(ε_m′)] of `on_shell_form`. The particle block is then G₁₁ of Ref.
# [5], Eq. (34), and the hole block its mirror G₁₁(-𝐪, -ω). The dropped
# counter-rotating terms are suppressed by |Σ_anom|/(ω + ε), so they are O(1/s)
# in the intensities and not small near a soft mode whose rotation generator has
# a quadratic part (e.g. at ±K of the triangular lattice), where the 1/ε
# divergences cancel only in the full Nambu inverse. There the scheme
# overestimates the intensity at ω ≫ ε. Near Γ of a Heisenberg magnet,
# total-spin conservation bounds the integrated intensity of the trace measure,
# (∫S dω)² ≤ f(𝐪) χ/2 with first moment f ∝ q², so ∫S dω must vanish linearly.
# Under `:one_loop` it does, by a cancellation of the continuum against the
# interference. The projected bare continuum K_dd instead tends to a constant,
# as does the full response of the projected model, which is why the projection
# keeps the published assembly.
#
# `:one_loop_on_shell` linearizes the full Nambu Dyson equation at each pole, as
# in Ref. [5], Eq. (35). The residues, which do change at first order, are
# deliberately not corrected, keeping the scheme a statement about the
# dispersion. The details, including which channel sets each width, are at
# `dyson_model`.
#
# The projection and the linearization both need the anomalous coupling to be
# small against 2ε. Near a Goldstone mode in two dimensions it is not, which is
# the origin of the breakdown near Goldstone wavevectors described in the
# docstring; it extends a distance of order 1/s. Two soft internal lines feed a
# static bubble that diverges like 1/δ at distance δ, in the stiff direction
# directly and in the Goldstone direction through its O(δ) overlap with that
# one, so that the Goldstone block goes as δ(b s δ - c) instead of b s δ². The
# Ward identity protects only δ = 0. Inside δ ≲ c/(b s), the resummation of
# `:one_loop` has complex frequencies, and |Σ_anom|/2ε exceeds one, so
# `:one_loop_on_shell` there is a first-order formula outside its range.
#
# No scheme dresses the internal lines by default: Ref. [5], Sec. IV.C notes
# that self-consistent dressing gaps the Goldstone modes. A `SpinWaveTheory`
# whose Gaussian state carries a correction does so, with a counterterm.

"""
    corrected_intensities(swt::SpinWaveTheory, qpts, scheme; energies, η, kernel=nothing, kT=0,
                          bz=nothing, mark_breakdown=true, threaded=false, verbose=false)

Dynamical spin structure factor at temperature ``T = 0``, including one-loop
spin-wave corrections, i.e., relative order ``1/s`` in dipole mode. Magnon
energies are shifted by the mean fields and the cubic self-energy, magnons that
can decay into two magnons are broadened, and the response includes the
two-magnon continuum.

The Green function is evaluated at ``ω + iη``, so `η` both regulates the loop
integrals and broadens the result into Lorentzians. Choose it small compared to
the linewidths of interest, but no smaller: the wavevector grid of the loop
integrals grows as `1/η` in each dispersing dimension. An optional `kernel`,
e.g. for instrumental resolution, is applied to the result. Only `kT = 0` is
currently supported.

`bz` specifies a relative tolerance or an explicit grid for integration; see
[`BZIntegration`](@ref). The static mean fields are summed on the grid of the
loop integrals, which keeps Goldstone modes exactly gapless.

The bosons are expanded about the vacuum of `swt`, by default the harmonic
vacuum of LSWT; see [`to_hartree_fock_corrected_vacuum`](@ref). A vacuum from
[`to_self_consistent_vacuum`](@ref) is used on the grid on which it was solved.
Its correction to the quadratic Hamiltonian is subtracted as a counterterm, so
that the result differs from the ``1/s`` expansion only at the order neglected.
A correction that was added to pin a Goldstone mode is therefore removed again;
for the spectrum of the corrected quadratic Hamiltonian itself, use
[`intensities`](@ref).

The `scheme` selects how the one-loop self-energy is resummed:

- `:one_loop` is the recommended starting point. It is the complete one-loop
  response, including the interference between one- and two-magnon channels. It
  keeps Goldstone modes gapless and satisfies the sum rules at order ``1/s``.
- `:two_magnon_ladder` additionally lets the two magnons of each pair interact,
  to all orders. Use it for two-magnon bound states and resonances in gapped
  magnets. This scheme neglects mixing with states of three or more magnons, so
  it is most accurate for features that lie well below the three-magnon
  continuum. Ladder resummation is unreliable near a Goldstone mode, where it
  can bind pairs of soft magnons spuriously.
- `:one_loop_on_shell` replaces each magnon by a single Lorentzian, with the
  energy and width corrected at order ``1/s`` but the intensity of linear spin
  wave theory, and adds the bare two-magnon continuum. Use it to read off
  renormalized dispersions and decay rates.
- `:one_loop_projected` keeps only the particle block of the self-energy and
  neglects the interference between one- and two-magnon channels. The spectrum
  is then positive by construction, but near a Goldstone mode the intensity is
  overestimated. Use it mainly to compare with historical calculations.

A large one-loop correction can make the resummed propagator unstable, signaling
that perturbation theory has broken down at that wavevector. This is most common
in low dimensions, near a Goldstone wavevector of a state with cubic magnon
vertices. With `:one_loop` and `:two_magnon_ladder`, Sunny marks these
wavevectors by setting the intensity within `±η` of ``ω = 0`` to `NaN`; use
`mark_breakdown=false` to keep the raw data. The other two schemes cannot detect
a breakdown, but are no more reliable there.

The returned data is the structure factor at ``T = 0``, as for
[`intensities`](@ref). This is ``\\mathcal{S}(𝐪, ω) = θ(ω) χ''(𝐪, ω)``, where
the spectral function ``χ''(𝐪, ω) = -\\mathrm{Im}\\, χ(𝐪, ω + iη)/π`` of the
retarded response ``χ`` is odd, and ``θ(0) = 1/2``. Intensities vanish at ``ω <
0``. At ``ω > 0`` they include the Lorentzian tails of the negative-frequency
poles, which cancel the high-frequency tail of a Goldstone mode. The `kernel` is
applied after the thermal factor, so `energies` should be uniformly spaced, and
extend a kernel width beyond the range of interest.

Set `threaded=true` to parallelize over `qpts`, and `verbose=true` to print a
progress bar and diagnostics.

!!! tip "Tradeoffs between the schemes"

    Every scheme starts from the one-loop self-energy in the Dyson–Beliaev
    equations [1, 2]. They differ in details of how the self-energy is resummed.

    `:one_loop` works in the full particle–hole (Nambu) space and keeps every
    term at order ``1/s``. Its guarantees, namely gapless Goldstone modes at
    all frequencies, conservation of total spin, and the frequency sum rule,
    rest on cancellations between the normal and anomalous self-energies, and
    between the one- and two-magnon channels. The individual channels need not
    be positive; their sum is, wherever the resummed propagator is stable.

    `:two_magnon_ladder` goes beyond one-loop order, in the spirit of the
    truncated Hilbert space exact diagonalization of Ref. [7], but in the full
    Nambu space. It is exact within the space of one- and two-magnon states; in
    a noncollinear state, the leading omitted coupling is the cubic vertex
    between a pair and three-magnon states. A static counterterm keeps the
    magnon dispersion at ``ω = 0``, and so every Goldstone mode, equal to that
    of `:one_loop`. Near a Goldstone mode, however, the bare vertices lack the
    Adler zero required by the broken symmetry, which is why the ladder can bind
    soft pairs spuriously.

    `:one_loop_on_shell` keeps only the first-order shift and width of each
    pole, ``ε̃ - iΓ = ε + Σ(ε)`` [5, Eq. (35)], evaluating the width at the
    shifted energy. This gives a clean dispersion and decay rate for each band,
    but no off-shell structure: intensities remain those of linear spin wave
    theory, and each line is Lorentzian, even where decay makes the true line
    shape strongly asymmetric. Peak energies are singular where a magnon
    crosses a van Hove singularity of the two-magnon continuum.

    `:one_loop_projected` drops the anomalous self-energy, which shifts an
    isolated magnon pole only at second order, and adds the bare two-magnon
    continuum without interference [3–6]. The spectrum is then positive by
    construction, and magnon energies are correct at order ``1/s`` away from
    soft modes. The cost is the cancellations above: near a Goldstone mode, it
    overestimates the intensity above the magnon branch and in the continuum,
    and it violates the conservation of total spin. As in Ref. [6], the
    non-resonant "source" self-energy is evaluated on shell, which removes a
    spurious branch near the ordering wavevector.

## References

1. [F. J. Dyson, _General theory of spin-wave interactions_, Phys. Rev. **102**,
   1217 (1956)](https://doi.org/10.1103/PhysRev.102.1217).
2. S. T. Beliaev, _Application of the methods of quantum field theory to a
   system of bosons_, Sov. Phys. JETP **7**,
   [289](https://jetp.ras.ru/cgi-bin/e/index/e/7/2/p289?a=list) (1958); _Energy
   spectrum of a non-ideal Bose gas_, Sov. Phys. JETP **7**,
   [299](https://jetp.ras.ru/cgi-bin/e/index/e/7/2/p299?a=list) (1958).
3. [A. B. Harris, D. Kumar, B. I. Halperin, and P. C. Hohenberg, _Dynamics of an
   antiferromagnet at low temperatures: Spin-wave damping and hydrodynamics_,
   Phys. Rev. B **3**, 961 (1971)](https://doi.org/10.1103/PhysRevB.3.961).
4. [A. L. Chernyshev and M. E. Zhitomirsky, _Spin waves in a triangular lattice
   antiferromagnet: Decays, spectrum renormalization, and singularities_, Phys.
   Rev. B **79**, 144416 (2009)](https://doi.org/10.1103/PhysRevB.79.144416).
5. [M. E. Zhitomirsky and A. L. Chernyshev, _Colloquium: Spontaneous magnon
   decays_, Rev. Mod. Phys. **85**, 219
   (2013)](https://doi.org/10.1103/RevModPhys.85.219).
6. [M. Mourigal, W. T. Fuhrman, A. L. Chernyshev, and M. E. Zhitomirsky,
   _Dynamical structure factor of triangular-lattice antiferromagnet_, Phys.
   Rev. B **88**, 094407 (2013)](https://doi.org/10.1103/PhysRevB.88.094407).
7. [Zhang et al. (2025) [arXiv:2508.21142]](https://arxiv.org/abs/2508.21142).
"""
function corrected_intensities(swt::SpinWaveTheory, qpts, scheme; energies, η, kernel=nothing, kT=0,
                               bz=nothing, mark_breakdown=true, threaded=false, verbose=false)
    num_observables(swt.measure) == 0 && error("No observables! Construct SpinWaveTheory with a `measure` argument.")
    iszero(kT) || error("Corrected intensities at kT > 0 are not yet supported.")
    (; cryst, qpts, energies, transverse, cross, direct, breakdown) =
        corrected_channels(swt, qpts, scheme; energies, η, bz, threaded, verbose)
    # The structure factor at T = 0 from the spectral function χ'', which is odd
    θ = [ω > 0 ? 1.0 : iszero(ω) ? 0.5 : 0.0 for ω in energies]
    data = reshape(θ .* (transverse + cross + direct), length(energies), size(qpts.qs)...)
    res = Intensities(cryst, qpts, energies, data)
    isnothing(kernel) || (res = broaden(res; kernel))
    if mark_breakdown
        reshape(res.data, length(energies), :)[abs.(energies) .≤ η, breakdown] .*= NaN
    end
    return res
end

# Workhorse of `corrected_intensities`, returning the three terms of S
# separately as (energy × wavevector) matrices: `transverse` from w'Gw, `cross`
# from the terms linear in K_md, and `direct` from the rest. Also returns
# `disp`, the harmonic energies; `breakdown`, a mask over wavevectors at which a
# `:one_loop` pole has moved onto the imaginary axis; for `:one_loop_on_shell`,
# `bands`, the energies, half widths and intensities of the poles; and if
# `spectral=true` then `specfunc`, the particle block of the magnon spectral
# matrix (G' - G)/2πi.
#
# The bosons are expanded about the Gaussian state of `swt`, harmonic by
# default. Its quasi-particles are the internal lines of the loops and the basis
# of the Dyson equation, its mean fields are contractions in it, and `disp`
# reports its energies.
function corrected_channels(swt::SpinWaveTheory, qpts, scheme; energies, η,
                            bz=nothing, threaded=false, verbose=false, spectral=false)
    scheme in (:one_loop, :two_magnon_ladder, :one_loop_projected, :one_loop_on_shell) ||
        error("Unknown scheme :$scheme; use :one_loop, :two_magnon_ladder, :one_loop_projected or :one_loop_on_shell.")

    (; sys, measure) = swt
    cryst = orig_crystal(sys)
    L = nbands(swt)
    Nobs = num_observables(measure)
    Ncells = nsites(sys) / natoms(cryst)

    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")
    qpts = convert(AbstractQPoints, qpts)

    # Frequencies of the response. As a range, they let the Cauchy transforms of
    # the bath be evaluated by FFT.
    zs = energies .+ im*η
    if length(energies) > 1
        dω = (energies[end] - energies[begin]) / (length(energies) - 1)
        dω > 0 && all(≈(dω), diff(energies)) || error("`energies` must be equally spaced.")
        zs = range(energies[begin] + im*η; step=dω, length=length(energies))
        if dω > (η/2) * (1 + 1e-8)
            spacing = number_to_simple_string(dω; digits=2)
            @warn """Requested `energies` are spaced by $spacing, which will not \
                     resolve features of width η = $η. A spacing of η/2 or less \
                     adds little cost."""
        end
    end
    ol = OneLoop(swt; energies, η, bz)

    chans = (; transverse = zeros(eltype(measure), length(energies), length(qpts.qs)),
               cross = zeros(eltype(measure), length(energies), length(qpts.qs)),
               direct = zeros(eltype(measure), length(energies), length(qpts.qs)))
    specfunc = spectral ? zeros(ComplexF64, L, L, length(energies), length(qpts.qs)) : nothing
    disp = zeros(L, length(qpts.qs))
    # On-shell decay rate of each band, for the `verbose` report
    linewidths = fill(NaN, L, length(qpts.qs))
    # Poles of `:one_loop_on_shell`: energies, half widths and intensities
    bands = scheme != :one_loop_on_shell ? nothing :
        (; disp = zeros(L, length(qpts.qs)), widths = zeros(L, length(qpts.qs)),
           data = zeros(eltype(measure), L, length(qpts.qs)))
    # Not a `BitVector`, whose threaded writes would race within each 64-bit chunk
    breakdown = fill(false, length(qpts.qs))

    # Nambu indices of the magnon legs, and the rows of K for the direct amplitudes
    p = 1:2L
    d = 2L .+ (1:Nobs)
    Ĩ = Diagonal([ones(L); -ones(L)])

    function calc_iq!(iq)
        corr = zeros(ComplexF64, num_correlations(measure))
        q = qpts.qs[iq]
        q_global = cryst.recipvecs * q
        se = SelfEnergy(ol, q; ladder = scheme == :two_magnon_ladder)
        (; ε, w) = se
        view(disp, :, iq) .= view(ε, 1:L)
        E = Diagonal(abs.(ε))
        (; Σ, K, poles, stable) = dyson_model(se, scheme, zs)

        # Contracts an Nobs×Nobs χ through the measure, taking S = (χ' - χ)/2πi
        function accum_channel!(accum, iω, χ)
            map!(((μ, ν),) -> (conj(χ[ν, μ]) - χ[μ, ν]) / (2π*im*Ncells), corr, measure.corr_pairs)
            accum[iω, iq] = measure.combiner(q_global, corr)
        end

        for iω in eachindex(energies)
            Kz = view(K, :, :, iω)
            (Kmd, Kdm, Kdd) = (Kz[p, d], Kz[d, p], Kz[d, d])
            G = inv(zs[iω]*Ĩ - E - Σ - Kz[p, p])
            accum_channel!(chans.transverse, iω, w' * G * w)
            # The routes through the bath complete the resolvent of the full
            # model. The other schemes assemble only the magnon term and the
            # bare continuum.
            if scheme in (:one_loop, :two_magnon_ladder)
                accum_channel!(chans.cross, iω, Kdm * G * w + w' * G * Kmd)
                accum_channel!(chans.direct, iω, Kdd + Kdm * G * Kmd)
            else
                accum_channel!(chans.direct, iω, Kdd)
            end
            isnothing(specfunc) || (view(specfunc, :, :, iω, iq) .= ((G' - G) / (2π*im))[1:L, 1:L])
        end

        # The resummed model is stable when its static magnon matrix M is
        # positive semidefinite, and breaks down for some mode at this 𝐪 when
        # it is not, not necessarily for the observed one, whose weight in S
        # may be small. An indefinite M need not have complex frequencies, so
        # the inertia is checked rather than the frequencies. That is done in
        # the boson basis, where M is congruent to T⁻†MT⁻¹: there the
        # Goldstone zero mode stays at roundoff, whereas the 1/√ε columns of
        # T near it amplify noise to a size set by `regularization`. The other
        # schemes have no anomalous self-energy in their propagator, so their
        # frequencies are always real.
        if scheme in (:one_loop, :two_magnon_ladder)
            M = E + Σ + dyson_model(se, scheme, [im*η]).K[p, p, 1]
            H = Ĩ * se.T * Ĩ * ((M + M') / 2) * Ĩ * se.T' * Ĩ
            λ = eigvals(Hermitian((H + H') / 2))
            breakdown[iq] = !stable || first(λ) < -1e-6 * maximum(abs, λ)
        end

        if scheme == :one_loop_on_shell
            for n in 1:L
                a = poles.U[:, n]' * w[1:L, :]
                map!(((μ, ν),) -> conj(a[μ]) * a[ν] / Ncells, corr, measure.corr_pairs)
                (bands.disp[n, iq], bands.widths[n, iq]) = (poles.λ[n], poles.Γ[n])
                bands.data[n, iq] = measure.combiner(q_global, corr)
            end
        end

        if verbose
            # The resummed on-shell linewidth: for `:two_magnon_ladder`, through
            # the interacting bath K, since the bare one-loop decay above is not
            # what broadens its bound-state branches.
            if scheme == :two_magnon_ladder
                K2 = dyson_model(se, :two_magnon_ladder, ε[1:L] .+ im*η).K
                view(linewidths, :, iq) .= [-imag(K2[n, n, n]) for n in 1:L]
            else
                Kdec = cauchy_transform((se.decay,), ε[1:L] .+ im*η)
                view(linewidths, :, iq) .= [-imag(Kdec[n, n, n]) for n in 1:L]
            end
        end
    end

    if verbose
        println("""
            corrected_intensities
              loop grid       $(join(ol.loop_grid, "×")) = $(prod(ol.loop_grid)) points""")
    end

    desc = verbose ? "  wavevectors     " : nothing
    elapsed = @elapsed foreach_chunked((_, iq) -> calc_iq!(iq), Returns(nothing), eachindex(qpts.qs);
                                       threaded, warn_blas=true, desc)

    if verbose
        r2 = x -> number_to_simple_string(x; digits=2)
        nthreads = threaded ? min(Threads.nthreads(), length(qpts.qs)) : 1
        per_q = 1000 * elapsed / length(qpts.qs)
        Γs = sort!(filter(isfinite, vec(linewidths)))
        report = isempty(Γs) ? "none" :
            "median $(r2(Γs[cld(end, 2)])), 90th pct $(r2(Γs[ceil(Int, 0.9end)])), against η = $(r2(η))"
        println("  elapsed         $(r2(elapsed)) s on $nthreads \
                 thread$(nthreads == 1 ? "" : "s"), $(r2(per_q)) ms per 𝐪")
        label = scheme == :two_magnon_ladder ? "on-shell Γ" : "on-shell -Im Σ"
        println("  $(rpad(label, 16))$report")
        scheme in (:one_loop, :two_magnon_ladder) && println("  breakdown       $(count(breakdown)) of $(length(qpts.qs)) wavevectors (NaN near ω = 0)")
    end

    return (; cryst, qpts, energies, chans..., specfunc, disp, breakdown, bands)
end

# The quadratic model of magnons and bath that the `scheme` propagates, given
# the one-loop self-energy `se`: a static Nambu matrix `Σ` and the bath
# transform `K` at each of the frequencies `zs`. `:one_loop_on_shell` also
# returns its `poles`. See the discussion at the head of this file.
function dyson_model(se::SelfEnergy, scheme, zs)
    (; ε, η, Σstat, decay, source) = se
    L = length(ε) ÷ 2
    p = 1:2L
    if scheme == :two_magnon_ladder
        # The bath of `:one_loop` made interacting, see `ladder_transform`. Its
        # change to the magnon block at ω = 0 is cancelled by a static
        # counterterm, which pins the stability matrix of the magnons, and with
        # it every Goldstone mode, to that of `:one_loop`, while keeping the
        # model Hermitian.
        (; K, δK0, stable) = ladder_transform(se, zs)
        return (; Σ=Σstat - δK0, K, poles=nothing, stable)
    end

    Kdec = cauchy_transform((decay,), zs)
    Ksrc = cauchy_transform((source,), zs)
    scheme == :one_loop && return (; Σ=Σstat, K=Kdec+Ksrc, poles=nothing, stable=true)

    if scheme == :one_loop_projected
        # The rotating-wave truncation, see the head of this file. The masks
        # keep decay on the particle legs and source on the hole legs. Each
        # non-resonant channel enters at real frequency in the van Vleck form
        # ½[K(ε_m) + K(ε_m′)], as in `on_shell_form`.
        Kdec[L+1:2L, :, :] .= 0
        Kdec[:, L+1:2L, :] .= 0
        Ksrc[1:L, :, :] .= 0
        Ksrc[:, 1:L, :] .= 0
        frozen(ρ) = [sum(((j, v),) -> bin_entry(v, m, m′) * (1/(ε[m] - j * ρ.Δ) + 1/(ε[m′] - j * ρ.Δ)) / 2, ρ.bins; init=0im)
                     for m in p, m′ in p]
        Pp = [m ≤ L && m′ ≤ L for m in p, m′ in p]
        Ph = [m > L && m′ > L for m in p, m′ in p]
        Σ = Σstat .* (Pp .| Ph) + frozen(source) .* Pp + frozen(decay) .* Ph
        return (; Σ, K=Kdec+Ksrc, poles=nothing, stable=true)
    end

    @assert scheme == :one_loop_on_shell
    # The particle and hole blocks of the self-energy frozen on shell, see
    # `on_shell_form`; the anomalous blocks are dropped. Each block is
    # diagonalized in its Hermitian part, which mixes bands only at first order
    # where they are nearly degenerate. The hole block at 𝐪 is the particle
    # block at -𝐪 by the Nambu symmetry of K, so the hole poles are the exact
    # mirror of the particle poles.
    #
    # Each pole takes its width from the channel that resonates at its own
    # shifted energy: decay at λ > 0 for a particle, source at -λ for a hole.
    # That channel's measure is positive semidefinite, so Γ ≥ 0 by
    # construction. Reading the width at the shifted energy rather than at the
    # harmonic one puts the threshold of the continuum where the pole sits.
    K = Kdec + Ksrc
    M = on_shell_form(se)
    (λp, Up) = eigen(Hermitian((M[1:L, 1:L] + M[1:L, 1:L]') / 2))
    (λh, Uh) = eigen(Hermitian((M[L+1:2L, L+1:2L] + M[L+1:2L, L+1:2L]') / 2))
    Kp = cauchy_transform((decay,), λp .+ im*η)
    Kh = cauchy_transform((source,), -λh .+ im*η)
    Γp = [-imag(dot(Up[:, n], view(Kp, 1:L, 1:L, n), Up[:, n])) for n in 1:L]
    Γh = [imag(dot(Uh[:, n], view(Kh, L+1:2L, L+1:2L, n), Uh[:, n])) for n in 1:L]
    # A frequency-independent Nambu matrix with exactly these poles, so that
    # the magnons propagate as a sum of unit-weight Lorentzians
    Σ = cat(Up * Diagonal(λp - im*Γp) * Up', Uh * Diagonal(λh + im*Γh) * Uh'; dims=(1, 2)) - Diagonal(abs.(ε))
    K[p, p, :] .= 0
    return (; Σ, K, poles=(; λ=λp, Γ=Γp, U=Up), stable=true)
end
