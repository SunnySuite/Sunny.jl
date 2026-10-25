# Dynamical structure factor at O(1/s), S = (χ' - χ)/2πi at z = ω + iη.
#
# Each `dyson` scheme propagates the magnons with the 2L×2L Nambu resolvent G
# of a quadratic model of magnons coupled to the binned two-magnon bath of
# Corrections.jl. The schemes make three independent choices: which blocks of
# the self-energy to keep, whether to invert the resolvent or take its
# first-order poles, and how much of the model's response to assemble.
#
# `:nambu` keeps the whole model and assembles the Schur complement of its
# resolvent onto the observables,
#
#     χ = w'Gw + K_dm G w + w'G K_md + K_dm G K_md + K_dd,
#
# which inherits the frequency sum rule and is positive at ω > 0, up to the
# O(η) tails of the mirror poles, whenever the model is stable. `:ladder` is
# the same model with the pair interaction of H₄ added to the bath (see
# `ladder_transform`), which is solved exactly and so keeps all of this. The
# other two schemes assemble only the magnon term w'Gw and the bare continuum
# K_dd.
#
# `:particle` is the rotating-wave truncation. Each bath couples only to the
# legs it resonates with, decay pairs to particles and source pairs to holes,
# and the anomalous blocks are dropped. Each block keeps its non-resonant
# channel frozen at the mean on-shell energy of its two legs, as in PRB 88,
# 094407 (2013), which removes a spurious branch pushed up from negative
# frequency near ±Q. The particle block is then G₁₁ of RMP 85, 219, Eq. (34)
# and PRB 79, 144416, Eq. (86), and the hole block its mirror G₁₁(-𝐪, -ω). The
# dropped counter-rotating terms are suppressed by |Σ_anom|/(ω + ε), so they
# are O(1/s) in the intensities and not small near a soft mode whose rotation
# generator has a quadratic part (e.g. at ±K of the triangular lattice), where
# the 1/ε divergences cancel only in the full Nambu inverse. There the scheme
# overestimates the intensity at ω ≫ ε.
#
# `:on_shell` linearizes the full Nambu Dyson equation at each pole, as in RMP
# Eq. (35) and PRB Eq. (55). The anomalous blocks couple poles at ±ε, a gap of
# 2ε, and so first shift a pole at second order: the same poles follow from the
# particle block alone. Each pole is a unit-weight Lorentzian; its residue,
# which does change at first order, is deliberately not corrected, keeping the
# scheme a statement about the dispersion. The details, including which channel
# sets each width, are at `dyson_model`. The scheme is controlled where the
# anomalous coupling is small against 2ε. Near a Goldstone mode in two
# dimensions it is not: two soft internal lines feed a static bubble that
# diverges like 1/δ at distance δ, in the stiff direction directly and in the
# Goldstone direction through its O(δ) overlap with that one, so that the
# Goldstone block goes as δ(b s δ - c) instead of b s δ². The Ward identity
# protects only δ = 0. Inside δ ≲ c/(b s), the resummation of `:nambu` has
# complex frequencies, and |Σ_anom|/2ε exceeds one, so `:on_shell` there is a
# first-order formula outside its range.
#
# No scheme dresses the internal lines by default: RMP Sec. IV.C notes that
# self-consistent dressing gaps the Goldstone modes. A `MagnonVacuum` with a
# correction does so, with a counterterm.

"""
    corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, kernel=nothing, dyson=:nambu,
                          vacuum=MagnonVacuum(swt), grid=auto_bzgrid(; η, vacuum, tol=0.01),
                          mark_breakdown=true, threaded=false, verbose=false)

Dynamical spin structure factor at temperature ``T = 0``, including one-loop
spin-wave corrections, i.e., relative order ``1/s`` in dipole mode. Magnon
energies are shifted by the mean fields and the cubic self-energy, magnons that
can decay into two magnons are broadened, and the response includes the
two-magnon continuum.

The Green function is evaluated at ``ω + iη``, so `η` both regulates the loop
integrals and broadens the result into Lorentzians. Choose it small compared to
the linewidths of interest, but no smaller: the wavevector grid of the loop
integrals grows as `1/η` in each dispersing dimension. An optional `kernel`,
e.g. for instrumental resolution, can be used to postprocess the result.

Loop integrals over the Brillouin zone are approximated as a discrete sum over
the points of the [`BZGrid`](@ref). The [`auto_bzgrid`](@ref) default resolves
`η` in the dispersion of `vacuum` to a relative accuracy of `tol=0.01`. The
static mean fields are summed on the same grid, which keeps Goldstone modes
exactly gapless. A self-consistent `vacuum` must be given the grid on which it
was solved.

The `dyson` option selects how the one-loop self-energy is resummed:

- `:nambu` (default) solves the Dyson equation in the full particle/hole (Nambu)
  space and assembles the complete one-loop response, including the interference
  between one- and two-magnon channels. It preserves Goldstone modes at all
  frequencies and satisfies the frequency sum rule. The intensity is positive at
  ``ω > 0`` wherever the resummed propagator is stable.
- `:ladder` extends `:nambu` with the interaction between the two magnons of
  each pair, summed to all orders. This produces two-magnon bound states and
  resonances, as in the truncated Hilbert space exact diagonalization of [Zhang
  et al., arXiv:2508.21142](https://arxiv.org/abs/2508.21142), but in the full
  Nambu space. A static counterterm keeps the magnon dispersion at ``ω = 0``,
  and so every Goldstone mode, exactly that of `:nambu`. In a gapped magnet
  this counterterm removes a physical shift, so that the result departs
  slightly from the exact ladder. The ladder is beyond one-loop order. Near a
  Goldstone mode the off-shell quartic vertex lacks its Adler zero, and in two
  dimensions it binds pairs of soft magnons below ``ω = 0`` as the loop grid is
  refined. Where the cubic vertex vanishes, i.e. for collinear order in dipole
  mode, the pair bath is rotated to restore the Adler zero, which keeps it
  stable on any grid. Otherwise, including collinear order in `:SUN` mode, the
  bare ladder is used. In two dimensions its instability is then reported as a
  breakdown; in three, results near the Goldstone mode are sensitive to the loop
  grid. Gapped magnets have no Adler zero to protect.
- `:particle` keeps only the particle block of the self-energy, and collects the
  magnon response plus the bare two-magnon continuum without interference. The
  intensity is never negative at ``ω > 0``. Pole shifts are correct at order
  ``1/s`` away from soft modes. Near a Goldstone mode it can put spurious
  intensity well above the magnon energy.
- `:on_shell` linearizes the Dyson equation at each magnon pole, so that each
  magnon becomes one Lorentzian of shifted energy and nonnegative width, with
  its harmonic intensity. Band energies and widths are correct at order ``1/s``,
  band intensities are not. The continuum is added as in `:particle`. Peak
  energies have logarithmic singularities where a magnon crosses a saddle point
  of the two-magnon continuum.

Sunny uses the retarded response, i.e., includes the tails of the
negative-frequency poles. These cancel the high-frequency Lorentzian tail of a
Goldstone mode.

With `:nambu`, a large one-loop correction can push a magnon pole past ``ω = 0``
and onto the imaginary axis. This condition signals that the perturbative
expansion has broken down at that wavevector. Sunny marks these points by
setting all intensity within `±η` of ``ω = 0`` to `NaN`. Use
`mark_breakdown=false` to keep the raw data. Note that a different Dyson
resummation scheme cannot rescue perturbation theory at the same one-loop order;
schemes `:particle` and `:on_shell` are similarly uncontrolled, even though they
do not show imaginary poles. This breakdown is most severe in low dimensions. In
two dimensions, for example, if the ordered state has cubic magnon vertices, one
expects perturbation theory to fail within a distance of order ``1/s`` of each
Goldstone wavevector.

Set `threaded=true` to parallelize over `qpts`, and `verbose=true` to print a
progress bar and diagnostics.

!!! tip "Origins and accuracy of the schemes"

    Interacting spin waves were treated diagrammatically by [Dyson, Phys. Rev.
    **102**, 1217 (1956)](https://doi.org/10.1103/PhysRev.102.1217). The coupled
    normal and anomalous Dyson equations of an interacting Bose system are due to
    Beliaev, Sov. Phys. JETP **7**,
    [289](https://jetp.ras.ru/cgi-bin/e/index/e/7/2/p289?a=list) and
    [299](https://jetp.ras.ru/cgi-bin/e/index/e/7/2/p299?a=list) (1958). The
    `:nambu` scheme solves these Dyson-Beliaev equations at full one-loop order.

    In the harmonic quasiparticle basis, the anomalous self-energy couples a pole at
    ``+ε`` to its mirror at ``-ε``, and so shifts an isolated pole only at second
    order. This observation has historically motivated projection onto the particle
    block [Chernyshev and Zhitomirsky, PRB **79**, 144416
    (2009)](https://doi.org/10.1103/PhysRevB.79.144416) and [Zhitomirsky and
    Chernyshev, RMP **85**, 219 (2013)](https://doi.org/10.1103/RevModPhys.85.219).
    These prior works collected the renormalized magnon response and the bare
    two-magnon continuum independently. Neglecting interference between the two is
    necessary to ensure a positive spectrum, but leaves the scheme incomplete at
    one-loop order. Sunny's `:particle` scheme follows this Chernyshev and
    Zhitomirsky recipe precisely.

    The `:on_shell` poles are ``ε̃ - iΓ = ε + Σ(ε)``, as in Eq. (55) of the PRB.
    Because the anomalous self-energy shifts an isolated pole only at second order,
    these same poles follow from the full Dyson-Beliaev equation or from its
    projection to the particle block. Each width is read at the shifted energy, so
    that decay begins at the renormalized two-magnon threshold. Both the projection
    and the linearization assume that the anomalous coupling is small against the
    ``2ε`` that separates a pole from its mirror.
"""
function corrected_intensities(swt::SpinWaveTheory, qpts; energies, η, kernel=nothing, dyson=:nambu,
                               vacuum=MagnonVacuum(swt), grid::BZGrid=auto_bzgrid(; η, vacuum, tol=0.01),
                               mark_breakdown=true, threaded=false, verbose=false)
    (; cryst, qpts, energies, transverse, cross, direct, breakdown) =
        corrected_channels(swt, qpts; energies, η, dyson, vacuum, grid, threaded, verbose)
    data = reshape(transverse + cross + direct, length(energies), size(qpts.qs)...)
    res = Intensities(cryst, qpts, energies, data)
    isnothing(kernel) || (res = broaden(res; kernel))
    if mark_breakdown
        reshape(res.data, length(energies), :)[abs.(energies) .≤ η, breakdown] .= NaN
    end
    return res
end

"""
    corrected_intensities_bands(swt::SpinWaveTheory, qpts; η, vacuum=MagnonVacuum(swt),
                                grid=auto_bzgrid(; η, vacuum, tol=0.01),
                                threaded=false, verbose=false)

Magnon bands at temperature ``T = 0`` with one-loop corrections, i.e., relative
order ``1/s`` in dipole mode, for fitting a measured dispersion. These are the
poles of the `dyson=:on_shell` option of [`corrected_intensities`](@ref), which
describes the scheme and its limits; the two-magnon continuum is omitted.

Each band carries a shifted energy and a half width at half maximum in the
field `widths`, arising from decay into the two-magnon continuum, both correct
at order ``1/s``. The intensity of each band is that of linear spin wave
theory, for the observables corrected at order ``1/s``; it omits the
redistribution of weight at that order. Here `η` regularizes the loop
integrals only, and does not broaden the result. The widths are nonnegative,
and the bands of ``𝐪`` mirror the hole poles at ``-𝐪`` exactly.
"""
function corrected_intensities_bands(swt::SpinWaveTheory, qpts; η, vacuum=MagnonVacuum(swt),
                                     grid::BZGrid=auto_bzgrid(; η, vacuum, tol=0.01),
                                     threaded=false, verbose=false)
    (; cryst, qpts, bands) = corrected_channels(swt, qpts; energies=Float64[], η, dyson=:on_shell,
                                                vacuum, grid, threaded, verbose)
    sz = (size(bands.disp, 1), size(qpts.qs)...)
    return BandIntensities(cryst, qpts, reshape(bands.disp, sz), reshape(bands.data, sz), reshape(bands.widths, sz))
end

# Workhorse of `corrected_intensities`, returning the three terms of S
# separately as (energy × wavevector) matrices: `transverse` from w'Gw, `cross`
# from the terms linear in K_md, and `direct` from the rest. Also returns
# `disp`, the harmonic energies; `breakdown`, a mask over wavevectors at which a
# `:nambu` pole has moved onto the imaginary axis; for `:on_shell`, `bands`, the
# energies, half widths and intensities of the poles; and if `spectral=true`
# then `specfunc`, the particle block of the magnon spectral matrix (G' -
# G)/2πi.
#
# The bosons are expanded about `vacuum`, harmonic by default. Its
# quasi-particles are the internal lines of the loops and the basis of the Dyson
# equation, its mean fields are contractions in its vacuum, and `disp` reports
# its energies.
function corrected_channels(swt::SpinWaveTheory, qpts; energies, η, dyson=:nambu,
                            vacuum=MagnonVacuum(swt), grid::BZGrid=auto_bzgrid(; η, vacuum, tol=0.01),
                            threaded=false, verbose=false, spectral=false)
    dyson in (:nambu, :ladder, :particle, :on_shell) ||
        error("Unknown `dyson=:$dyson`; use :nambu, :ladder, :particle or :on_shell.")

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
            @warn """Requested `energies` are spaced by $(round(dω, sigdigits=2)) on \
                     average, which will not resolve features of width η = $η. A spacing \
                     of η/2 or less adds little cost."""
        end
    end
    ol = OneLoop(swt; energies, η, vacuum, grid)

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
        se = SelfEnergy(ol, q; ladder = dyson == :ladder)
        (; ε, w) = se
        view(disp, :, iq) .= view(ε, 1:L)
        E = Diagonal(abs.(ε))
        (; Σ, K, poles, stable) = dyson_model(se, dyson, zs)

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
            if dyson in (:nambu, :ladder)
                accum_channel!(chans.cross, iω, Kdm * G * w + w' * G * Kmd)
                accum_channel!(chans.direct, iω, Kdd + Kdm * G * Kmd)
            else
                accum_channel!(chans.direct, iω, Kdd)
            end
            isnothing(specfunc) || (view(specfunc, :, :, iω, iq) .= ((G' - G) / (2π*im))[1:L, 1:L])
        end

        # A pole pushed through ω = 0 collides with its mirror, and the pair
        # moves onto the imaginary axis: a mode frequency of Ĩ·M at ω = 0 has
        # an imaginary part. This is a breakdown of the one-loop resummation
        # for some mode at this 𝐪, not necessarily for the observed one; its
        # weight in S may be small. The other schemes have no anomalous
        # self-energy in their propagator, so their frequencies are always
        # real. The check is made just above ω = 0.
        if dyson in (:nambu, :ladder)
            M = E + Σ + dyson_model(se, dyson, [im*η]).K[p, p, 1]
            breakdown[iq] = !stable || any(z -> abs(imag(z)) > 1e-8 * opnorm(M), eigvals(Ĩ * (M + M') / 2))
        end

        if dyson == :on_shell
            for n in 1:L
                a = poles.U[:, n]' * w[1:L, :]
                map!(((μ, ν),) -> conj(a[μ]) * a[ν] / Ncells, corr, measure.corr_pairs)
                (bands.disp[n, iq], bands.widths[n, iq]) = (poles.λ[n], poles.Γ[n])
                bands.data[n, iq] = measure.combiner(q_global, corr)
            end
        end

        if verbose
            # The resummed on-shell linewidth: for `:ladder`, through the
            # interacting bath K, since the bare one-loop decay above is not
            # what broadens its bound-state branches.
            if dyson == :ladder
                K2 = dyson_model(se, :ladder, ε[1:L] .+ im*η).K
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

    t0 = time()
    desc = verbose ? "  wavevectors     " : nothing
    foreach_chunked((_, iq) -> calc_iq!(iq), Returns(nothing), eachindex(qpts.qs);
                    threaded, warn_blas=true, desc)
    elapsed = time() - t0

    if verbose
        r2 = x -> round(x; sigdigits=2)
        nthreads = threaded ? min(Threads.nthreads(), length(qpts.qs)) : 1
        per_q = 1000 * elapsed / length(qpts.qs)
        Γs = sort!(filter(isfinite, vec(linewidths)))
        report = isempty(Γs) ? "none" :
            "median $(r2(Γs[cld(end, 2)])), 90th pct $(r2(Γs[ceil(Int, 0.9end)])), against η = $(r2(η))"
        println("  elapsed         $(round(elapsed; digits=1)) s on $nthreads \
                 thread$(nthreads == 1 ? "" : "s"), $(r2(per_q)) ms per 𝐪")
        label = dyson == :ladder ? "on-shell Γ" : "on-shell -Im Σ"
        println("  $(rpad(label, 16))$report")
        dyson in (:nambu, :ladder) && println("  breakdown       $(count(breakdown)) of $(length(qpts.qs)) wavevectors (NaN near ω = 0)")
    end

    return (; cryst, qpts, energies, chans..., specfunc, disp, breakdown, bands)
end

# The quadratic model of magnons and bath that the `dyson` scheme propagates,
# given the one-loop self-energy `se`: a static Nambu matrix `Σ` and the bath
# transform `K` at each of the frequencies `zs`. `:on_shell` also returns its
# `poles`. See the discussion at the head of this file.
function dyson_model(se::SelfEnergy, dyson, zs)
    (; ε, η, Σstat, decay, source) = se
    L = length(ε) ÷ 2
    p = 1:2L
    if dyson == :ladder
        # The bath of `:nambu` made interacting, see `ladder_transform`. Its
        # change to the magnon block at ω = 0 is cancelled by a static
        # counterterm, which pins the stability matrix of the magnons, and with
        # it every Goldstone mode, to that of `:nambu`, while keeping the model
        # Hermitian.
        (; K, δK0, stable) = ladder_transform(se, zs)
        return (; Σ=Σstat - δK0, K, poles=nothing, stable)
    end

    Kdec = cauchy_transform((decay,), zs)
    Ksrc = cauchy_transform((source,), zs)
    dyson == :nambu && return (; Σ=Σstat, K=Kdec+Ksrc, poles=nothing, stable=true)

    if dyson == :particle
        # The rotating-wave truncation of the auxiliary model. Each bath couples
        # only to the legs it resonates with, decay pairs to particles and
        # source pairs to holes, the anomalous blocks are dropped, and each
        # block keeps its non-resonant channel frozen on shell, at the mean
        # on-shell energy of its two legs, which is Hermitian.
        Kdec[L+1:2L, :, :] .= 0
        Kdec[:, L+1:2L, :] .= 0
        Ksrc[1:L, :, :] .= 0
        Ksrc[:, 1:L, :] .= 0
        frozen(ρ) = [sum(((j, v),) -> bin_entry(v, m, m′) / ((ε[m] + ε[m′])/2 - j * ρ.Δ), ρ.bins; init=0im)
                     for m in p, m′ in p]
        Pp = [m ≤ L && m′ ≤ L for m in p, m′ in p]
        Ph = [m > L && m′ > L for m in p, m′ in p]
        Σ = Σstat .* (Pp .| Ph) + frozen(source) .* Pp + frozen(decay) .* Ph
        return (; Σ, K=Kdec+Ksrc, poles=nothing, stable=true)
    end

    @assert dyson == :on_shell
    # The Dyson equation of `:nambu`, linearized at each pole, with the
    # self-energy frozen on shell as in `on_shell_form`. The anomalous
    # blocks couple poles at ±ε, a gap of 2ε, so they first shift a pole at
    # second order and are dropped. Each of the particle and hole blocks is
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
