# Reproduce Fig. 2(a) and Fig. 4(a) of Mourigal, Fuhrman, Chernyshev and Zhitomirsky,
# PRB 88, 094407 (2013) = arXiv:1306.1231v3: intensity maps of the magnon spectral
# function A₁₁(𝐪, ω) and of the total structure factor S^tot(𝐪, ω) for the
# triangular-lattice Heisenberg antiferromagnet, along the path K-Γ-M-Y₁, with the
# harmonic ε_𝐪 of their Eq. (11) overlaid.
#
#   julia --project=/tmp/fig2 -t auto -e 'include("/tmp/fig2/fig2.jl"); fig2()'
#
#   fig2()                 # the published figure: fwhm = 0.03, ~30 s
#   fig2(fwhm=0.12)        # quick look, a few seconds
#   fig2(s=3/2)            # their panel (b)
#   fig2(force=true)       # ignore the cache and recompute
#
# THE ONE KNOB IS `fwhm`. The paper integrates the self-energy with an artificial
# broadening δ = 0.03JS, which at s = 1/2 is 0.015, and Γ = fwhm/2 plays exactly that
# role here, so fwhm = 2δ = 0.03 reproduces the width of their branch. Everything else
# follows from it and is chosen automatically: the frequency grid needs dω ≲ Γ/2 to
# sample a peak of half-width Γ, and the loop grid needs its pair-energy spacing
# |∇ε|/nk ≲ Γ or the continuum comes out speckled. Overriding `nw` or `nk` below those
# values is what produces artifacts; the defaults warn rather than guess again.
#
# Their A₁₁ is the diagonal Nambu component of the *magnon* propagator: Eq. (14) applies
# the (u_𝐪 ± v_𝐪)² and Λ± factors outside it, so what is plotted is the particle block
# of the Bogoliubov-basis Green function, which is exactly the object
# `intensities_corrected` contracts with observable vectors. Their model has one branch
# per wavevector over the whole zone, while our three-site cell folds 𝐪 together with
# 𝐪 ± 𝐐, 𝐐 = [1/3, 1/3, 0]; the branch belonging to 𝐪 itself is the one whose harmonic
# energy is ε_𝐪, so it is picked out by matching.
#
# The Dyson assembly is duplicated from CorrectedIntensities.jl rather than called,
# because A₁₁ is not exposed by the public routine and evaluating it here lets the total
# structure factor be built from the same self-energy instead of paying for it twice.
# The duplication must be kept in step: the equation is solved in the L×L particle
# block, and the source channel of the self-energy is frozen at its on-shell frequency.
# Both are what keep A₁₁ non-negative and bounded by 1/πΓ; see the comments in
# CorrectedIntensities.jl and SelfEnergy.jl. `compute` cross-checks itself against
# `intensities_corrected` at two wavevectors, and with a shared loop grid the two agree
# to round-off, so any drift in the duplication shows up immediately.

using Sunny, LinearAlgebra, Printf, Statistics, Serialization
using CairoMakie
CairoMakie.activate!(type="png")

# Threading this calculation is pointless with the OpenBLAS that Julia ships. Its
# small-matrix calls do not merely fail to scale, they anti-scale: called from 18
# threads, `bogoliubov!` slows by a factor of 30 and a 36×6×6 `gemm` by a factor of 22,
# so the per-𝐪 cost grows in proportion to the number of threads and total throughput
# saturates at two cores' worth. It is contention inside the library, not in our code,
# and `BLAS.set_num_threads(1)` does not touch it. Apple's Accelerate is reentrant:
# scaling goes 0.6× → 12.6× for `bogoliubov!` and 2.3× → 11.6× for `vertex!`, at the
# price of a slightly slower serial eigensolve. On a machine with 18 threads that is a
# factor of four on the whole transverse calculation.
Sys.isapple() && @eval using AppleAccelerate

const DIR = @__DIR__
const OPTS = (; maxevals = 20000)   # momentum integrals for the mean fields

say(args...) = (println(args...); flush(stdout))

# Harmonic dispersion of the unfolded model in closed form, Eq. (11), in the RLU of
# the original one-site cell
γq(q) = (cos(2π*q[1]) + cos(2π*q[2]) + cos(2π*(q[1] + q[2]))) / 3
ε11(q, s) = 3s * sqrt(max(0, (1 - γq(q)) * (1 + 2γq(q))))

# 120° spiral on the triangular lattice, in the three-site magnetic cell
function htaf(s)
    cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(; s, g=2)], :dipole)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
    Q = cryst.recipvecs * [1/3, 1/3, 0]
    for site in eachsite(sys)
        θ = dot(Q, global_positions(sys)[site])
        set_dipole!(sys, [cos(θ), sin(θ), 0], site)
    end
    return (; swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false)), cryst)
end

# Everything the figure needs, for one set of parameters. Results are cached on disk
# keyed by those parameters, so re-rendering never recomputes.
function compute(; s=1/2, fwhm=0.03, ωmax=3.0, npath=241,
                 nw=1 + round(Int, 4ωmax/fwhm), nk=8*ceil(Int, 3/(4fwhm)), force=false)
    file = @sprintf("%s/cache/s%.2f_fw%.3f_np%d_nw%d_nk%d.jls", DIR, s, fwhm, npath, nw, nk)
    !force && isfile(file) && return deserialize(file)

    Γ = fwhm / 2
    dω = ωmax / (nw - 1)
    dω > Γ/2 && @warn "Frequency grid too coarse to sample the peaks" dω Γ
    nk < 6/fwhm && @warn "Loop grid too coarse for this broadening; expect speckle" nk

    (; swt, cryst) = htaf(s)
    sys = swt.sys
    (; measure) = swt
    L = Sunny.nbands(swt)
    Ncells = Sunny.nsites(sys) / Sunny.natoms(cryst)

    # K = ordering wavevector = zone corner; M = edge midpoint; Y₁ = M + 𝐐, the point
    # the paper singles out as the "blow-out" region of strong decay
    Kpt = [1/3, 1/3, 0]
    Mpt = [1/2, 0, 0]
    path = q_space_path(cryst, [Kpt, [0, 0, 0], Mpt, Mpt + Kpt], npath; labels=["K", "Γ", "M", "Y₁"])
    qs = collect(path.qs)
    energies = collect(range(0, ωmax, nw))
    kernel = lorentzian(; fwhm)

    say(@sprintf("computing s = %.2f, fwhm = %.3f: %d 𝐪 × %d ω, loop grid %d², %d threads",
                 s, fwhm, npath, nw, nk, Threads.nthreads()))
    t0 = time()

    # Mean-field corrections, all 𝐪-independent
    tad = Sunny.tadpole_correction(swt; OPTS...)
    terms2 = [Sunny.hartree_fock_correction(swt; OPTS...).terms2
              tad.terms2
              Sunny.anisotropy_correction(swt).terms2]
    δc = Sunny.observable_corrections(swt; v=tad.v, OPTS...)
    terms3 = Sunny.cubic_monomials(swt)
    ps = Sunny.loop_grid((nk, nk, 1))

    Asel = zeros(nw, npath)      # spectral function of the branch belonging to 𝐪
    Aall = zeros(nw, npath)      # trace over the folded bands, for the artifact detector
    Strans = zeros(nw, npath)    # transverse structure factor
    Slong = zeros(nw, npath)     # two-magnon continuum
    εref = [ε11(q, s) for q in qs]
    amin = fill(Inf, npath)      # least eigenvalue of A, which must be ≥ 0

    Threads.@threads for iq in eachindex(qs)
        q = qs[iq]
        q_reshaped = Sunny.to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q

        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        δH = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Sunny.num_observables(measure))
        Σ3 = zeros(ComplexF64, L, L, nw)
        corr = zeros(ComplexF64, Sunny.num_correlations(measure))

        ε = Sunny.excitations!(T, H, swt, q)
        Sunny.accum_quadratic!(δH, terms2, q_reshaped)
        Σstat = Diagonal([ones(L); -ones(L)]) * transpose(T' * δH * T)
        onshell = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
        Sunny.accum_cubic_self_energy!(Σ3, swt, terms3, q_reshaped, energies .+ im*Γ, ps, 0.0;
                                       source_freqs=onshell, bin_width=Γ/16)

        Sunny.set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
        Sunny.accum_observable_corrections!(u, swt, q_reshaped, q_global, δc)
        w = T' * u

        # The folded band belonging to 𝐪 itself, identified by its harmonic energy
        n = argmin(abs.(view(ε, 1:L) .- εref[iq]))

        for (iω, ω) in enumerate(energies)
            G = inv((ω + im*Γ)*I - Diagonal(view(ε, 1:L)) - view(Σstat, 1:L, 1:L) - view(Σ3, :, :, iω))
            A = (G - G') / 2im
            Asel[iω, iq] = -real(A[n, n]) / π
            Aall[iω, iq] = -real(tr(A)) / π
            amin[iq] = min(amin[iq], minimum(eigvals(Hermitian(-A/π))))
            map!(corr, measure.corr_pairs) do (μ, ν)
                -dot(view(w, 1:L, μ), A, view(w, 1:L, ν)) / (π * Ncells)
            end
            Strans[iω, iq] = real(measure.combiner(q_global, corr))
        end

        Slong[:, iq] = Sunny.intensities_two_magnon(swt, [q]; energies, kernel, grid=(nk, nk, 1)).data
    end
    t_spec = time() - t0

    # Guard on the duplicated assembly: same loop grid, so agreement is to round-off
    dev = maximum([1 + npath ÷ 3, 1 + 2npath ÷ 3]) do iq
        r = Sunny.intensities_corrected(swt, [qs[iq]]; energies, kernel, grid=(nk, nk, 1), OPTS...)
        mine = Strans[:, iq] + Slong[:, iq]
        maximum(abs, vec(r.data) - mine) / maximum(abs, mine)
    end

    say(@sprintf("  spectra %.0f s, cross-check %.0f s; \
                  least eigenvalue of A %+.2e, peak A₁₁ %.2f (bound 1/πΓ = %.2f), \
                  assembly matches to %.1e",
                 t_spec, time() - t0 - t_spec,
                 minimum(amin), maximum(Asel), 1/(π*Γ), dev))

    data = (; s, fwhm, nk, npath, nw, energies, Asel, Aall, Strans, Slong, εref,
              xticks = path.xticks, amin, dev)
    mkpath("$DIR/cache")
    serialize(file, data)
    return data
end

# Peak position of the magnon branch. It is the peak descended from the harmonic pole,
# so it is sought in a window bracketing ε_𝐪 rather than as a global maximum: elsewhere
# in the range the continuum, and the divergence at the Goldstone wavevectors, are both
# larger. A renormalization outside -10%..+55% would fall outside this window, so the
# number reported is only meaningful because it lands well inside it.
function peak_energies(d, window=(0.45, 1.10))
    map(eachindex(d.εref)) do iq
        m = findall(e -> window[1]*d.εref[iq] <= e <= window[2]*d.εref[iq] + 1e-9, d.energies)
        isempty(m) && return NaN
        return d.energies[m[argmax(view(d.Asel, m, iq))]]
    end
end

# Detector for the artifact that motivated the present assembly, kept as a regression
# check: it should find nothing. Inverting the full 2L Nambu denominator pushed a mirror
# pole up through ω = 0 at isolated wavevectors, and because the resulting root lay
# closer to the real axis than Γ the pole sat on the wrong side of it and the spectral
# density went negative. Two symptoms: weight piled up below the branch, where a magnon
# cannot be, and weight missing from the column altogether. The second is measured
# against the median of a seven-point neighbourhood, which a dispersing feature survives
# and an isolated column does not. Both tests are restricted to ε_𝐪 > 0.5, because near
# K and Γ the weight really does vary sharply: the occupation ⟨n̂_𝐪⟩ diverges there.
function suspect_columns(d)
    nq = length(d.εref)
    wcol = [sum(view(d.Aall, :, iq)) for iq in 1:nq]
    return filter(1:nq) do iq
        d.εref[iq] < 0.5 && return false
        lo = findall(<(0.35 * d.εref[iq]), d.energies)
        piled = !isempty(lo) && maximum(view(d.Aall, lo, iq)) > 1.0
        nb = median(wcol[[j for j in iq-3:iq+3 if 1 <= j <= nq && j != iq]])
        return piled || abs(wcol[iq] - nb) > 0.2 * nb
    end
end

# Styling follows the published figures so the two can be compared by eye:
#
#   * `:jet`, the MATLAB/Mathematica rainbow they used. Its dark navy at zero is what
#     makes the background of their panels read as empty, and it puts a lot of contrast
#     in the lowest tenth of the range, where the continuum lives.
#   * a *linear* colour scale, cut off well below the peak of the map: (0, 3) for A₁₁ as
#     in their Fig. 2, and (0, 2) for S^tot as in their Fig. 4, whose colourbar is
#     labelled ">2" at the top for that reason. Both maps have divergences that are real
#     and are not the subject of the figure — A₁₁ at Γ, where the Bogoliubov factors
#     entering the static self-energy diverge, and S^tot at K, the magnetic Bragg point
#     — and both are allowed to saturate. Only a percent or two of pixels clip.
#
# An earlier version used a square-root transform over the full data range. That showed
# the continuum better but compressed the branch into a single hue, and looked nothing
# like the paper. The linear scale only reads like theirs at their broadening, since at
# larger fwhm the branch spreads over more pixels at lower amplitude and little of it
# saturates.
function render(d; file="$DIR/fig2.png")
    nq = length(d.εref)
    Stot = d.Strans + d.Slong
    bad = suspect_columns(d)

    fig = Figure(size=(1150, 470), fontsize=15)
    panels = [(d.Asel, "A₁₁(𝐪, ω)   — cf. Mourigal Fig. 2(a)", 0:1:3, "A₁₁"),
              (Stot, "Sᵗᵒᵗ(𝐪, ω) = transverse + continuum   — cf. Fig. 4(a)", 0:0.5:2, "Sᵗᵒᵗ")]

    for (col, (data, title, cticks, clab)) in enumerate(panels)
        ax = Axis(fig[1, 2col-1]; title, ylabel = col == 1 ? "ω / J" : "",
                  xticks=d.xticks, xgridvisible=true, xgridcolor=(:white, 0.5),
                  ygridvisible=false)
        # The particle-block assembly makes the density non-negative by construction,
        # so the clip at zero is inert; the minimum is reported below to confirm it.
        hm = heatmap!(ax, 1:nq, d.energies, permutedims(max.(data, 0));
                      colormap=:jet, colorrange=(0, last(cticks)))
        lines!(ax, 1:nq, d.εref; color=(:white, 0.85), linestyle=:dash, linewidth=1.5)
        xlims!(ax, 1, nq)
        ylims!(ax, 0, maximum(d.energies))
        labs = [t == last(cticks) ? "≥$t" : string(t) for t in cticks]
        Colorbar(fig[1, 2col], hm; label=clab, width=12, ticks=(collect(cticks), labs))
        colgap!(fig.layout, 2col-1, 8)
    end

    Label(fig[0, 1:4],
          @sprintf("Triangular-lattice Heisenberg antiferromagnet, s = %.1f, 1/s-corrected \
                    (Lorentzian fwhm = %.2f J, self-energy grid %d×%d)", d.s, d.fwhm, d.nk, d.nk),
          fontsize=16, font=:bold)
    rowgap!(fig.layout, 6)
    save(file, fig; px_per_unit=2)
    say("wrote $file")

    # Renormalization of the branch, to compare with their Sec. III.1 figure of ~18% for
    # S = 1/2. Only wavevectors where the peak is clear of the elastic and Goldstone
    # regions are counted, and the spread over search windows says how well defined it is.
    ωpk = peak_energies(d)
    ok = [iq for iq in 1:nq if d.εref[iq] > 0.3 && !(iq in bad)]
    r = 1 .- ωpk[ok] ./ d.εref[ok]
    windows = ((0.35, 1.15), (0.45, 1.10), (0.55, 1.05), (0.60, 1.20), (0.70, 1.10))
    spread = map(w -> mean(iq -> 1 - peak_energies(d, w)[iq]/d.εref[iq], ok), windows)
    say(@sprintf("\nrenormalization over %d wavevectors (ε_𝐪 > 0.3, %d suspect columns dropped):\
                  \n   mean %.1f%%, range %.1f%%..%.1f%%, and %.1f%%..%.1f%% across search windows\
                  \n   paper, Sec. III.1: \"about 18%% for S=1/2\"",
                 length(ok), length(bad), 100mean(r), 100minimum(r), 100maximum(r),
                 100minimum(spread), 100maximum(spread)))
    say(@sprintf("suspect columns: %d of %d %s", length(bad), nq, string(bad)))
    say(@sprintf("least density: A₁₁ %.2e (peak %.2f, bound %.2f), Sᵗᵒᵗ %.2e; clipped pixels %.2f%% / %.2f%%",
                 minimum(d.Asel), maximum(d.Asel), 2/(π*d.fwhm), minimum(Stot),
                 100mean(>(3), d.Asel), 100mean(>(2), Stot)))

    say(@sprintf("\n%-4s %8s %8s %8s %10s %10s", "pt", "ε_har", "ω_peak", "shift", "contin.", "A₁₁ peak"))
    for (lab, i) in zip(d.xticks[2], d.xticks[1])
        i = clamp(round(Int, i), 1, nq)
        cont = sum(view(d.Slong, :, i)) / sum(view(Stot, :, i))
        shift = d.εref[i] > 0.3 ? @sprintf("%6.1f%%", 100*(1 - ωpk[i]/d.εref[i])) : "     —"
        say(@sprintf("%-4s %8.3f %8.3f %8s %9.1f%% %10.2f", lab, d.εref[i], ωpk[i], shift,
                     100*cont, maximum(view(d.Asel, :, i))))
    end
    return fig
end

fig2(; file="$DIR/fig2.png", kw...) = render(compute(; kw...); file)
