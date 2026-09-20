# Reproduce Fig. 3(a) of Maksimov, Zhitomirsky and Chernyshev, PRB 94, 140407(R) (2016)
# = arXiv:1607.08238: the dynamical structure factor 𝒮(𝐪, ω) of the easy-plane XXZ
# triangular-lattice antiferromagnet in an out-of-plane field, including the corrections
# of order 1/s, for s = 1/2, Δ = 0.9 and H = 0.2 Hs, along the path M-K′-Γ-K-M-Γ.
#
#   julia --project=/tmp/fig2 -t auto -e 'include("<this directory>/xxz_fig3.jl"); xxz_fig3()'
#
# The Mourigal scripts in this directory (fig2.jl, fig4b.jl) cover the zero-field
# Heisenberg triangular antiferromagnet, whose 120° spiral is coplanar. This model adds
# exchange anisotropy and a field, and its ground state is the *umbrella*: the three
# sublattices keep the 120° winding in the easy plane and cant uniformly out of it. The
# umbrella carries a staggered scalar chirality 𝐒ᵢ·(𝐒ⱼ×𝐒ₖ), which breaks inversion, so
# ε_𝐤 ≠ ε_{-𝐤} and the two classes of zone corner are no longer equivalent: the corners
# in the class of the ordering wavevector (K) are pushed down by the field and the others
# (K′) are pushed up. That asymmetry is what opens the decay channel K′ → K + K, and it
# is the whole point of the paper. Being a single-𝐐 state it is still one that the cubic
# vertex cannot connect across branches, so — like the spiral — it probes only the
# diagonal of Σ̂; the off-diagonal elements need inequivalent sublattices.
#
# Conventions. The paper's Eq. (1) is H = J Σ_⟨ij⟩ (SˣSˣ + SʸSʸ + Δ SᶻSᶻ) − H Σ Sᶻ with
# saturation field Hs = 6Js(Δ + 1/2), and its ordering wavevector is 𝐐 = (4π/3, 0), which
# is [2/3, -1/3, 0] in hexagonal r.l.u. Sunny's Zeeman energy is +𝐁·(g𝐒), so with g = 1 a
# field along -ẑ is what cants the moments toward +ẑ. Classical minimization gives the
# canting angle sinθ = H/Hs, which `xxz_umbrella` imposes directly; `minimize_energy!`
# from random spins confirms it is the global minimum (unlike the isotropic model in a
# field, the easy-plane anisotropy leaves no classical degeneracy to fall into).
#
# The regulator η is the paper's artificial width δ = 0.005 J, half of the 2δ = 0.01 J
# bar drawn in the inset of their Fig. 3(a). It is small, and the loop grid grows as 1/η
# in each dispersing direction, so this is a 555×555 grid per wavevector — the dominant
# cost, and the reason for `threaded=true`.

using Sunny, LinearAlgebra, Printf, Statistics
using CairoMakie
CairoMakie.activate!(type="png")

# See the comment in fig4b.jl: Julia's OpenBLAS anti-scales when `bogoliubov!` is called
# from many threads, and Apple's Accelerate does not.
BLAS.set_num_threads(1)
Sys.isapple() && @eval using AppleAccelerate

const DIR = @__DIR__

# The paper's labeled wavevectors, in r.l.u. of the one-site hexagonal cell. K is the
# ordering wavevector (4π/3, 0); K′ is the adjacent corner (2π/3, 2π/√3); M = (π, π/√3)
# is the midpoint of the zone edge that joins them.
const Kpt = [2/3, -1/3, 0]
const Kppt = [1/3, 1/3, 0]
const Mpt = [1/2, 0, 0]

# Umbrella state of the easy-plane XXZ triangular antiferromagnet, in the three-site
# magnetic cell. The winding is taken along 𝐐 = K rather than -K, which selects the
# chirality domain in which K is the corner that the field pushes down, as in the paper.
function xxz_umbrella(; s=1/2, Δ=0.9, hfrac=0.2)
    cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(; s, g=1)], :dipole)
    set_exchange!(sys, diagm([1.0, 1.0, Δ]), Bond(1, 1, [1, 0, 0]))
    Hs = 6s * (Δ + 1/2)
    set_field!(sys, [0, 0, -hfrac * Hs])
    sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
    θ = asin(hfrac)
    Q = cryst.recipvecs * Kpt
    for site in eachsite(sys)
        φ = dot(Q, global_positions(sys)[site])
        set_dipole!(sys, [cos(θ)cos(φ), cos(θ)sin(φ), sin(θ)], site)
    end
    return (; swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false)), cryst, θ, Hs)
end

function xxz_fig3(; s=1/2, Δ=0.9, hfrac=0.2, η=0.005, npath=400, tol=0.01, ωmax=2.0,
                    cmax=21.0, file="$DIR/xxz_fig3.png")
    (; swt, cryst, θ, Hs) = xxz_umbrella(; s, Δ, hfrac)

    # Harmonic cross-check against their Fig. 1. Folded into the three-site cell, the
    # bands at Γ are the extended-zone ε_𝐤 at 𝐤 = 0, ±K, i.e. the Goldstone mode
    # together with the two corner energies, which the field splits.
    ε0 = sort(vec(dispersion(swt, [[0, 0, 0]])))
    @printf("harmonic: ε(Γ) = %.4f, ε(K) = %.4f, ε(K′) = %.4f  (their Fig. 1: 0, 0.27, 1.18)\n",
            ε0[1], ε0[2], ε0[3])
    @printf("umbrella: θ = %.2f°, Hs = %.3f J, H = %.3f J\n", rad2deg(θ), Hs, hfrac*Hs)

    path = q_space_path(cryst, [Mpt, Kppt, [0, 0, 0], Kpt, Mpt, [0, 0, 0]], npath;
                        labels=["M", "K′", "Γ", "K", "M", "Γ"])
    energies = range(0, ωmax, 1 + ceil(Int, 2ωmax/η))

    t = @elapsed res = Sunny.intensities_corrected(swt, path; energies, η, tol, threaded=true, verbose=true)
    @printf("%.0f s for %d 𝐪 × %d ω on %d threads\n", t, npath, length(energies), Threads.nthreads())

    fig = Figure(size=(1000, 500), fontsize=15)
    # Their colour scale is a rainbow from zero, cut off at ≈21 so that the Bragg
    # divergence at the ordering wavevector is allowed to saturate.
    plot_intensities!(fig[1, 1], res; colormap=:jet, colorrange=(0, cmax), axis=(; ylabel="ω / J"),
                      title=@sprintf("XXZ triangular AFM, s = %.1f, Δ = %.2f, H = %.1f Hs, 𝒮ᵗᵒᵗ(𝐪, ω) to order 1/s \
                                      — cf. Maksimov et al., Fig. 3(a)", s, Δ, hfrac))
    @printf("intensity: max %.1f, %.2f%% of pixels clipped at %.1f\n",
            maximum(res.data), 100mean(>(cmax), res.data), cmax)

    # Their inset: the line shape at K′, whose peak they find at ω ≈ 1.0 J with a width
    # "modest compared with" the 2δ = 2η bar. Report both, since the peak position is a
    # renormalization of the harmonic 1.18 and is the sharpest number in the figure. The
    # search is restricted to a window about that harmonic energy, because the tallest
    # feature of this column is elsewhere: K′ also carries the Goldstone mode of the
    # three-site cell, whose weight diverges as ω → 0.
    iq = argmin(iq -> norm(path.qs[iq] - Kppt), eachindex(path.qs))
    cut = res.data[:, iq]
    win = findall(ω -> 0.7 ≤ ω ≤ 1.4, energies)
    (peak, i) = findmax(view(cut, win))
    iω = win[i]
    # The half-maximum crossings are interpolated: `energies` is spaced by η/2, which is
    # coarse compared with a width of a few η and biases a nearest-sample width upward.
    crossing(i, j) = (energies[i]*(cut[j] - peak/2) + energies[j]*(peak/2 - cut[i])) / (cut[j] - cut[i])
    lo = findlast(<(peak/2), view(cut, 1:iω))
    hi = findfirst(<(peak/2), view(cut, iω:length(cut)))
    fwhm = isnothing(lo) || isnothing(hi) ? NaN : crossing(iω+hi-1, iω+hi-2) - crossing(lo, lo+1)
    @printf("K′ line shape: peak %.1f at ω = %.3f J, FWHM %.4f J = 2(η + Γ) with Γ = %.4f J \
             (theirs: ≈24 at ω ≈ 1.00, Γ ≈ 0.0025)\n", peak, energies[iω], fwhm, fwhm/2 - η)

    # Decorations sit outside the white plot box, over the dark background of the map, so
    # they are drawn in white.
    ax = Axis(fig[1, 1]; width=Relative(0.36), height=Relative(0.26), halign=0.45, valign=0.98,
              xlabel="ω / J", xlabelpadding=0, ylabelpadding=0, backgroundcolor=:white,
              xgridvisible=false, ygridvisible=false, xlabelcolor=:white,
              xticklabelcolor=:white, yticklabelcolor=:white, xtickcolor=:white, ytickcolor=:white)
    translate!(ax.blockscene, 0, 0, 100)
    lines!(ax, energies, cut; color=:blue)
    xlims!(ax, energies[iω] - 0.05, energies[iω] + 0.045)
    ylims!(ax, 0, 1.3peak)
    text!(ax, 0.95, 0.9; text="𝐪 = K′", space=:relative, align=(:right, :top))

    save(file, fig; px_per_unit=2)
    println("wrote $file")
    return fig
end

xxz_fig3()
