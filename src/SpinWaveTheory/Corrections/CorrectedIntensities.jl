# Dynamical spin structure factor including the O(1/s) corrections assembled by
# the rest of this module. Conventions, including the Dyson equation solved here,
# are collected in Corrections.jl.
#
# The transverse channel is resummed rather than corrected term by term. Solving
# the Dyson equation at each frequency shifts the magnon poles by the static mean
# fields and the real part of the cubic self-energy, gives them the width supplied
# by its imaginary part, and — because the same self-energy appears in the
# resummation — moves the weight that a decaying magnon loses into the continuum
# where it decays. An additive treatment would instead double count it, once in the
# unit-area line shape of the pole and once in the continuum.
#
# The longitudinal channel of TwoMagnon.jl is a separate observable, even in the
# boson number rather than odd, so it is still added. Its weight is of order s⁰
# already, and correcting it would be a higher-order calculation.
#
# Only what the truncation can keep analytic is resummed, which fixes both the
# equation solved and the treatment of the source channel. The equation is projected
# onto the L×L particle block, as in Eq. (12) of Mourigal et al., rather than
# inverting the full Nambu denominator and keeping the particle block of the
# solution; the source channel of the self-energy is frozen on shell, as SelfEnergy.jl
# explains. What remains of the frequency dependence is then a sum of terms
# R/(ω - x + iΓ) with R ⪰ 0 and x real, whose imaginary part is negative
# semidefinite, and the resulting spectral function A = -Im D⁻¹/π satisfies two
# properties exactly, at any s and on any wavevector grid:
#
#   * Im D = ηI - Im Σ ⪰ ηI ≻ 0, so D is nonsingular for every real frequency, A is
#     positive semidefinite, and ‖A‖ ≤ 1/(πη). No feature can be sharper or taller
#     than the regulator allows.
#   * D → ωI at large frequency, so ∫dω A = I, and the transverse weight of each 𝐪
#     is exactly the static weight Σ_n |ũ_n|² of the corrected observables. Weight is
#     conserved identically, rather than up to the order worked to.
#
# Neither survives the full Nambu inversion at s = 1/2: the discarded blocks carry
# poles at ω = -ε_{-𝐪n}, which a correction comparable to ε can push up through zero,
# and the near-singular direction then reaches the particle block through the
# anomalous blocks of the self-energy.
#
# The η above is the regulator of the analytic continuation, and it is a numerical
# parameter rather than a physical one. Because a retarded function is analytic in the
# upper half plane, evaluating it at ω + iη is the same as convolving its spectrum
# with a Lorentzian of half-width η, and the shift is applied to the self-energy as
# well as to the Dyson denominator. Its role is to give the Dirac deltas of the loop
# integrand a finite width, so that a finite wavevector grid can resolve them; the
# grid and the bin width of the pair energy are then both set relative to η.
# Instrumental resolution is a separate convolution, left to the caller.

"""
    intensities_corrected(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                          mean_field_maxevals=100_000, threaded=false, verbose=false)

Dynamical spin structure factor at temperature ``T = 0``, including corrections of
relative order ``1/s`` to linear spin wave theory. Three things distinguish the
result from [`intensities`](@ref). The magnon poles are shifted in energy, by the
mean fields of [`hartree_fock_correction`](@ref) and [`tadpole_correction`](@ref)
together with the real part of [`cubic_self_energy`](@ref). They are broadened by
minus its imaginary part, which is to say that magnons able to decay into two
magnons have a finite lifetime, and the weight they lose appears in the continuum
into which they decay. Added to this is the longitudinal two-magnon continuum of
[`intensities_two_magnon`](@ref).

The regulator `η`, with units of energy, is required. It gives every Dirac delta a
finite width, so that the momentum integrals below can be performed on a finite
grid: the Green function is evaluated at ``ω + iη``, which by analyticity is the
same as convolving the spectrum with `lorentzian(fwhm=2η)`. Choose it small
compared to the magnon linewidths being calculated, but no smaller, because the
wavevector grid below must resolve it and so grows as `1/η` in each dimension that
disperses. Instrumental resolution is a separate, and usually larger, broadening;
apply it to the result afterwards, computing `energies` over a range padded beyond
the one to be displayed so that the convolution is not truncated.

Everything else controls convergence, and `tol` sets it all. It is a target for the
relative accuracy of the momentum integrals, of which there are two kinds. The
frequency-dependent ones, the self-energy and the two-magnon continuum, are
performed on a uniform grid of the magnetic Brillouin zone, whose dimensions follow
from `η` and `tol` unless `loop_grid` is given explicitly as a tuple; here `tol` is
a calibration rather than a guarantee, and halving it doubles the work in two
dimensions. The static mean fields are performed instead by adaptive cubature, which
gives up after `mean_field_maxevals` evaluations of the integrand and warns if `tol`
was not reached by then. Spacing the `energies` more finely than `η` costs almost
nothing, and a warning is issued if they are spaced more coarsely.

Set `threaded=true` to parallelize over `qpts`, and `verbose=true` to print the
selected parameters together with the linewidths that resulted.
"""
function intensities_corrected(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                               mean_field_maxevals=100_000, threaded=false, verbose=false)
    check_corrections_supported(swt)
    η > 0 || error("Regulator `η` must be positive.")

    (; sys, measure) = swt
    cryst = orig_crystal(sys)
    L = nbands(swt)
    Nobs = num_observables(measure)
    # Number of chemical cells in the magnetic cell
    Ncells = nsites(sys) / natoms(cryst)

    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")
    qpts = convert(AbstractQPoints, qpts)

    loop_grid = @something loop_grid auto_loop_grid(swt, η, tol)
    ps = loop_wavevectors(loop_grid)
    # Discretization of the pair energy, whose error is O((bin_width/η)²). The cap of
    # η/16 is what `intensities_two_magnon` defaults to, and is already negligible.
    bin_width = η * min(1/16, sqrt(tol))

    # Sampling the frequency axis is cheap compared to the wavevector loop, which is
    # shared by every frequency, so there is no reason to undersample the Lorentzian.
    if length(energies) > 1
        dω = (energies[end] - energies[begin]) / (length(energies) - 1)
        dω > η/2 && @warn """Requested `energies` are spaced by $(round(dω, sigdigits=2)) on \
                             average, which will not resolve features of width η = $η. A spacing \
                             of η/2 or less adds little cost."""
    end

    # Self-consistent HF can be enabled with maxiters > 1. Note, however, that
    # this resummation is an uncontrolled approximation, e.g. violates Ward
    # identity and can gap Goldstone modes.
    #
    #   hartree_fock_correction(swt; maxiters=100, damping=0.5, rtol=tol, ...)

    tad = tadpole_correction(swt; rtol=tol, maxevals=mean_field_maxevals)
    terms2 = [hartree_fock_correction(swt; maxiters=1, rtol=tol, maxevals=mean_field_maxevals).terms2
              tad.terms2
              anisotropy_correction(swt).terms2]
    δc = observable_corrections(swt; v=tad.v, rtol=tol, maxevals=mean_field_maxevals)
    terms3 = cubic_monomials(swt)

    Ĩ = Diagonal([ones(L); -ones(L)])
    ret = zeros(eltype(measure), length(energies), length(qpts.qs))
    # Decay rate of each magnon at its own pole, for the `verbose` report
    linewidths = fill(NaN, L, length(qpts.qs))

    # Buffers are allocated per wavevector rather than reused, which is what makes the
    # loop below safe to thread. The cost is negligible beside the wavevector loop of
    # `accum_cubic_self_energy!` inside.
    function calc_iq!(iq)
        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        δH = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Nobs)
        Σ3 = zeros(ComplexF64, L, L, length(energies))
        corr = zeros(ComplexF64, num_correlations(measure))

        q = qpts.qs[iq]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        ε = excitations!(T, H, swt, q)

        # The longitudinal channel, to which the transverse one is added below
        res2 = intensities_two_magnon(swt, [q]; energies, kernel=lorentzian(fwhm=2η), grid=loop_grid, bin_width)
        view(ret, :, iq) .= vec(res2.data)

        accum_quadratic!(fill!(δH, 0), terms2, q_reshaped)
        Σstat = Ĩ * transpose(T' * δH * T)
        # Frequencies at which to freeze the source channel. Averaging the two
        # external legs keeps the frozen contribution Hermitian, and reduces on the
        # diagonal to the ω = ε_𝐪n of Mourigal et al.; the off-diagonal choice is an
        # ambiguity of relative order 1/s.
        onshell = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
        accum_cubic_self_energy!(Σ3, swt, terms3, q_reshaped, energies .+ im*η, ps, 0.0; source_freqs=onshell, bin_width)

        # Conjugated amplitudes conj(ũ) = T† u, including the 1/s correction to the
        # observables themselves. Their harmonic part, conj(ũ[n, μ]), is the
        # amplitude that `intensities_bands` squares.
        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
        accum_observable_corrections!(u, swt, q_reshaped, q_global, δc)
        w = T' * u

        for (iω, ω) in enumerate(energies)
            # Dyson equation for the block that propagates physical magnons. The
            # metric Ĩ is the identity there, so it does not appear.
            G = inv((ω + im*η)*I - Diagonal(view(ε, 1:L)) - view(Σstat, 1:L, 1:L) - view(Σ3, :, :, iω))
            A = (G - G') / 2im
            map!(corr, measure.corr_pairs) do (μ, ν)
                -dot(view(w, 1:L, μ), A, view(w, 1:L, ν)) / (π * Ncells)
            end
            ret[iω, iq] += measure.combiner(q_global, corr)
        end

        for n in 1:L
            if energies[begin] ≤ ε[n] ≤ energies[end]
                iω = argmin(iω -> abs(energies[iω] - ε[n]), eachindex(energies))
                linewidths[n, iq] = -imag(Σ3[n, n, iω])
            end
        end
    end

    if threaded
        Threads.@threads for iq in eachindex(qpts.qs)
            calc_iq!(iq)
        end
    else
        for iq in eachindex(qpts.qs)
            calc_iq!(iq)
        end
    end

    if verbose
        # A regulator much larger than the calculated linewidths is dominating the line
        # shapes, and one much smaller than them is being paid for needlessly. The upper
        # figure is a quantile rather than the maximum, which is set by the divergence of
        # the vertices at an ordering wavevector and says nothing about the rest.
        Γs = sort!(filter(isfinite, vec(linewidths)))
        r2 = x -> round(x; sigdigits=2)
        report = isempty(Γs) ? "none of the LSWT poles lie within `energies`" :
            "median $(r2(Γs[cld(end, 2)])), 90th pct $(r2(Γs[ceil(Int, 0.9end)])), against η = $(r2(η))"
        println("""
            intensities_corrected with tol = $tol
              loop grid       $(join(loop_grid, "×")) = $(length(ps)) wavevectors
              bin width       $(r2(bin_width))
              on-shell -Im Σ  $report""")
    end

    return Intensities(cryst, qpts, energies, reshape(ret, length(energies), size(qpts.qs)...))
end

