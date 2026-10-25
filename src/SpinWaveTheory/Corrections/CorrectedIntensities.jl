# Dynamical spin structure factor including the O(1/s) corrections assembled by
# the rest of this module. Conventions, including the Dyson equation solved
# here, are collected in Corrections.jl.
#
# The transverse channel is resummed rather than corrected term by term, which
# is what moves the weight a decaying magnon loses into the continuum where it
# decays; an additive treatment would double count it, once in the unit-area
# line shape of the pole and once in the continuum. Nor is the two-magnon
# channel simply added, a pair being reachable both from the longitudinal part
# of the observable directly and from its transverse part through the cubic
# vertex. The structure factor is therefore
#
#     S(𝐪, ω) = ũ† A ũ + Σ_bins lor(ω - x_bin) [2 Re(ũᵗG √18v β̄) + |β|²],
#
# the spectral function A = -Im G/π of Corrections.jl contracted with the
# corrected observable amplitudes, plus the part of the pair weight that the
# resummation does not already carry. `corrected_channels` returns those as
# three channels: `transverse`, the first term, holding the quasiparticle peak
# together with the weight a decaying magnon sheds into the continuum; `direct`,
# the two-magnon continuum of TwoMagnon.jl; and `cross`, their interference. The
# magnon-magnon block of the binned measure is absent from that formula because
# broadening it to ω gives exactly -Im Σ/π, which the Dyson denominator of A
# already contains; only the columns carrying a direct amplitude are needed
# below.
#
# Resumming only what the truncation keeps analytic fixes both the equation
# solved — the projection onto the particle block, described in Corrections.jl —
# and the treatment of the source channel, frozen on shell as SelfEnergy.jl
# explains. The remaining frequency dependence is then a sum of R/(ω - x + iΓ)
# with R ⪰ 0 and x real, so Im D = ηI - Im Σ ⪰ ηI ≻ 0 and the spectral function
# A = -Im D⁻¹/π is positive semidefinite with ‖A‖ ≤ 1/πη, while D → ωI at large
# ω gives ∫dω A = I: the transverse weight of each 𝐪 stays exactly the static
# Σ_n |ũ_n|² of the corrected observables, at any s and on any grid.
#
# The η above regulates the analytic continuation and is numerical rather than
# physical. Because a retarded function is analytic in the upper half plane,
# evaluating it at ω + iη convolves its spectrum with a Lorentzian of half-width
# η; the shift is applied to the self-energy as well as to the Dyson
# denominator. Its role is to give the Dirac deltas of the loop integrand a
# finite width, so that a finite wavevector grid can resolve them, and both that
# grid and the pair-energy bin width are set relative to it. Instrumental
# resolution is a separate convolution, left to the caller.

"""
    intensities_corrected(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                          mean_field_maxevals=100_000, threaded=false, verbose=false)

Dynamical spin structure factor at temperature ``T = 0``, including corrections
of relative order ``1/s`` to linear spin wave theory. Three things distinguish
the result from [`intensities`](@ref). The magnon poles are shifted in energy,
by the mean fields of [`hartree_fock_correction`](@ref) and
[`tadpole_correction`](@ref) together with the real part of
[`cubic_self_energy`](@ref). They are broadened by minus its imaginary part,
which is to say that magnons able to decay into two magnons have a finite
lifetime, and the weight they lose appears in the continuum into which they
decay. Included as well is the two-magnon continuum of
[`intensities_two_magnon`](@ref), which the observable creates directly,
together with its interference with the continuum a decaying magnon feeds. The
two are not separately observable, a pair of magnons being reachable either way.

The regulator `η`, with units of energy, is required. It gives every Dirac delta
a finite width, so that the momentum integrals below can be performed on a
finite grid: the Green function is evaluated at ``ω + iη``, which by analyticity
is the same as convolving the spectrum with `lorentzian(fwhm=2η)`. Choose it
small compared to the magnon linewidths being calculated, but no smaller,
because the wavevector grid below must resolve it and so grows as `1/η` in each
dimension that disperses. Instrumental resolution is a separate, and usually
larger, broadening; apply it to the result afterwards, computing `energies` over
a range padded beyond the one to be displayed so that the convolution is not
truncated.

Everything else controls convergence, and `tol` sets it all. It is a target for
the relative accuracy of the momentum integrals, of which there are two kinds.
The frequency-dependent ones, the self-energy and the two-magnon continuum, are
performed on a uniform grid of the magnetic Brillouin zone, whose dimensions
follow from `η` and `tol` unless `loop_grid` is given explicitly as a tuple;
here `tol` is a calibration rather than a guarantee, and halving it doubles the
work in two dimensions. The static mean fields are performed instead by adaptive
cubature, which gives up after `mean_field_maxevals` evaluations of the
integrand and warns if `tol` was not reached by then. Spacing the `energies`
more finely than `η` costs almost nothing, and a warning is issued if they are
spaced more coarsely.

Set `threaded=true` to parallelize over `qpts`, and `verbose=true` to print the
selected parameters together with the linewidths that resulted.

The three contributions summed here can be obtained separately from
`Sunny.corrected_channels`, which takes the same arguments and returns
`transverse`, `cross` and `direct`. That is the way to compare against published
calculations, which omit the interference `cross`: their convention is
`transverse + direct`.
"""
function intensities_corrected(swt::SpinWaveTheory, qpts; kwargs...)
    (; cryst, qpts, energies, transverse, cross, direct) = corrected_channels(swt, qpts; kwargs...)
    data = transverse + cross + direct
    return Intensities(cryst, qpts, energies, reshape(data, length(energies), size(qpts.qs)...))
end

# Workhorse of `intensities_corrected`, which returns the sum of the channels
# described above. They are kept apart here because each is separately
# meaningful: `transverse` is the resummed magnon pole, `direct` is what
# `intensities_two_magnon` returns on its own, and `cross` is their
# interference. Each is a matrix over (energy, wavevector). Also returned are
# `specfunc`, the magnon spectral matrix before contraction with observables,
# and `disp`, the harmonic energies.
#
# Keeping them apart is also how `cross` is toggled: it has no counterpart in
# the published 1/s calculations, so `transverse + direct` is their convention
# and the sum that `intensities_corrected` forms is ours. See Corrections.jl.
function corrected_channels(swt::SpinWaveTheory, qpts; energies, η, tol=0.01, loop_grid=nothing,
                            mean_field_maxevals=100_000, threaded=false, verbose=false,
                            spectral=false)
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
    # Discretization of the pair energy, whose error is O((bin_width/η)²). The
    # cap of η/16 is what `intensities_two_magnon` defaults to, and is already
    # negligible.
    bin_width = η * min(1/16, sqrt(tol))

    # Sampling the frequency axis is cheap compared to the wavevector loop,
    # which is shared by every frequency, so there is no reason to undersample
    # the Lorentzian.
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

    chans = (; transverse = zeros(eltype(measure), length(energies), length(qpts.qs)),
               cross = zeros(eltype(measure), length(energies), length(qpts.qs)),
               direct = zeros(eltype(measure), length(energies), length(qpts.qs)))
    # Decay rate of each magnon at its own pole, for the `verbose` report
    linewidths = fill(NaN, L, length(qpts.qs))
    # Magnon spectral matrix A = -Im G_pp/π in the Bogoliubov basis, whose
    # diagonal is the A₁₁ of arXiv:1306.1231, which applies the observable
    # factors outside it. Opt-in, because holding every frequency of every
    # wavevector costs L² times the channels themselves.
    specfunc = spectral ? zeros(ComplexF64, L, L, length(energies), length(qpts.qs)) : nothing
    disp = zeros(Float64, L, length(qpts.qs))

    # Buffers are allocated per wavevector rather than reused, which is what
    # makes the loop below safe to thread. The cost is negligible beside the
    # wavevector loop of `accum_cubic_self_energy!` inside. The loop grid is
    # also rebuilt per wavevector, its offset depending on 𝐪 for the reason
    # `loop_wavevectors` explains.
    function calc_iq!(iq)
        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        δH = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Nobs)
        Σ3 = zeros(ComplexF64, L, L, length(energies))
        pref = zeros(ComplexF64, Nobs, L)
        corr = zeros(ComplexF64, num_correlations(measure))

        q = qpts.qs[iq]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        ε = excitations!(T, H, swt, q)
        view(disp, :, iq) .= view(ε, 1:L)

        accum_quadratic!(fill!(δH, 0), terms2, q_reshaped)
        # Static self-energy of Corrections.jl, Σ̂ = Ĩ (T†δH T)ᵗ, restricted to
        # the particle block that the Dyson equation below solves. The metric Ĩ
        # is the identity there, so only the transpose survives.
        Σstat = transpose(view(T' * δH * T, 1:L, 1:L))
        # Frequencies at which to freeze the source channel. Averaging the two
        # external legs keeps the frozen contribution Hermitian, and reduces on the
        # diagonal to the ω = ε_𝐪n of arXiv:1306.1231; the off-diagonal choice is an
        # ambiguity of relative order 1/s.
        onshell = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
        ps = loop_wavevectors(loop_grid, q_reshaped)
        # One pass over the loop wavevectors serves the self-energy and both parts of
        # the pair amplitude, which must share Bogoliubov matrices; see Corrections.jl.
        pair_amplitude_prefactors!(pref, swt, q_reshaped, q_global)
        ρ = Matrix{ComplexF64}[]
        Σsrc = accum_pair_measure!(ρ, swt, terms3, q_reshaped, ps; source_freqs=onshell, bin_width, pref)
        pair_self_energy!(Σ3, ρ, Σsrc, energies .+ im*η, bin_width)

        # Conjugated amplitudes conj(ũ) = T† u, including the 1/s correction to the
        # observables themselves. Their harmonic part, conj(ũ[n, μ]), is the
        # amplitude that `intensities_bands` squares.
        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
        accum_observable_corrections!(u, swt, q_reshaped, q_global, δc)
        w = T' * u

        # Pair-creation blocks of the binned measure after broadening to a
        # single frequency: the interference of the magnon-mediated amplitude
        # with the direct one, and the direct one squared. Only the columns that
        # carry a direct amplitude are needed, the magnon-magnon block of the
        # measure already being resummed into Σ3 and so reaching the transverse
        # channel through G.
        Mω = zeros(ComplexF64, L + Nobs, Nobs)
        md = view(Mω, 1:L, :)
        dd = view(Mω, L+1:L+Nobs, :)
        # Bins the measure actually occupies; the rest contribute nothing at any
        # ω.
        live = findall(!iszero, ρ)

        # Contracts an amplitude product over observable pairs into one channel
        function accum_channel!(accum, iω, f)
            map!(f, corr, measure.corr_pairs)
            accum[iω, iq] = measure.combiner(q_global, corr)
        end

        for (iω, ω) in enumerate(energies)
            # Dyson equation for the block that propagates physical magnons. The
            # metric Ĩ is the identity there, so it does not appear.
            G = inv((ω + im*η)*I - Diagonal(view(ε, 1:L)) - Σstat - view(Σ3, :, :, iω))
            # Magnon spectral matrix A = -Im G/π, i.e. the anti-Hermitian part
            # of G. Contracting it with the observable amplitudes gives the
            # whole transverse channel at once: the quasiparticle peak and the
            # weight a decaying magnon sheds into the continuum are the same
            # resummed pole, not two terms.
            A = (G' - G) ./ (2im * π)
            isnothing(specfunc) || (view(specfunc, :, :, iω, iq) .= A)
            # Amplitude for observable μ to create a magnon of band n which then
            # propagates to frequency ω, conj(z[n, μ]) = (ũ_μᵗ G)[n]. The
            # interference below pairs it against the direct pair-creation
            # amplitude.
            z = G' * view(w, 1:L, :)

            fill!(Mω, 0)
            for bin in live
                Mω .+= ((η/π) / ((ω - (bin - 1) * bin_width)^2 + η^2)) .* view(ρ[bin], :, L+1:L+Nobs)
            end

            accum_channel!(chans.transverse, iω, ((μ, ν),) ->
                dot(view(w, 1:L, μ), A, view(w, 1:L, ν)) / Ncells)
            # The cross term of the squared pair amplitude, read plainly.
            accum_channel!(chans.cross, iω, ((μ, ν),) ->
                (dot(view(z, :, μ), view(md, :, ν)) +
                 conj(dot(view(z, :, ν), view(md, :, μ)))) / Ncells)
            accum_channel!(chans.direct, iω, ((μ, ν),) -> dd[μ, ν] / Ncells)
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
        # A regulator much larger than the calculated linewidths is dominating
        # the line shapes, and one much smaller than them is being paid for
        # needlessly. The upper figure is a quantile rather than the maximum,
        # which is set by the divergence of the vertices at an ordering
        # wavevector and says nothing about the rest.
        Γs = sort!(filter(isfinite, vec(linewidths)))
        r2 = x -> round(x; sigdigits=2)
        report = isempty(Γs) ? "none of the LSWT poles lie within `energies`" :
            "median $(r2(Γs[cld(end, 2)])), 90th pct $(r2(Γs[ceil(Int, 0.9end)])), against η = $(r2(η))"
        println("""
            intensities_corrected with tol = $tol
              loop grid       $(join(loop_grid, "×")) = $(prod(loop_grid)) wavevectors
              bin width       $(r2(bin_width))
              on-shell -Im Σ  $report""")
    end

    return (; cryst, qpts, energies, chans..., specfunc, disp)
end

