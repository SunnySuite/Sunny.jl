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
#   * Im D = ΓI - Im Σ ⪰ ΓI ≻ 0, so D is nonsingular for every real frequency, A is
#     positive semidefinite, and ‖A‖ ≤ 1/(πΓ). No feature can be sharper or taller
#     than the instrumental resolution allows.
#   * D → ωI at large frequency, so ∫dω A = I, and the transverse weight of each 𝐪
#     is exactly the static weight Σ_n |ũ_n|² of the corrected observables. Weight is
#     conserved identically, rather than up to the order worked to.
#
# Neither survives the full Nambu inversion at s = 1/2: the discarded blocks carry
# poles at ω = -ε_{-𝐪n}, which a correction comparable to ε can push up through zero,
# and the near-singular direction then reaches the particle block through the
# anomalous blocks of the self-energy.
#
# The instrumental resolution enters as the single broadening parameter. Because a
# retarded function is analytic in the upper half plane, convolving the spectrum
# with a Lorentzian of half-width Γ is the same as evaluating it at ω + iΓ, and
# that shift is applied to the self-energy as well as to the Dyson denominator. It
# doubles as the regulator of the loop integral, whose wavevector grid must
# therefore resolve the instrumental width.

"""
    intensities_corrected(swt::SpinWaveTheory, qpts; energies, kernel, grid, opts...)

Dynamical spin structure factor at temperature ``T = 0``, including corrections of
relative order ``1/s`` to linear spin wave theory. Three things distinguish the
result from [`intensities`](@ref). The magnon poles are shifted in energy, by the
mean fields of [`hartree_fock_correction`](@ref) and [`tadpole_correction`](@ref)
together with the real part of [`cubic_self_energy`](@ref). They are broadened by
minus its imaginary part, which is to say that magnons able to decay into two
magnons have a finite lifetime, and the weight they lose appears in the continuum
into which they decay. Added to this is the longitudinal two-magnon continuum of
[`intensities_two_magnon`](@ref).

The `kernel` must be `lorentzian(; fwhm)`, representing instrumental resolution.
It is applied by analytic continuation rather than by explicit convolution, so a
magnon that cannot decay appears with the resolution width alone. It also
regularizes the momentum-space integral of the self-energy.

The two frequency-dependent corrections, the self-energy and the two-magnon
continuum, are both integrated over the magnetic Brillouin zone on a uniform
`grid` of the given dimensions. That grid must be fine enough to resolve `fwhm`;
see [`cubic_self_energy`](@ref), and note that a linewidth is only meaningful
once the grid is converged.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of the momentum-space integrals of the remaining, static corrections.
"""
function intensities_corrected(swt::SpinWaveTheory, qpts; energies, kernel::AbstractBroadening, grid, opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")
    check_corrections_supported(swt)
    isa(kernel, Broadening) && !isnan(kernel.fwhm) && kernel.fwhm > 0 &&
        kernel(0.0, 1.0) ≈ lorentzian(; fwhm=kernel.fwhm)(0.0, 1.0) ||
        error("Keyword `kernel` must be `lorentzian(; fwhm)`, which the Dyson equation applies by analytic continuation.")

    (; sys, measure) = swt
    cryst = orig_crystal(sys)
    L = nbands(swt)
    Nobs = num_observables(measure)
    # Number of chemical cells in the magnetic cell
    Ncells = nsites(sys) / natoms(cryst)
    Γ = kernel.fwhm / 2

    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")
    qpts = convert(AbstractQPoints, qpts)

    tad = tadpole_correction(swt; opts...)
    terms2 = [hartree_fock_correction(swt; opts...).terms2
              tad.terms2
              anisotropy_correction(swt).terms2]
    δc = observable_corrections(swt; v=tad.v, opts...)
    terms3 = cubic_monomials(swt)
    ps = loop_grid(grid)

    # The longitudinal channel, to which the transverse one is added below
    ret = intensities_two_magnon(swt, qpts; energies, kernel, grid).data

    Ĩ = Diagonal([ones(L); -ones(L)])
    T = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    δH = zeros(ComplexF64, 2L, 2L)
    u = zeros(ComplexF64, 2L, Nobs)
    Σ3 = zeros(ComplexF64, L, L, length(energies))
    corr = zeros(ComplexF64, num_correlations(measure))

    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        ε = excitations!(T, H, swt, q)

        accum_quadratic!(fill!(δH, 0), terms2, q_reshaped)
        Σstat = Ĩ * transpose(T' * δH * T)
        # Frequencies at which to freeze the source channel. Averaging the two
        # external legs keeps the frozen contribution Hermitian, and reduces on the
        # diagonal to the ω = ε_𝐪n of Mourigal et al.; the off-diagonal choice is an
        # ambiguity of relative order 1/s.
        onshell = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
        # The bin width fwhm/32 is the one `intensities_two_magnon` also defaults to
        accum_cubic_self_energy!(fill!(Σ3, 0), swt, terms3, q_reshaped, energies .+ im*Γ, ps, 0.0; source_freqs=onshell, bin_width=Γ/16)

        # Conjugated amplitudes conj(ũ) = T† u, including the 1/s correction to the
        # observables themselves. Their harmonic part, conj(ũ[n, μ]), is the
        # amplitude that `intensities_bands` squares.
        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
        accum_observable_corrections!(u, swt, q_reshaped, q_global, δc)
        w = T' * u

        for (iω, ω) in enumerate(energies)
            # Dyson equation for the block that propagates physical magnons. The
            # metric Ĩ is the identity there, so it does not appear.
            G = inv((ω + im*Γ)*I - Diagonal(view(ε, 1:L)) - view(Σstat, 1:L, 1:L) - view(Σ3, :, :, iω))
            A = (G - G') / 2im
            map!(corr, measure.corr_pairs) do (μ, ν)
                -dot(view(w, 1:L, μ), A, view(w, 1:L, ν)) / (π * Ncells)
            end
            ret[iω, iq] += measure.combiner(q_global, corr)
        end
    end

    return Intensities(cryst, qpts, energies, ret)
end

