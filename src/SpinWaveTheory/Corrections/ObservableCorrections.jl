# Corrections of relative order 1/s to the amplitudes with which an observable
# creates one magnon or a pair. Its boson words are those of
# Bosons/ObservableWords.jl, and conventions are collected in Corrections.jl.
#
# In a local frame where the classical dipole points along ẑ, an observable
# Â_i = 𝐎 ⋅ 𝐒_i splits into a transverse part and a longitudinal one,
#
#     Â_i = c₊ S⁺_i + c₋ S⁻_i + O₃ (s - b†b),   c± = (O₁ ∓ i O₂)/2,
#
# which are respectively odd and even in the boson number. Only the odd part can
# create a single magnon and only the even part can create a pair, so at leading
# order the one-magnon bands of `intensities_bands` come entirely from the
# transverse components and the direct two-magnon continuum entirely from the
# longitudinal one.
#
# LSWT keeps only the leading part of S⁺ = σ(b - b†bb/4s). Retaining the cubic
# word costs nothing at leading order, because it changes the boson number by
# one just as b does, and its two Wick contractions
#
#     b†bb → 2⟨b†b⟩ b + ⟨bb⟩ b†
#
# leave a one-boson operator behind. That is a correction of relative order 1/s
# to the amplitude for creating one magnon. What remains of the cubic word
# creates three magnons, whose weight is smaller by 1/s again.
#
# A second correction of the same order comes from the tadpole. Zero-point
# fluctuations displace the bosons by b → b + v, of order s^(-1/2), and the
# longitudinal word then contributes a one-boson term of its own,
#
#     O₃ (s - b†b) → O₃ (s - |v|² - b†b) - O₃ (v̄ b + v b†),
#
# which is smaller than the transverse amplitude σ by 1/s, again the target
# order. It is the tilt of the ordered moment, in the displaced description of
# the tadpole in StaticCorrections.jl.
#
# The first of these corrections supplies the entire O(1/s) deficit in the
# static transverse weight, and hence in the quantum sum rule; the second is
# invisible to a trace measure, because summing O₃ c₊ over three orthonormal
# observable directions gives zero. By completeness Σ_f ⟨0|Â†|f⟩⟨f|Â|0⟩ = ⟨Â†Â⟩,
# so the sum of the one-magnon weights over bands and wavevectors is a static
# expectation value, and the truncated substitution above makes S⁻S⁺ =
# n̂(2s+1-n̂) exact. Writing n̂ = b†b, the transverse and longitudinal weights
# of one site are therefore
#
#     ⟨(Sˣ)² + (Sʸ)²⟩ = s + 2s⟨n̂⟩ - ⟨n̂²⟩,   ⟨(Sᶻ)²⟩ = (s - ⟨n̂⟩)² + Var(n̂),
#
# which sum to s(s+1) identically, for any state. The two terms of order s
# cancel between the transverse weight and the ordered moment, which is why LSWT
# already satisfies the sum rule to relative order 1/s; the entire content of
# the next order is the -⟨n̂²⟩ above, and the contraction of the cubic word
# produces it.

# Correction of relative order 1/s to the amplitude for each observable to
# create a pair of magnons, as two-boson words of each observable μ and site i,
# without the Fourier prefactor: one contraction of its four-boson words. It is
# the pair analogue of `observable_corrections`. It vanishes for a spin
# component, which has no four-boson word, and in mode :SUN.
#
# Its effect on the weight is beyond the order kept, and is only part of that
# order. For an observable that also creates a single magnon, the pair weight is
# itself one order below the one-magnon weight. For one that does not, e.g. a
# quadrupole that changes Sᶻ by two in the local frame, the pair weight is
# leading and exact, and this is part of its first correction. The rest is the
# shift of ⟨W₂†W₂⟩ that the vacuum takes at first order: a Gaussian part, from
# the mean field of H₄ on the two lines of the pair; a four-magnon admixture,
# from H₄ coupling the pair to two further magnons; and in a canted state the
# second order of H₃ and the tilt of the pair amplitude, one displaced leg of
# the three-boson words. The part kept here costs nothing and reduces the error.
# For 𝒪₂² on a collinear Néel dimer, as in the sum-rule test of the derivation
# tier, the first correction is -0.093, -0.035 and -0.016 of the leading weight
# at s = 2, 4 and 8, of which this gives -0.140, -0.056 and -0.025, leaving an
# error of about +0.08/s.
function observable_pair_corrections(swt::SpinWaveTheory, g)
    return [wick_reduce(observable_monomials(swt, μ, i, Val{4}()), g, nothing, Val{2}())
            for μ in 1:num_observables(swt.measure), i in 1:nsites(swt.sys)]
end

# The two-boson words of each observable μ and site i, as a matrix over (μ, i)
observable_words2(swt::SpinWaveTheory) =
    [observable_monomials(swt, μ, i, Val{2}()) for μ in 1:num_observables(swt.measure), i in 1:nsites(swt.sys)]

# The even words of each observable at one wavevector, as
# [`pair_amplitude`](@ref) consumes them: one list per observable, carrying the
# Fourier phase and form factor of its site. The words `words2[μ, i]` of
# observable μ at site i are by default those of `observable_words2`, and may
# include the corrections of `observable_pair_corrections`. The prefactor is
# conjugated because the amplitude sought is that of creating a pair, whereas
# `observable_prefactor` describes the observable as `observable_vectors!`
# applies it, to the amplitude for the adjoint process.
function observable_pair_monomials(swt::SpinWaveTheory, q_global; words2=observable_words2(swt))
    (; sys, measure) = swt
    return map(1:num_observables(measure)) do μ
        [BosonMonomial(conj(observable_prefactor(measure, μ, i, q_global, sys)) * c, as, ns)
         for i in 1:nsites(sys) for (; c, as, ns) in words2[μ, i]]
    end
end

# Amplitude for an observable to create the pair of magnons (𝐩 a, 𝐪-𝐩 b)
# directly, given its even part as the `BosonMonomial{2}` list of
# `observable_pair_monomials`, already carrying its Fourier phase. The two terms
# symmetrize over which line takes which slot of the word, and the 1/√2 is the
# norm of the symmetrized two-boson state. An observable is onsite, so both
# slots carry the same cell offset and no phase.
#
# Both lines are expressed in the matrices T(𝐩) and T(𝐪-𝐩) that the cubic
# vertex shares, as the interference with it requires. The sign of the word,
# that of Sᶻ = s - b†b in dipole mode, cancels in the |β|² of the direct channel
# but is the whole sign of the interference.
function pair_amplitude(words, T1, T2, a, b)
    return sum(words; init=zero(ComplexF64)) do (; c, as)
        c * (T1[as[2], a]*T2[as[1], b] + T1[as[1], a]*T2[as[2], b])
    end / √2
end

"""
    accum_observable_corrections!(u, swt::SpinWaveTheory, q_global, δc)

Adds the coefficients `δc` of [`observable_corrections`](@ref) to the linearized
observables `u` of `observable_vectors!`, which stores the coefficient of
the boson operator labeled `a` in `u[ā, μ]`, where `ā = mod1(a+L, 2L)`.
"""
function accum_observable_corrections!(u, swt::SpinWaveTheory, q_global, δc)
    return accum_observable_vectors!(u, swt.sys, swt.measure, reshape(δc, size(δc)..., 1), q_global)
end
