# Corrections of relative order 1/s to the observables whose correlations
# `intensities` measures. Conventions are collected in Corrections.jl.
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
# LSWT keeps only the leading part of S⁺ = σ(b - b†bb/4s). Retaining the cubic word
# costs nothing at leading order, because it changes the boson number by one just
# as b does, and its two Wick contractions
#
#     b†bb → 2⟨b†b⟩ b + ⟨bb⟩ b†
#
# leave a one-boson operator behind. That is a correction of relative order 1/s to
# the amplitude for creating one magnon. What remains of the cubic word creates
# three magnons, whose weight is smaller by 1/s again.
#
# A second correction of the same order comes from the tadpole. Zero-point
# fluctuations displace the bosons by b → b + v, of order s^(-1/2), and the
# longitudinal word then contributes a one-boson term of its own,
#
#     O₃ (s - b†b) → O₃ (s - |v|² - b†b) - O₃ (v̄ b + v b†),
#
# which is smaller than the transverse amplitude σ by 1/s, again the target order.
# Equivalently, the displacement tilts the local frame, and re-expanding the
# transverse components about the tilted axis shifts their amplitude by the same
# amount: with a tilt (θ₁, θ₂) = σ(Re v, Im v)/s away from ẑ, the rotated
# components are O₁ - θ₁O₃ and O₂ - θ₂O₃, so
#
#     σ c₊' = σ c₊ - σ O₃ (θ₁ - iθ₂)/2 = σ c₊ - O₃ v̄,
#
# using σ² = 2s. The displaced description is used here because every other part of
# the calculation is performed at the classical energy minimum, where LSWT is well
# defined; see Tadpole.jl.
#
# The first of these corrections supplies the entire O(1/s) deficit in the static
# transverse weight, and hence in the quantum sum rule; the second is invisible to
# a trace measure, because summing O₃ c₊ over three orthonormal observable
# directions gives zero. By completeness
# Σ_f ⟨0|Â†|f⟩⟨f|Â|0⟩ = ⟨Â†Â⟩, so the sum of the one-magnon weights over bands and
# wavevectors is a static expectation value, and the truncated substitution above
# makes S⁻S⁺ = n̂(2s+1-n̂) exact. Writing n̂ = b†b, the transverse and longitudinal
# weights of one site are therefore
#
#     ⟨(Sˣ)² + (Sʸ)²⟩ = s + 2s⟨n̂⟩ - ⟨n̂²⟩,   ⟨(Sᶻ)²⟩ = (s - ⟨n̂⟩)² + Var(n̂),
#
# which sum to s(s+1) identically, for any state. The two terms of order s cancel
# between the transverse weight and the ordered moment, which is why LSWT already
# satisfies the sum rule to relative order 1/s; the entire content of the next
# order is the -⟨n̂²⟩ above, and the contraction of the cubic word produces it.

# Monomials of `K` bosons in the expansion of observable μ at site i, in the local
# frame and without the Fourier phase or form factor that
# `observable_prefactor` supplies. Only K = 2, the longitudinal word, and K = 3,
# the leading correction to the transverse ones, are needed; the one-boson word is
# the LSWT amplitude that `set_swt_observable_vectors!` builds.
#
# In mode :SUN an observable is a matrix on the same footing as a term of the
# Hamiltonian, so its words are those of `local_words`. In the dipole modes the
# expansion above gives the two families directly.
function observable_words(swt::SpinWaveTheory, μ, i, ::Val{K}) where K
    (; sys, data) = swt
    L = nbands(swt)
    o = zero(Vec3)

    if sys.mode == :SUN
        @assert num_parts_per_unit(swt.measure) == 1  # Entangled units are rejected
        A = (data::SWTDataSUN).observables[μ, i, 1]
        return local_words(A, i, o, Val{K}(), nflavors(swt), L)
    end

    @assert sys.mode in (:dipole, :dipole_uncorrected)
    (; sqrtS, observables) = data::SWTDataDipole
    O = observables[μ, i]
    if K == 2
        # Sᶻ = s - b†b
        return [BosonMonomial(ComplexF64(-O[3]), (L+i, i), (o, o))]
    elseif K == 3
        # Both cubic words carry -σ/4s = -1/2σ, that of S⁺ being b†bb and that of
        # S⁻ its adjoint b†b†b.
        σ = √2 * sqrtS[i]
        (cp, cm) = ((O[1] - im*O[2])/2, (O[1] + im*O[2])/2)
        return [BosonMonomial(-cp / 2σ, (L+i, i, i), (o, o, o)),
                BosonMonomial(-cm / 2σ, (L+i, L+i, i), (o, o, o))]
    end
end

# Every K-boson word of observable μ, over the sites of the magnetic cell.
observable_words(swt::SpinWaveTheory, μ, ::Val{K}) where K =
    reduce(vcat, observable_words(swt, μ, i, Val{K}()) for i in 1:nsites(swt.sys))

# The even words of each observable at one wavevector, as [`pair_amplitude`](@ref)
# consumes them: one list per observable, carrying the Fourier phase and form
# factor of its site. The prefactor is conjugated because the amplitude sought is
# that of creating a pair, whereas `observable_prefactor` describes the observable
# as `set_swt_observable_vectors!` applies it, to the amplitude for the adjoint
# process.
function observable_pair_words(swt::SpinWaveTheory, q_reshaped, q_global)
    (; sys, measure) = swt
    return map(1:num_observables(measure)) do μ
        [BosonMonomial(conj(observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)) * c, as, ns)
         for i in 1:nsites(sys) for (; c, as, ns) in observable_words(swt, μ, i, Val{2}())]
    end
end

"""
    observable_corrections(swt::SpinWaveTheory; v=nothing, tol, maxevals)

Correction of relative order ``1/s`` to the amplitude for a magnon to be created
by each observable. Two effects contribute at this order: the cubic term of the
Holstein-Primakoff expansion of the transverse spin components, and the tilt of
the ordered structure by zero-point fluctuations. The latter requires the boson
displacement `v` of [`tadpole_correction`](@ref), and is omitted if `v` is
`nothing`. Returns the coefficients `δc[a, μ]` of the one-boson operators, labeled
as in [`accum_observable_corrections!`](@ref), which is what applies them.

The onsite correlations are integrated over the Brillouin zone by adaptive
cubature. At least one of `tol` (a relative accuracy target) or `maxevals` (a
budget of integrand evaluations) is required to control it.
"""
function observable_corrections(swt::SpinWaveTheory; v=nothing, tol=nothing, maxevals=nothing)
    isnothing(tol) && isnothing(maxevals) && error("Must specify `tol` or `maxevals` to control momentum-space integration.")
    check_corrections_supported(swt)

    L = nbands(swt)
    Nobs = num_observables(swt.measure)
    words2 = [observable_words(swt, μ, Val{2}()) for μ in 1:Nobs]
    words3 = [observable_words(swt, μ, Val{3}()) for μ in 1:Nobs]

    # Contracting two legs of the cubic word is the same operation that generates
    # the tadpole, so `tadpole_vector` performs it. The contracted pair acts on the
    # same site as the surviving operator, so only onsite correlations are needed
    # and no wavevector dependence survives.
    ckeys = correlation_keys(L, reduce(vcat, words3))
    gs = nambu_correlations(swt, ckeys, BosonMonomial{2}[]; tol, maxevals)
    g = correlation_lookup(ckeys, gs, L)
    noise = max(@something(tol, 1e-3), 1e-8)

    # Nambu packing of the displacement, as `tadpole_correction` forms it
    w = isnothing(v) ? zeros(ComplexF64, 2L) : [v; conj(v)]

    δc = zeros(ComplexF64, 2L, Nobs)
    for μ in 1:Nobs
        view(δc, :, μ) .= tadpole_vector(words3[μ], g, L, noise)
        # The tilt, as the displacement of one leg of the longitudinal word.
        for (; c, as) in words2[μ]
            δc[as[1], μ] += c * w[as[2]]
            δc[as[2], μ] += c * w[as[1]]
        end
    end

    return δc
end

"""
    accum_observable_corrections!(u, swt::SpinWaveTheory, q_reshaped, q_global, δc)

Adds the coefficients `δc` of [`observable_corrections`](@ref) to the linearized
observables `u` of `set_swt_observable_vectors!`, which stores the coefficient of
the boson operator labeled `a` in `u[ā, μ]`, where `ā = mod1(a+L, 2L)`.
"""
function accum_observable_corrections!(u, swt::SpinWaveTheory, q_reshaped, q_global, δc)
    (; sys, measure) = swt
    L = nbands(swt)
    Nf = nflavors(swt)
    for μ in 1:num_observables(measure), a in 1:L
        i = boson_site(a, L, Nf)
        pref = observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)
        u[a, μ]   += pref * δc[L+a, μ]
        u[L+a, μ] += pref * δc[a, μ]
    end
    return u
end
