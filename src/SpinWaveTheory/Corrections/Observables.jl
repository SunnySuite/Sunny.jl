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
# transverse components and the two-magnon continuum of TwoMagnon.jl entirely from
# the longitudinal one.
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
# order is the -⟨n̂²⟩ above, and the contraction of the cubic word is what produces
# it. See "1/s correction to the observable amplitudes".

"""
    observable_corrections(swt::SpinWaveTheory; v=nothing, opts...)

Correction of relative order ``1/s`` to the amplitude for a magnon to be created
by each observable. Two effects contribute at this order: the cubic term of the
Holstein-Primakoff expansion of the transverse spin components, and the tilt of
the ordered structure by zero-point fluctuations. The latter requires the boson
displacement `v` of [`tadpole_correction`](@ref), and is omitted if `v` is
`nothing`. Returns the coefficients `δc[a, μ]` of the one-boson operators, labeled
as in [`accum_observable_corrections!`](@ref), which is what applies them.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration.
"""
function observable_corrections(swt::SpinWaveTheory; v=nothing, opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")
    check_corrections_supported(swt)

    (; measure, data) = swt
    (; sqrtS, observables) = data::SWTDataDipole
    L = nbands(swt)
    Nobs = num_observables(measure)

    # The contracted pair acts on the same site as the surviving operator, so only
    # the onsite correlations ⟨b†ᵢbᵢ⟩ and ⟨bᵢbᵢ⟩ are needed, and no wavevector
    # dependence survives.
    ckeys = [[(L+i, i, (0, 0, 0)) for i in 1:L]; [(i, i, (0, 0, 0)) for i in 1:L]]
    gs = nambu_correlations(swt, ckeys, BosonMonomial{2}[]; opts...)

    δc = zeros(ComplexF64, 2L, Nobs)
    for i in 1:L
        (n, Δ) = (gs[i], gs[L+i])
        # Both cubic words carry -σ/4s = -1/2σ, and the adjoint word b†b†b of S⁻
        # contracts to 2⟨b†b⟩ b† + ⟨b†b†⟩ b.
        σ = √2 * sqrtS[i]
        vi = isnothing(v) ? zero(ComplexF64) : v[i]
        for μ in 1:Nobs
            O = observables[μ, i]
            (cp, cm) = ((O[1] - im*O[2])/2, (O[1] + im*O[2])/2)
            δc[i, μ]   = -(2n*cp + conj(Δ)*cm) / 2σ - O[3] * conj(vi)
            δc[L+i, μ] = -(Δ*cp + 2n*cm) / 2σ - O[3] * vi
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
    for μ in 1:num_observables(measure), i in 1:L
        pref = observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)
        u[i, μ]   += pref * δc[L+i, μ]
        u[L+i, μ] += pref * δc[i, μ]
    end
    return u
end
