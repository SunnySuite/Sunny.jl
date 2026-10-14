# Holstein-Primakoff expansion of the Hamiltonian in the dipole modes, as boson
# monomials at every order: H₂ for LSWT, and the higher and sub-leading terms
# for the 1/s corrections. The words of each local operator are those of
# StevensWords.jl, and the monomial representation is that of Monomials.jl.

# In a local frame where the classical dipole points along ẑ, the
# Holstein-Primakoff expansion reads
#
#     Sᶻ = s - b†b,   S⁺ = √(2s) (b - b†bb/4s + …),   S⁻ = (S⁺)†.
#
# Note that Sᶻ is exact, and that S⁺ contains only odd powers of the bosons.
#
# THE TRUNCATION RULE. Every term of the Hamiltonian is expanded in 1/s and each
# n-boson sector Hₙ is kept at its leading order only, with the sole exception of
# H₀, H₁ and H₂, where the first sub-leading order is kept as well because that is
# precisely the O(1/s) correction to LSWT being computed. In particular the cubic
# coefficient -σ/4s above is the leading Taylor coefficient of √(1 - n̂/2s) and
# must not be replaced by σ(√(1-1/2s) - 1), the exact coefficient of the word
# b†bb, even though the latter represents S⁺ better and is exact at s = 1/2.
#
# The rule is what makes continuous symmetries exact, not bookkeeping hygiene. The
# truncated substitution above satisfies [S⁺, S⁻] = 2Sᶻ only up to O(1/s²), so it
# represents su(2) — and hence commutes with the generator of a symmetry the ordered
# structure breaks — only through the order retained. Since the energy of a rotated
# structure is exactly invariant, each coefficient of its 1/s expansion is separately
# invariant, and a calculation keeping every term of one order and none of the next
# inherits that exactly. A partial set of higher-order terms does not: it perturbs the
# quadratic form at a protected zero mode, where a gap grows as the square root of the
# perturbation, so an O(1/s²) error becomes an O(1/s) gap — as large as the correction.
#
# Both failure modes have been measured. The exact b†bb coefficient destroys the O(1/s)
# band shift of the square-lattice antiferromagnet, a thirty-fold cancellation between the
# longitudinal and transverse parts of H₄, and spoils its exact independence of 𝐪.
# Expanding an onsite anisotropy exactly gaps out the Goldstone mode of an easy-plane
# ferromagnet; see the test of `anisotropy_correction`.

# THE PRODUCER. Every term of the Hamiltonian is a product of local operators,
# each a Stevens polynomial on one site, and every local operator expands into
# the words of `stevens_words`, a leading part and a sub-leading part. The
# K-boson term of a product is every split of the K bosons between the words of
# its factors, as in `append_products!`. Each factor carries a single power of
# x = √2s in its leading part, so each product carries a single power as well,
# and the sub-leading part of a product is the sum of the products of one
# sub-leading factor with the leading others. Then the truncation rule of this
# file is that the term of K bosons keeps the leading part only, for K ≥ 3, and
# also the sub-leading part, for K ≤ 2.
#
# The terms are
#
#   * the Zeeman coupling +𝐁⋅𝐒 and the bilinear exchange 𝐒ᵢ⋅J𝐒ⱼ, products of
#     spin components, whose sub-leading words of at most two bosons vanish;
#   * the onsite anisotropy, a single Stevens polynomial;
#   * the biquadratic exchange, the product of the quadrupoles of its two sites.
#
# LSWT is the leading part of the term of two bosons, and is evaluated at the
# coefficients as stored, which in mode :dipole carry the renormalization λ_k of
# `rcs_factors`. That choice makes the classical energy the exact expectation
# value in a spin coherent state, and LSWT a small-tilt expansion of it. The 1/s
# corrections instead expand the operator the user supplied, dividing λ_k back
# out order by order. In mode :dipole_uncorrected the stored coefficients are
# used directly, which lifts the user's classical polynomial to the operator
# with the same Stevens coefficients. That lift is a convention, a classical
# polynomial determining an operator only up to the relative order 1/s computed
# here.
#
# Two consequences for the sub-leading words, both gated against exact
# diagonalization in the test suite. Mode :dipole leaves only the anomalous word
# b†b† of an onsite anisotropy: a boson coherent state is not a spin coherent
# state, their amplitudes on one deviation agreeing but those on two differing
# by √(1 - 1/2s). In mode :dipole_uncorrected every word of at most two bosons
# is also rescaled, by the leading deviation of λ_k from unity (see
# `rcs_factor_first_order`). In either mode the one-boson word is proportional
# to the classical energy gradient, so an onsite anisotropy sources
# `tadpole_correction` only away from the classical minimum.


# Monomials of the term of K bosons of the Hamiltonian in the dipole modes,
# about the structure whose local frames are `frames`, beyond LSWT if
# `beyond_lswt`. For K ≥ 3, which LSWT does not contain, that is the whole term,
# kept at leading order; for K ≤ 2 it is the sub-leading part only, the O(1/s)
# correction to LSWT. Without `beyond_lswt` the K = 2 term is the LSWT
# Hamiltonian itself, which `QuadraticTable` compiles. For `twist`, see
# `local_interactions`.
#
# Some cubic words survive only away from a classical energy minimum. The cubic
# part of a transverse component, multiplying the classical Sᶻ of its partner or
# the Zeeman field, sums over the bonds of a site to its transverse effective
# field, which vanishes at a minimum exactly as H₁ does. These are retained
# because tadpole relaxation moves the structure off that minimum.
function hamiltonian_monomials(sys::System, frames::Vector{Mat3}, ::Val{K}; beyond_lswt=false, twist=nothing) where K
    L = nsites(sys)
    rcs = sys.mode == :dipole

    # The leading and sub-leading part, at each number of bosons, of the local
    # operator `O` at site i. The number of bosons is kept static, as in
    # `append_products!`.
    counts = ntuple(p -> Val{p-1}(), Val{K+1}())
    function words(O, i)
        (lead, sub) = stevens_words(O, √(2sys.κs[i]), rcs)
        return (map(v -> local_word_monomials(lead, i, L, v), counts),
                map(v -> local_word_monomials(sub, i, L, v), counts))
    end
    leading_only(w) = (w[1], map(_ -> BosonMonomial{0}[], w[2]))
    axes = (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1))
    spins = [map(r -> leading_only(words(spin_component(r), i)), axes) for i in 1:L]
    quads = [ntuple(q -> words(quadrupole(q), i), 5) for i in 1:L]

    # The kept part of the term of K bosons of the product of factors (A, B), on
    # sites whose cells are offset by n, with coefficient c
    terms = BosonMonomial{K}[]
    keep_lead = !beyond_lswt || K >= 3
    keep_sub = beyond_lswt && K <= 2
    function add_product!(c, (Al, As), (Bl, Bs), n)
        keep_lead && append_products!(terms, c, Al, Bl, n)
        if keep_sub
            append_products!(terms, c, As, Bl, n)
            append_products!(terms, c, Al, Bs, n)
        end
    end
    function add_local!((Al, As))
        keep_lead && append!(terms, Al[K+1])
        keep_sub && append!(terms, As[K+1])
    end

    for i in 1:L
        int = local_interactions(sys, frames, i; twist)
        B = frames[i]' * (sys.gs[i]' * sys.extfield[i])
        add_local!(leading_only(words(spin_component(B), i)))
        iszero(int.onsite) || add_local!(words(int.onsite, i))

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break
            (j, n) = bond_target(sys, i, bond)
            # Component a of 𝐒_i times the spin component (J 𝐒_j)_a
            J = coupling.bilin
            if !iszero(J)
                for a in 1:3
                    r = (J::Mat3)[a, :]
                    iszero(r) || add_product!(1, spins[i][a], leading_only(words(spin_component(r), j)), n)
                end
            end
            # Quadrupole p of site i times quadrupole q of site j
            Kq = coupling.biquad
            if !iszero(Kq)
                for p in 1:5, q in 1:5
                    iszero(Kq[p, q]) || add_product!(Kq[p, q], quads[i][p], quads[j][q], n)
                end
            end
        end
    end

    # Merged because the cubic list is contracted once per loop wavevector, the
    # innermost loop of the whole calculation, and because the quartic list is
    # Wick decoupled term by term
    return merge_monomials(terms)
end

# Monomials of the sub-leading part of the term of K ≤ 2 bosons, the O(1/s)
# correction to what LSWT treats at leading order. K = 1 enters through
# `tadpole_correction`, K = 2 through `anisotropy_correction`, and K = 0 is
# the correction to the classical energy. Empty in mode :SUN; see
# ExpansionSUN.jl.
function lswt_correction_monomials(swt::SpinWaveTheory, k::Val{K}) where K
    @assert K <= 2
    return hamiltonian_monomials(swt.sys, swt.frames, k; beyond_lswt=true)
end

"""
    anisotropy_correction(swt::SpinWaveTheory)

Correction of relative order ``1/s`` arising because linear spin wave theory
expands a classical energy function, whereas an onsite anisotropy is a quantum
operator. Returns `(; terms2, δE)`, where `terms2` is a correction to the
quadratic Hamiltonian that can be passed to [`static_self_energy`](@ref) and `δE`
is a correction to the classical energy per site. Both vanish in the absence of
anisotropy.

In mode `:dipole` the renormalization of Stevens coefficients derived in
[arXiv:2304.03874](https://arxiv.org/abs/2304.03874) already makes the classical
energy and the diagonal part of the quadratic Hamiltonian exact. There `δE`
vanishes and `terms2` is purely anomalous, i.e. of the form ``b†b†`` and ``bb``.

Both parts vanish in mode `:SUN`, where an onsite coupling enters the boson
Hamiltonian exactly and LSWT is already expanding an operator rather than a
classical energy function.
"""
function anisotropy_correction(swt::SpinWaveTheory)
    check_corrections_supported(swt)
    δE = sum(t -> real(t.c), lswt_correction_monomials(swt, Val{0}()); init=0.0)
    return (; terms2 = lswt_correction_monomials(swt, Val{2}()),
              δE = δE / nsites(uncontracted_system(swt.sys)))
end

# Monomials of H₃ and H₄, cached in `swt` and not to be mutated
cubic_monomials(swt::SpinWaveTheory) = cached!(swt.monomials, :cubic) do
    hamiltonian_monomials(swt.sys, swt.frames, Val{3}(); beyond_lswt=true)
end
quartic_monomials(swt::SpinWaveTheory) = cached!(swt.monomials, :quartic) do
    hamiltonian_monomials(swt.sys, swt.frames, Val{4}(); beyond_lswt=true)
end

# True when the cubic vertex is numerically zero, so that everything it
# generates — the self-energy, the tadpole, and their interference with the
# direct pair amplitude — vanishes and need not be computed. The usual cause is
# a collinear structure in a dipole mode, where the monomials cancel on merging
# rather than being absent: a two-sublattice antiferromagnet has 26 monomials
# whose coefficients sum to O(1e-11), so the test must be on magnitude and not
# on `isempty`. The scale is the quadratic Hamiltonian, which makes it
# dimensionless; measured ratios are 9e-13 for that antiferromagnet against
# 0.046 once a field cants it and 0.083 for the triangular-lattice 120°
# structure, so the threshold sits in a ten-order gap.
#
# Collinearity buys nothing in mode :SUN, and the same easy-axis Néel state
# gives 0.0497 there: the cancellation is between the two transverse components
# of a dipole, and the flavor-changing words that promote a single-ion level
# have no such partner. So this is a test of the assembled vertex rather than of
# the structure, and returns what it finds in either mode.
function cubic_vertex_vanishes(swt::SpinWaveTheory, terms3)
    isempty(terms3) && return true
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    dynamical_matrix!(H, swt, zero(Vec3))
    return maximum(abs(t.c) for t in terms3) < 1e-8 * norm(H)
end
