# Holstein-Primakoff expansion of the Hamiltonian in the dipole modes, as boson
# monomials beyond the quadratic (LSWT) order. Overall conventions are collected
# in Corrections.jl and the monomial representation in Monomials.jl.

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

# Truncated Laurent polynomials in x = √(2s), which is the variable that grades the
# Holstein-Primakoff expansion: a spin operator of degree k contributes to the word
# of n bosons at leading order x^{2k-n}, with corrections descending in powers of
# x^{-2} = 1/2s. Coefficients outside `LAURENT_RANGE` are discarded, which is
# harmless as long as the window extends well below every power read back out.
const LAURENT_RANGE = -14:12
const NLAUR = length(LAURENT_RANGE)
const LAURENT_P0 = 1 - first(LAURENT_RANGE)  # index of x⁰

struct Laurent <: Number
    c::SVector{NLAUR, ComplexF64}
end

# Coefficient of xᵖ, zero outside the retained window, and the monomial z xᵖ.
coefficient(a::Laurent, p::Int) = p in LAURENT_RANGE ? a.c[p + LAURENT_P0] : 0.0im
monomial(p::Int, z=1) = Laurent(setindex(zero(SVector{NLAUR, ComplexF64}), z, p + LAURENT_P0))

Laurent(z::Number) = monomial(0, z)
Base.zero(::Type{Laurent}) = Laurent(zero(SVector{NLAUR, ComplexF64}))
Base.one(::Type{Laurent}) = monomial(0)
Base.promote_rule(::Type{Laurent}, ::Type{<:Union{Real, Complex}}) = Laurent
Base.:+(a::Laurent, b::Laurent) = Laurent(a.c + b.c)
Base.:-(a::Laurent, b::Laurent) = Laurent(a.c - b.c)
Base.:-(a::Laurent) = Laurent(-a.c)
Base.:/(a::Laurent, z::Number) = Laurent(a.c / z)
Base.:/(a::Laurent, b::Laurent) = a * inv(b)

# Powers add, so the coefficient of xᵖ is the convolution Σ_q a_q b_{p-q}, summed
# over the pairs of powers that both lie inside the retained window. Indices are
# used directly rather than powers because this is the innermost loop of the
# anisotropy expansion.
function Base.:*(a::Laurent, b::Laurent)
    function f(i)
        acc = 0.0im
        for t in max(1, i+LAURENT_P0-NLAUR):min(i+LAURENT_P0-1, NLAUR)
            acc += a.c[t] * b.c[i+LAURENT_P0-t]
        end
        return acc
    end
    return Laurent(SVector(ntuple(f, Val{NLAUR}())))
end

# Reciprocal. Factoring out the leading monomial leaves 1 + u with u carrying only
# negative powers, so the geometric series Σ (-u)ᵐ for its inverse falls out of the
# retained window after finitely many terms. Treating the highest retained power as
# the leading one is what makes the result a descending series in 1/s, as wanted
# here, rather than an ascending one.
function Base.inv(a::Laurent)
    p0 = something(findlast(!iszero, a.c), 0) - LAURENT_P0
    -p0 in LAURENT_RANGE || error("Cannot invert a Laurent series of leading power $p0")
    lead = monomial(-p0, inv(coefficient(a, p0)))
    u = lead * a - one(Laurent)
    r = one(Laurent)
    for _ in 1:NLAUR
        r = one(Laurent) - u * r
    end
    return lead * r
end

# Spin matrices in the basis of boson number n = s - m_z, running over 0:nmax.
# Mirrors `spin_matrices_of_dim` with s kept symbolic, the transverse elements
# becoming infinite series through √(j(2s+1-j)) = x √j √(1 + (1-j)/x²).
function graded_spin_matrices(nmax::Int)
    off = map(1:nmax) do j
        w = zero(Laurent)
        b = √j / 2
        for m in 0:div(1 - first(LAURENT_RANGE), 2)
            w += monomial(1 - 2m, b)
            b *= (1/2 - m) * (1 - j) / (m + 1)
        end
        return w
    end
    Jx = diagm(1 => off, -1 => off)
    Jy = diagm(1 => -im*off, -1 => im*off)
    Jz = diagm([monomial(2, 1/2) - n for n in 0:nmax])
    return (Jx, Jy, Jz)
end

# Graded spin matrices, and the reciprocals of the `rcs_factors`
# renormalizations. Only matrix elements between boson numbers ≤ 4 are ever
# read, and a Stevens operator of order k ≤ 6 shifts the boson number by at most
# k while its diagonal factors reach one further, so seven states suffice.
const GRADED_SPINS = graded_spin_matrices(6)
const GRADED_INV_RCS = map(inv, rcs_factors(monomial(2, 1/2)))

# Holstein-Primakoff words of a single-site operator `𝒪`, a graded matrix over
# boson numbers, as coefficients of the normal-ordered words (b†)^{d+j} bʲ keyed
# by (d, j), for d + 2j ≤ 4. The word of n = d + 2j bosons is read off at the
# single power x^p(n) of x = √(2s), which is how the truncation rule at the top
# of this file is imposed.
#
# A single-site operator is block diagonal in the boson number up to a shift,
#
#     𝒪 = g_0(n̂) + Σ_{d > 0} [(b†)^d g_d(n̂) + h.c.],   n̂ = b†b,
#
# and expanding each g_d in the Newton basis n̂(n̂-1)…(n̂-j+1) = (b†)ʲ bʲ gives the
# words. The coefficients follow from finite differences of the matrix elements
# of 𝒪, which are assembled over `Laurent` so that each word is obtained as a
# series in 1/s rather than as a number. Sunny lists the eigenvalues of Sᶻ in
# descending order, so the boson number of basis state m is m-1, and
# ⟨m+d|(b†)^d|m⟩ = √((m+d)!/m!) relates 𝒪 to g_d.
function local_words!(words, 𝒪, x, p)
    g(d, m) = 𝒪[m+d+1, m+1] / √prod(m+1:m+d; init=1)
    for d in 0:4, j in 0:div(4-d, 2)
        W = sum(m -> (-1)^(j-m) * binomial(j, m) * g(d, m), 0:j) / factorial(j)
        words[(d, j)] = get(words, (d, j), 0.0im) + coefficient(W, p(d + 2j)) * x^p(d + 2j)
    end
    return words
end

# Monomials of exactly K bosons at site i among `words`: each word (b†)^{d+j} bʲ
# and, unless it is self-adjoint, its adjoint (b†)ʲ b^{d+j}.
function local_word_monomials(words, i, L, ::Val{K}) where K
    ns = ntuple(_ -> zero(Vec3), Val{K}())
    terms = BosonMonomial{K}[]
    for d in K:-2:0
        j = div(K - d, 2)
        c = word(words, (d, j))
        iszero(c) && continue
        push!(terms, BosonMonomial(c, ntuple(t -> t <= d+j ? L+i : i, Val{K}()), ns))
        iszero(d) && continue
        push!(terms, BosonMonomial(conj(c), ntuple(t -> t <= j ? L+i : i, Val{K}()), ns))
    end
    return terms
end

# The normal-ordered words (b†)^{d+j} bʲ, keyed by (d, j), of at most four
# bosons, and their numbers of bosons. A set of words is a vector over these.
const WORD_KEYS = Tuple((d, j) for d in 0:4 for j in 0:div(4-d, 2))
const WORD_NB = map(((d, j),) -> d + 2j, WORD_KEYS)
const Words = SVector{length(WORD_KEYS), ComplexF64}

word(w::Words, key) = w[findfirst(==(key), WORD_KEYS)]

# The words of each Stevens operator 𝒪_k^q, for k ≤ 6 and ordered as in
# `stevens_matrices`, as a pair (lead, sub): the word of n bosons at its leading
# power x^{2k-n}, and at the next, x^{2k-n-2}, each at x = 1. A word scales with
# x as its power, so these depend on nothing, and are built once here.
const STEVENS_WORDS = map(0:6) do k
    map(stevens_abstract_polynomials(; J=GRADED_SPINS, k)) do 𝒪
        function words(p)
            w = local_words!(Dict{Tuple{Int, Int}, ComplexF64}(), 𝒪, 1.0, p)
            return Words(map(key -> w[key], WORD_KEYS))
        end
        (words(n -> 2k - n), words(n -> 2k - n - 2))
    end
end

# The leading correction to 1/λ_k, i.e. 1/λ_k = 1 + RCS_SHIFT[k]/2s + O(1/s²),
# for the renormalizations λ_k of `rcs_factors`.
const RCS_SHIFT = map(k -> real(coefficient(GRADED_INV_RCS[k], -2)), 0:6)

# The pair (lead, sub) of words, as in `STEVENS_WORDS`, of the local operator
# Σ_kq c_kq 𝒪_k^q at a site of spin s, given x = √2s and an iterable `cs` of
# pairs (k, c_k), with c_k the coefficients ordered as in `stevens_matrices`.
# Where `rcs` holds, the stored c_k carry the renormalization λ_k of mode
# :dipole, which the sub-leading words divide back out. A vacant site, x = 0, has
# no bosons and so no words.
function stevens_words(cs, x, rcs::Bool)
    iszero(x) && return (zero(Words), zero(Words))
    # The orders of `cs` may hold coefficient vectors of different types, so
    # they are visited by recursion rather than iteration, which stays static
    return mapreduce(.+, cs) do (k, c)
        lead = sub = zero(Words)
        powers = map(n -> x^(2k - n), Words(WORD_NB))
        for (q, ckq) in enumerate(c)
            iszero(ckq) && continue
            (wl, ws) = STEVENS_WORDS[k+1][q]
            lead += ckq * wl .* powers
            sub += ckq * (ws + rcs * RCS_SHIFT[k+1] * wl) .* powers / x^2
        end
        (lead, sub)
    end
end

# Stevens coefficients of the spin component 𝐫 ⋅ 𝐒, since (𝒪_1^1, 𝒪_1^0, 𝒪_1^-1) =
# (Sˣ, Sᶻ, Sʸ), and the leading words of that component. Its sub-leading words of
# at most two bosons vanish: Sᶻ = s - b†b is exact, as are the one-boson words
# of S±.
spin_component_coefs(r) = ((1, (r[1], r[3], r[2])),)
spin_component_words(r, x) = first(stevens_words(spin_component_coefs(r), x, false))

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
# diagonalization in the test suite. Mode :dipole leaves only the anomalous
# word b†b† of an onsite anisotropy: a boson coherent state is not a spin
# coherent state, their amplitudes on one deviation agreeing but those on two
# differing by √(1 - 1/2s). In mode :dipole_uncorrected every word of at most two
# bosons is also rescaled, by -binomial(k, 2)/2s, the leading deviation of λ_k
# from unity. In either mode the one-boson word is proportional to the
# classical energy gradient, so an onsite anisotropy sources
# `tadpole_correction` only away from the classical minimum.

# Each local operator of the Hamiltonian as Stevens coefficients, in the local
# frame: the onsite anisotropy of each site, which is time-reversal even, and
# the quadrupoles 𝒪_2^q, ordered as in `stevens_matrices`.
onsite_stevens_coefs(stvexp::StevensExpansion) =
    ((0, stvexp.c0), (2, stvexp.c2), (4, stvexp.c4), (6, stvexp.c6))
quadrupole_coefs(q) = ((2, ntuple(p -> p == q ? 1.0 : 0.0, 5)),)

# Monomials of the term of K bosons of the Hamiltonian in the dipole modes,
# beyond LSWT if `beyond_lswt`. For K ≥ 3, which LSWT does not contain, that is
# the whole term, kept at leading order; for K ≤ 2 it is the sub-leading part
# only, the O(1/s) correction to LSWT. Without `beyond_lswt` the K = 2 term is
# the LSWT Hamiltonian itself, which `QuadraticTable` compiles.
#
# Some cubic words survive only away from a classical energy minimum. The cubic
# part of a transverse component, multiplying the classical Sᶻ of its partner or
# the Zeeman field, sums over the bonds of a site to its transverse effective
# field, which vanishes at a minimum exactly as H₁ does. These are retained
# because tadpole relaxation moves the structure off that minimum.
function monomials_dipole(sys::System, data::SWTDataDipole, ::Val{K}; beyond_lswt) where K
    (; local_rotations, stevens_coefs, sqrtS, interactions) = data
    L = nsites(sys)
    rcs = sys.mode == :dipole

    # The leading and sub-leading part, at each number of bosons, of the local
    # operator with Stevens coefficients `cs` at site i. The number of bosons is
    # kept static, as in `append_products!`.
    counts = ntuple(p -> Val{p-1}(), Val{K+1}())
    function words(cs, i)
        (lead, sub) = stevens_words(cs, √2 * sqrtS[i], rcs)
        return (map(v -> local_word_monomials(lead, i, L, v), counts),
                map(v -> local_word_monomials(sub, i, L, v), counts))
    end
    leading_only(w) = (w[1], map(_ -> BosonMonomial{0}[], w[2]))
    axes = (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1))
    spins = [map(r -> leading_only(words(spin_component_coefs(r), i)), axes) for i in 1:L]
    quads = [ntuple(q -> words(quadrupole_coefs(q), i), 5) for i in 1:L]

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

    for (i, int) in enumerate(interactions)
        B = local_rotations[i]' * (sys.gs[i]' * sys.extfield[i])
        add_local!(leading_only(words(spin_component_coefs(B), i)))
        iszero(int.onsite) || add_local!(words(onsite_stevens_coefs(stevens_coefs[i]), i))

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break
            (j, n) = bond_target(sys, i, bond)
            # Component a of 𝐒_i times the spin component (J 𝐒_j)_a
            J = coupling.bilin
            if !iszero(J)
                for a in 1:3
                    r = (J::Mat3)[a, :]
                    iszero(r) || add_product!(1, spins[i][a], leading_only(words(spin_component_coefs(r), j)), n)
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
    # Wick decoupled term by term. The LSWT Hamiltonian is merged instead by
    # `QuadraticTable`, which `SpinWaveTheory` builds on every construction.
    return beyond_lswt ? merge_monomials(terms) : terms
end

monomials_dipole(swt::SpinWaveTheory, k::Val; beyond_lswt=true) =
    monomials_dipole(swt.sys, swt.data::SWTDataDipole, k; beyond_lswt)

# Monomials of the sub-leading part of the term of K ≤ 2 bosons, the O(1/s)
# correction to what LSWT treats at leading order. K = 1 enters through
# `tadpole_correction`, K = 2 through `anisotropy_correction`, and K = 0 is
# the correction to the classical energy. Empty in mode :SUN, where an onsite
# coupling is not a Stevens polynomial to be re-expanded but a matrix that
# `monomials_SUN` promotes exactly, leaving no remainder at any boson number.
# See ExpansionSUN.jl.
function lswt_correction_monomials(swt::SpinWaveTheory, k::Val{K}) where K
    @assert K <= 2
    swt.sys.mode == :SUN && return BosonMonomial{K}[]
    return monomials_dipole(swt, k)
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

# Monomials of H₂, the LSWT Hamiltonian, which `QuadraticTable` compiles for
# `dynamical_matrix!`, and of H₃ and H₄, cached in `swt` and not to be mutated.
# Mode :SUN expands in 1/M rather than 1/s and is handled by ExpansionSUN.jl,
# which shares everything below the monomial representation itself.
quadratic_monomials(sys::System, data) =
    sys.mode == :SUN ? monomials_SUN(sys, data, Val{2}()) : monomials_dipole(sys, data, Val{2}(); beyond_lswt=false)
cubic_monomials(swt::SpinWaveTheory) = cached!(swt.monomials, :cubic) do
    swt.sys.mode == :SUN ? monomials_SUN(swt, Val{3}()) : monomials_dipole(swt, Val{3}())
end
quartic_monomials(swt::SpinWaveTheory) = cached!(swt.monomials, :quartic) do
    swt.sys.mode == :SUN ? monomials_SUN(swt, Val{4}()) : monomials_dipole(swt, Val{4}())
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
