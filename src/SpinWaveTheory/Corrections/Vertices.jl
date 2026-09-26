# Vertices of the Holstein-Primakoff expansion beyond the quadratic (LSWT)
# order. Overall conventions are collected in Corrections.jl.
#
# The expansion is most transparent in real space. Each monomial
#
#     c Σ_𝐫 ∏ₛ b^{σₛ}_{iₛ, 𝐫+𝐧ₛ}
#
# is a product of K boson operators, summed over the N magnetic cells 𝐫.
# Operator s acts on sublattice iₛ of the cell displaced by the integer offset
# 𝐧ₛ. Following Sunny's Nambu packing, the operator is labeled by an index a ∈
# 1:2L that selects b_i when a = i and b†_i when a = L+i. Operators within a
# monomial are stored in the order they are to be multiplied, which matters only
# when two of them act on the same site of the same cell.
struct BosonMonomial{K}
    c::ComplexF64
    as::NTuple{K, Int}
    ns::NTuple{K, Vec3}
end

# Sums the coefficients of monomials that describe the same operator product,
# which typically shrinks a decoupled term list several-fold. Worth doing
# because the resulting list is contracted once per point of a momentum-space
# integration. Offsets are keyed by integers, which are exactly comparable.
#
# Returned sorted by offset tuple, which lets `vertex!` compute the phases of
# one tuple once for the whole run of terms sharing it.
function merge_monomials(terms::Vector{BosonMonomial{K}}) where K
    ret = Dict{Tuple{NTuple{K, Int}, NTuple{K, NTuple{3, Int}}}, BosonMonomial{K}}()
    for term in terms
        key = (term.as, map(n -> round.(Int, Tuple(n)), term.ns))
        prev = get(ret, key, nothing)
        ret[key] = isnothing(prev) ? term : BosonMonomial(prev.c + term.c, term.as, term.ns)
    end
    return sort!(collect(values(ret)); by = t -> map(n -> round.(Int, Tuple(n)), t.ns))
end

# All permutations of (1, …, K), used to symmetrize a vertex over its slots.
function slot_permutations(K::Int)
    K == 1 && return [(1,)]
    return [(p[1:i-1]..., K, p[i:end]...) for p in slot_permutations(K-1) for i in 1:K]
end

# Cached for the K of interest, because `vertex!` needs them in the innermost loop
# of a momentum-space integration, where rebuilding the list costs more than the
# rest of the call.
const SLOT_PERMUTATIONS = ntuple(slot_permutations, 4)

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
# ferromagnet; see `anisotropy_words`.
#
# Writing a bilinear coupling in the raising/lowering basis,
#
#     Σ_{ab∈{x,y}} J_{ab} Sᵃ_i Sᵇ_j = (1/4) Σ_{pq∈{+,-}} c_{pq} Sᵖ_i S^q_j,
#     c_{pq} = J_{xx} - i η_q J_{xy} - i η_p J_{yx} - η_p η_q J_{yy},
#
# with η₊ = 1 and η₋ = -1, organizes the quartic term below. As a check on the
# convention, its leading part reproduces Sunny's quadratic coefficients:
# c₋₊ σ_i σ_j / 4 is `Q` and c₊₊ σ_i σ_j / 4 is `P`, where σ_i = √(2s_i).
function pm_coefficient(J, ηp, ηq)
    return J[1, 1] - im*ηq*J[1, 2] - im*ηp*J[2, 1] - ηp*ηq*J[2, 2]
end

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

# Graded Stevens operators 𝒪_k^q for even k, ordered as in `stevens_matrices`, and
# the reciprocals of the `rcs_factors` renormalizations. These depend on nothing,
# and building them dominates the cost of `anisotropy_words`, so they are built once
# here. Only matrix elements between boson numbers ≤ 4 are ever read, and a Stevens
# operator of order k ≤ 6 shifts the boson number by at most k while its diagonal
# factors reach one further, so seven states suffice.
const GRADED_STEVENS = let J = graded_spin_matrices(6)
    map(k -> stevens_abstract_polynomials(; J, k), (0, 2, 4, 6))
end
const GRADED_INV_RCS = map(inv, rcs_factors(monomial(2, 1/2)))

# Holstein-Primakoff expansion of the onsite anisotropy of site i, as coefficients
# of the normal-ordered words (b†)^{d+j} bʲ, keyed by (d, j) and truncated to
# d + 2j ≤ 4. A single-site operator is block diagonal in the boson number up to a
# shift,
#
#     𝒪 = g_0(n̂) + Σ_{d > 0} [(b†)^d g_d(n̂) + h.c.],   n̂ = b†b,
#
# and expanding each g_d in the Newton basis n̂(n̂-1)…(n̂-j+1) = (b†)ʲ bʲ gives those
# words. The coefficients follow from finite differences of the matrix elements of
# 𝒪, which are assembled over `Laurent` so that each word is obtained as a series
# in 1/s rather than as a number. Note that the transverse components generate
# words of every odd length, since S⁺ = σ√(1 - n̂/2s) b.
#
# Returned is only what LSWT does not already contain, and only to the order that
# the truncation rule at the top of this file prescribes. A word of n bosons has
# leading power x^{2k-n} for each Stevens order k, and
#
#   * for n ≤ 2 that leading power is exactly what LSWT keeps — the classical
#     energy, its gradient, and the coefficients A1 and A2 of
#     `swt_hamiltonian_dipole!` — so the sub-leading power x^{2k-n-2} is returned;
#   * for n = 3, 4 nothing is already present, so the leading power is returned.
#
# Reading off a fixed power per word keeps this construction consistent with the truncated
# vertices of `cubic_monomials` and `quartic_monomials`, so no separate subtraction of the
# LSWT terms is needed.
#
# Sunny stores an anisotropy as coefficients of the Stevens operators 𝒪_k^q. In mode :dipole
# those carry the renormalization `rcs_factors`, chosen so that the classical energy
# reproduces the exact expectation value in a spin coherent state; dividing it out recovers
# the operator the user supplied, and because the factor is itself a series in 1/s it must
# be divided out order by order. Mode :dipole_uncorrected instead uses the stored
# coefficients directly, lifting the user's classical polynomial to the operator with the
# same Stevens coefficients. That lift is a convention, a classical polynomial determining
# an operator only up to the relative order 1/s computed here.
#
# For n ≤ 2 the result is therefore two scalars times the leading word L that LSWT already
# holds, gated against `swt_hamiltonian_dipole!` in the test suite:
#
#   * The lift rescales every such word by -binomial(k, 2)/2s, the leading deviation of
#     `rcs_factors` from unity, in mode :dipole_uncorrected only. Mode :dipole is the more
#     accurate choice precisely because λ_k removes it, leaving the classical energy, its
#     gradient, and A1 with no correction at all.
#   * A boson coherent state is not a spin coherent state: their amplitudes on one spin
#     deviation agree, but on two they differ by √(1 - 1/2s). LSWT reads the anomalous A2
#     off the classical energy and so gets it too small by that factor, which both modes
#     correct by a further +L/4s.
#
# In mode :dipole that leaves A2 as the only correction below three bosons. The one-boson
# word is proportional to the classical energy gradient in both modes, so an onsite
# anisotropy sources `tadpole_correction` only away from the classical minimum, and never in
# mode :dipole.
function anisotropy_words(swt::SpinWaveTheory, i::Int)
    (; sys, data) = swt
    stvexp = data.stevens_coefs[i]
    x = √2 * data.sqrtS[i]

    words = Dict((d, j) => 0.0im for d in 0:4 for j in 0:div(4-d, 2))
    for (k, c) in zip((0, 2, 4, 6), (stvexp.c0, stvexp.c2, stvexp.c4, stvexp.c6))
        iszero(c) && continue
        𝒪 = c' * GRADED_STEVENS[div(k, 2) + 1]
        if sys.mode == :dipole
            𝒪 = 𝒪 * GRADED_INV_RCS[k]
        end
        # Sunny lists the eigenvalues of Sᶻ in descending order, so the boson
        # number of basis state m is m-1, and ⟨m+d|(b†)^d|m⟩ = √((m+d)!/m!)
        # relates 𝒪 to g_d.
        g(d, m) = 𝒪[m+d+1, m+1] / √prod(m+1:m+d; init=1)
        for (dj, _) in words
            (d, j) = dj
            W = sum(m -> (-1)^(j-m) * binomial(j, m) * g(d, m), 0:j) / factorial(j)
            p = 2k - (d + 2j) - (d + 2j <= 2 ? 2 : 0)
            words[dj] += coefficient(W, p) * x^p
        end
    end

    return words
end

# Monomials of the onsite anisotropies containing exactly K bosons, beyond what
# LSWT already accounts for. K = 3 and 4 complete the cubic and quartic
# vertices, while K = 1 and 2 correct terms that LSWT treats at leading order
# and enter through `tadpole_correction` and `anisotropy_correction`
# respectively.
#
# Empty in mode :SUN, where an onsite coupling is not a Stevens polynomial to be
# re-expanded but a matrix that `sun_monomials` promotes exactly, leaving no
# remainder at any boson number. See VerticesSUN.jl.
function anisotropy_monomials(swt::SpinWaveTheory, ::Val{K}) where K
    L = nbands(swt)
    terms = BosonMonomial{K}[]
    swt.sys.mode == :SUN && return terms
    ns = ntuple(_ -> zero(Vec3), Val{K}())

    for i in 1:L
        iszero(swt.sys.interactions_union[i].onsite) && continue
        words = anisotropy_words(swt, i)
        for d in K:-2:0
            j = div(K - d, 2)
            c = words[(d, j)]
            iszero(c) && continue
            # The word (b†)^{d+j} bʲ and, unless it is self-adjoint, its adjoint
            # (b†)ʲ b^{d+j}.
            push!(terms, BosonMonomial(c, ntuple(t -> t <= d+j ? L+i : i, Val{K}()), ns))
            iszero(d) && continue
            push!(terms, BosonMonomial(conj(c), ntuple(t -> t <= j ? L+i : i, Val{K}()), ns))
        end
    end

    return terms
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
    δE = sum(1:nbands(swt); init=0.0) do i
        (swt.sys.mode == :SUN || iszero(swt.sys.interactions_union[i].onsite)) ? 0.0 :
            real(anisotropy_words(swt, i)[(0, 0)])
    end
    return (; terms2 = anisotropy_monomials(swt, Val{2}()),
              δE = δE / nsites(uncontracted_system(swt.sys)))
end

# Monomials of H₃ and H₄, the three- and four-boson terms. Mode :SUN expands in
# 1/M rather than 1/s and is handled by VerticesSUN.jl, which shares everything
# below the monomial representation itself.
cubic_monomials(swt::SpinWaveTheory) =
    swt.sys.mode == :SUN ? sun_monomials(swt, Val{3}()) : cubic_monomials_dipole(swt)
quartic_monomials(swt::SpinWaveTheory) =
    swt.sys.mode == :SUN ? sun_monomials(swt, Val{4}()) : quartic_monomials_dipole(swt)

# Monomials of H₃, the three-boson term. A cubic term arises either as the
# product of a transverse component on one site of a bond with the fluctuating
# part of Sᶻ on the other, or as the cubic part of a single transverse component
# multiplying the classical Sᶻ of its partner. Both scale as s^(1/2), i.e. they
# are smaller than H₂ by s^(-1/2). Every such monomial changes the boson number by
# ±1, so the source vertex ⟨3 magnons|H₃|0⟩ appears only after the Bogoliubov
# transformation mixes b with b†. An onsite anisotropy contributes as well, and
# does produce ±3 monomials directly.
function cubic_monomials_dipole(swt::SpinWaveTheory)
    (; sys, data) = swt
    (; local_rotations, sqrtS) = data
    L = nbands(swt)

    terms = BosonMonomial{3}[]
    add!(c, as, ns) = iszero(c) || push!(terms, BosonMonomial(ComplexF64(c), as, ns))
    o = zero(Vec3)

    for (i, int) in enumerate(sys.interactions_union)
        # The Zeeman coupling is +𝐁⋅𝐒, as implied by the quadratic Hamiltonian.
        # Its transverse part survives at a classical minimum, where exchange
        # rather than the local frame cancels it, and enters through the cubic
        # part of S⁺ and S⁻.
        R = local_rotations[i]
        B = sys.gs[1, 1, 1, i]' * sys.extfield[1, 1, 1, i]
        add!(-dot(B, R[:, 1] - im*R[:, 2]) / (4√2 * sqrtS[i]), (L+i, i, i), (o, o, o))
        add!(-dot(B, R[:, 1] + im*R[:, 2]) / (4√2 * sqrtS[i]), (L+i, L+i, i), (o, o, o))

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break
            @assert i == bond.i
            iszero(coupling.bilin) && continue

            j = bond.j
            n = Vec3(bond.n)
            J = coupling.bilin::Mat3
            si = sqrtS[i]^2
            sj = sqrtS[j]^2

            # Transverse component on one site times the -b†b part of Sᶻ on the
            # other. Because the rotated exchange matrix is real, the amplitude
            # to create a transverse deviation is the conjugate of the amplitude
            # to destroy one.
            V31 = -√(si/2) * (J[1, 3] - im*J[2, 3])
            V32 = -√(sj/2) * (J[3, 1] - im*J[3, 2])
            add!(V31, (i, L+j, j), (o, n, n))
            add!(conj(V31), (L+i, L+j, j), (o, n, n))
            add!(V32, (L+i, i, j), (o, o, n))
            add!(conj(V32), (L+i, i, L+j), (o, o, n))

            # Cubic part of that same transverse component, multiplying the
            # classical Sᶻ = s of the other site. Relative to the terms above
            # this carries a factor s_other / 4s_self.
            #
            # Like the Zeeman terms, these are proportional to the transverse
            # effective field hˣ - i hʸ on the site in question, so summed over
            # all bonds they cancel whenever the structure is a classical energy
            # minimum, exactly as H₁ does. They are retained because tadpole
            # relaxation moves the structure off that minimum.
            add!(V31 * sj/4si, (L+i, i, i), (o, o, o))
            add!(conj(V31) * sj/4si, (L+i, L+i, i), (o, o, o))
            add!(V32 * si/4sj, (L+j, j, j), (n, n, n))
            add!(conj(V32) * si/4sj, (L+j, L+j, j), (n, n, n))
        end
    end

    # Merged because this list is contracted once per loop wavevector, which is the
    # innermost loop of the whole calculation.
    return merge_monomials([terms; anisotropy_monomials(swt, Val{3}())])
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

# Monomials of H₄, the four-boson term, which is smaller than H₂ by s^(-1).
# Since Sᶻ is exact and the transverse components have no four-boson part, only
# two families survive: the product of the longitudinal fluctuations on the two
# sites of a bond, and the product of a linear transverse component on one site
# with the cubic part of a transverse component on the other. A Zeeman coupling,
# being linear in 𝐒, contributes nothing here, but an onsite anisotropy does.
function quartic_monomials_dipole(swt::SpinWaveTheory)
    (; sys, data) = swt
    (; sqrtS) = data
    L = nbands(swt)

    terms = BosonMonomial{4}[]
    add!(c, as, ns) = iszero(c) || push!(terms, BosonMonomial(ComplexF64(c), as, ns))
    o = zero(Vec3)

    for (i, int) in enumerate(sys.interactions_union)
        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break
            @assert i == bond.i
            iszero(coupling.bilin) && continue

            j = bond.j
            n = Vec3(bond.n)
            J = coupling.bilin::Mat3

            # Longitudinal: J_zz (b†b)_i (b†b)_j
            add!(J[3, 3], (L+i, i, L+j, j), (o, o, n, n))

            # Transverse: one site supplies Sᵖ at linear order, σ_i b or σ_i b†,
            # and the other supplies the cubic part of S^q, -Θ^q/2σ_j with
            # Θ⁺ = b†bb and Θ⁻ = b†b†b. Both assignments of the two roles occur.
            for (ηp, p) in ((1, i), (-1, L+i)), (ηq, q) in ((1, j), (-1, L+j))
                V = -pm_coefficient(J, ηp, ηq) / 8
                Θi = ηp > 0 ? (L+i, i, i) : (L+i, L+i, i)
                Θj = ηq > 0 ? (L+j, j, j) : (L+j, L+j, j)
                add!(V * sqrtS[i]/sqrtS[j], (p, Θj...), (o, n, n, n))
                add!(V * sqrtS[j]/sqrtS[i], (Θi..., q), (o, o, o, n))
            end
        end
    end

    return [terms; anisotropy_monomials(swt, Val{4}())]
end

# Fourier transforming a monomial and summing over cells yields
#
#     Σ_𝐫 ∏ₛ b^{σₛ}_{iₛ, 𝐫+𝐧ₛ} = N^{1-K/2} Σ_{Σ𝐤=0} e^{2πi Σₛ 𝐤ₛ⋅𝐧ₛ} ∏ₛ x_{𝐤ₛ}[aₛ],
#
# where x_𝐤 = [b_𝐤; b†_{-𝐤}] is the Nambu vector, so that every component of
# x_𝐤 removes momentum 𝐤 and the momenta of a monomial sum to zero. This phase
# convention agrees with `swt_hamiltonian_dipole!`, whose b†_i b_j coefficient
# carries `cis(2π 𝐪⋅𝐧)` for a bond offset 𝐧. Writing x_𝐤 = T_𝐤 y_𝐤 with
# y_𝐤 = [α_𝐤; α†_{-𝐤}] then gives, for example,
#
#     H₃ = N^(-1/2) Σ_{𝐤₁+𝐤₂+𝐤₃=0} Σ_{n₁n₂n₃} U₃(𝐤; n) y_{𝐤₁}[n₁] y_{𝐤₂}[n₂] y_{𝐤₃}[n₃].
#
# Bands range over the full Nambu space 1:2L, so a single tensor holds every
# channel: a slot taken from 1:L annihilates a quasi-particle and a slot taken
# from L+1:2L creates one. In particular the decay vertex Γ₁ and the source
# vertex Γ₂ are both read off from U₃.
#
# Averaging over the K! assignments of a monomial's operators to slots makes the
# result symmetric under simultaneous permutation of the (momentum, band) pairs.
# This is exact up to commutators, which reorder operators only into terms of
# lower boson number, and so contribute to the linear term addressed by tadpole
# relaxation rather than to the vertex itself.
#
# The work is arranged as an accumulation followed by a change of basis, rather
# than as a sum of rank-one tensors. Each (monomial, permutation) contributes a
# single coefficient to the Nambu operator product it names, so the whole term
# list collapses into one (2L)^K tensor C at a cost independent of L;
# transforming its slots to the Bogoliubov basis is then K matrix products.
# Summing rank-one tensors instead costs K! per monomial times the (2L)^K size
# of the output, which for the cubic vertex of a three-band system is seventy
# times more arithmetic. `scratch` holds C and the intermediates; the innermost
# loop of a momentum-space integration should pass one in to be reused.
function vertex!(U::Array{ComplexF64, K}, terms::Vector{BosonMonomial{K}},
                 qs::NTuple{K, Vec3}, Ts::NTuple{K, Matrix{ComplexF64}},
                 scratch::Array{ComplexF64, K}=similar(U)) where K
    @assert all(x -> abs(x - round(x)) < 1e-12, sum(qs)) "Vertex momenta must sum to zero"
    # A change of basis leaves a vanishing tensor vanishing, so an empty term
    # list skips the matrix products entirely. Worth a branch because callers
    # that have discarded a negligible vertex still run the loop for its other
    # consumers.
    isempty(terms) && return fill!(U, 0)
    N = size(U, 1)
    # The cache is indexed by a value rather than a type, so the lookup must be
    # annotated for the loop below to be type stable. That loop is the innermost
    # one of a momentum-space integration.
    perms = SLOT_PERMUTATIONS[K]::Vector{NTuple{K, Int}}

    # Coefficient of the operator product ∏ₜ x_{𝐤ₜ}[aₜ]. Zero offsets are
    # common and carry no phase, so they are given none. Many terms share one
    # tuple of offsets — a triangular-lattice cubic vertex has 52 of them over
    # 13 distinct tuples — and `merge_monomials` groups those together, so the
    # K² phases need recomputing only when the tuple changes. That reuse is
    # worth having because this is the innermost loop of a momentum-space
    # integration, and the `cis` calls dominate it.
    fill!(scratch, 0)
    local ph
    for (i, term) in enumerate(terms)
        if i == 1 || term.ns != terms[i-1].ns
            ph = ntuple(Val{K}()) do t
                ntuple(s -> iszero(term.ns[s]) ? 1.0+0im : cis(2π * dot(qs[t], term.ns[s])), Val{K}())
            end
        end
        for p in perms
            c = term.c / length(perms)
            for t in 1:K
                c *= ph[t][p[t]]
            end
            scratch[CartesianIndex(ntuple(t -> term.as[p[t]], Val{K}()))] += c
        end
    end

    # Transform one slot at a time. With the slot being transformed leading, the
    # contraction is a matrix product, and writing its result last cycles the next
    # slot into the leading position, so no index permutation is ever needed.
    (A, B) = (scratch, U)
    for t in 1:K
        mul!(reshape(B, N^(K-1), N), transpose(reshape(A, N, N^(K-1))), Ts[t])
        (A, B) = (B, A)
    end
    A === U || copyto!(U, A)

    return U
end

# Bogoliubov transformations at each of the given wavevectors.
function bogoliubov_matrices(swt::SpinWaveTheory, qs::NTuple{K, Vec3}) where K
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    return ntuple(K) do t
        T = zeros(ComplexF64, 2L, 2L)
        dynamical_matrix!(H, swt, qs[t])
        bogoliubov!(T, H)
        return T
    end
end

# Vertex for the given momenta, which must sum to zero.
function vertex(swt::SpinWaveTheory, terms::Vector{BosonMonomial{K}}, qs::NTuple{K, Vec3}) where K
    L = nbands(swt)
    U = zeros(ComplexF64, ntuple(_ -> 2L, K))
    return vertex!(U, terms, qs, bogoliubov_matrices(swt, qs))
end
