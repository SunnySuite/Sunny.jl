# Holstein-Primakoff words of single-site Stevens operators in the dipole
# modes. The truncation rule they implement is explained in ExpansionDipole.jl.

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

# The exact words of the local operator Σ_kq c_kq 𝒪_k^q at a site of spin s,
# given as in `stevens_words`. Words of at most four bosons read only matrix
# elements between the lowest seven boson states, which are closed-form at any
# s, so no series is needed: for half-integer s these are the words of the
# (2s+1)×(2s+1) matrix itself, and otherwise their continuation in s. An
# observable never enters the Hamiltonian, so unlike a coupling it is not
# truncated; see ObservableWords.jl.
function exact_stevens_words(cs, s)
    iszero(s) && return zero(Words)
    off = [√complex(j * (2s + 1 - j)) / 2 for j in 1:6]
    J = (diagm(1 => off, -1 => off), diagm(1 => -im*off, -1 => im*off), diagm(ComplexF64[s - n for n in 0:6]))
    𝒪 = sum(cs; init=zeros(ComplexF64, 7, 7)) do (k, c)
        sum(((ckq, O),) -> ckq * O, zip(c, stevens_abstract_polynomials(; J, k)); init=zeros(ComplexF64, 7, 7))
    end
    g(d, m) = 𝒪[m+d+1, m+1] / √prod(m+1:m+d; init=1)
    return Words(map(WORD_KEYS) do (d, j)
        sum(m -> (-1)^(j-m) * binomial(j, m) * g(d, m), 0:j) / factorial(j)
    end)
end

# Stevens coefficients of the spin component 𝐫 ⋅ 𝐒, since (𝒪_1^1, 𝒪_1^0, 𝒪_1^-1) =
# (Sˣ, Sᶻ, Sʸ), and the leading words of that component. Its sub-leading words of
# at most two bosons vanish: Sᶻ = s - b†b is exact, as are the one-boson words
# of S±.
spin_component_coefs(r) = ((1, (r[1], r[3], r[2])),)
spin_component_words(r, x) = first(stevens_words(spin_component_coefs(r), x, false))

# The quadrupole 𝒪_2^q as Stevens coefficients, as in `stevens_words`
quadrupole_coefs(q) = ((2, ntuple(p -> p == q ? 1.0 : 0.0, 5)),)
