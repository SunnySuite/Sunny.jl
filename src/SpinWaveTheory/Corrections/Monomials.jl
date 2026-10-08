# Boson monomials, the representation in which every term of the
# Holstein-Primakoff expansion is held, and their contraction into Nambu
# matrices and vertex tensors. Overall conventions are collected in
# Corrections.jl.
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

# Nambu index of the adjoint operator, exchanging b_i and b†_i.
nambu_conj(a, L) = mod1(a + L, 2L)

# Accumulates quadratic monomials into the Nambu matrix H, in the convention
# H₂ = (1/2) x†_𝐪 H_𝐪 x_𝐪 + const of `swt_hamiltonian_dipole!`. Fourier
# transforming a monomial as in `vertex!` gives
#
#     c Σ_𝐫 O_{a₁}(𝐫+𝐧₁) O_{a₂}(𝐫+𝐧₂) = Σ_𝐪 c φ x_𝐪[a₁] x_𝐪[ā₂]†,
#
# where φ = exp(2πi 𝐪⋅(𝐧₁-𝐧₂)) and ā = mod1(a+L, 2L) labels the adjoint. Up to a
# constant from the commutator this is c φ x†[ā₂] x[a₁]. Since Σ_𝐪 x†_𝐪[α] x_𝐪[β]
# is invariant under (α, β, φ) → (β̄, ᾱ, conj(φ)), the weight can be split evenly
# between the two forms; that is what feeds a single b†b monomial into both the
# (1,1) and the (2,2) Nambu block. Applied to the monomials of H₂ itself, the rule
# reproduces every entry written by `swt_hamiltonian_dipole!`.
#
# The rule is covariant under taking adjoints, so a monomial list describing a
# Hermitian operator yields a Hermitian matrix. Two monomials can nonetheless
# describe the same operator while being written differently, in which case
# Hermiticity holds only to the accuracy with which their coefficients were
# computed; `hermitianpart!` projects out that residual.
function accum_quadratic!(H, terms::Vector{BosonMonomial{2}}, q_reshaped)
    L = div(size(H, 1), 2)
    for (; c, as, ns) in terms
        φ = cis(2π * dot(q_reshaped, ns[1] - ns[2]))
        H[nambu_conj(as[2], L), as[1]] += c * φ
        H[nambu_conj(as[1], L), as[2]] += c * conj(φ)
    end
    hermitianpart!(H)
end

# Site owning boson `a` of the Nambu labeling, bosons being laid out as (flavor,
# atom) with flavor fastest, `Nf` of them per site. The identity in dipole mode.
boson_site(a, L, Nf) = div(mod1(a, L) - 1, Nf) + 1

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
