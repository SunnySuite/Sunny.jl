# Vertices of the Holstein-Primakoff expansion in mode :SUN, which plays the
# role that Vertices.jl plays for the dipole modes. Overall conventions are
# collected in Corrections.jl and the monomial representation in Vertices.jl.
#
# The expansion parameter is not the spin magnitude but the number of boxes M of
# the symmetric SU(N) representation, Sunny's :SUN mode being M = 1. Writing the
# condensate flavor as N and n̂ = Σ_{m<N} b†_m b_m for the number of bosons that
# have left it, a local operator A — an N×N matrix in the local frame that
# `swt_data!` rotates into — is promoted to the M-box representation as
#
#     Â = Σ_{mn<N} A[m,n] b†_m b_n + Σ_{m<N} (A[m,N] b†_m √(M-n̂) + A[N,m] √(M-n̂) b_m)
#         + A[N,N] (M - n̂).
#
# This is exact, being nothing but the Schwinger representation of the
# generators with the condensate flavor eliminated by Σ_α b†_α b_α = M.
# Expanding the square root grades the words by boson number, each word carrying
# a *single* power M^{1-nb/2}:
#
#     nb = 0:  A[N,N]
#     nb = 1:  A[m,N] b†_m + h.c.
#     nb = 2:  (A[m,n] - δ_mn A[N,N]) b†_m b_n
#     nb = 3:  -(1/2) A[m,N] b†_m n̂ + h.c.
#
# Three things follow, and together they make this file far shorter than
# Vertices.jl. There is no series at fixed word length, so none of the `Laurent`
# apparatus that the Stevens expansion of dipole mode requires is needed. There
# is no four-boson word at all, the square root contributing only odd ones. And
# the words of nb ≤ 2 are exactly what LSWT already holds —
# `swt_hamiltonian_SUN!` writes the nb = 2 word above verbatim, and the
# classical energy is A[N,N] — so unlike dipole mode there is no sub-leading
# remainder at low boson number, which is why `anisotropy_monomials` has nothing
# to return here.
#
# THE TRUNCATION RULE, in the form the discussion in Vertices.jl takes here:
# each n-boson sector is kept at its leading order in 1/M, which is its single
# word above. The first omission is the five-boson word -(1/8)A[m,N] b†_m n̂²,
# one full order below the three-boson one. Like the three-boson word it is
# proportional to A[m,N], whose sum over interactions is the gradient of the
# classical energy and so vanishes at a classical minimum; the omission
# therefore costs nothing there, and only O(1/M²) once tadpole relaxation has
# moved the structure off one. The same statement makes the three-boson vertex
# *exact* at a classical minimum, verified against exact diagonalization in the
# M-box representation for M = 4…7.
#
# That verification is sensitive to which minimum. The M-box classical energy is
# M Σᵢ onsiteᵢ[N,N] + M² Σ A[N,N] B[N,N], the two terms carrying different
# powers, so the stationary reference state depends on M and only at M = 1 — the
# one case the code is used in — is it the state Sunny's own minimizer finds.
# Evaluated instead at the M = 1 minimum, the cubic residual is nonzero and
# falls as 1/M, which is precisely the omitted five-boson word.
#
# Two consequences for the physics are worth flagging. Biquadratic exchange
# needs no special treatment, `swt_data!` having already decomposed every
# coupling into the (A, B) tensor pairs of `coupling.general.data` that the
# assembly below consumes. And the cubic vertex does *not* vanish for a
# collinear structure, as it does in dipole mode: flavor-changing words let a
# single-ion level decay into two magnons of different flavors, a channel with
# no dipole-mode counterpart.

# Words containing exactly `nb` bosons in the expansion of the local operator `A`
# at site `i`, each carrying the cell offset `off`. Sunny's Nambu labeling is used,
# so a label a ≤ L annihilates boson a and L+a creates it, and bosons are laid out
# as (flavor, atom) with flavor fastest.
#
# `M` is the number of boxes, which is 1 for Sunny's own :SUN mode. It is a
# parameter only so that the expansion can be tested where it is asymptotic rather
# than where it is being used, the M = 1 case having no room for a site to hold two
# bosons. Each word carries the single power M^(1-nb/2) that expanding √(M-n̂)
# gives it, so a pair coupling, being a product of two such factors, carries
# M^(2-K/2) at every split of its K bosons.
function local_words(A, i, off, ::Val{nb}, Nf, L, M=1) where nb
    N = Nf + 1
    terms = BosonMonomial{nb}[]
    ns = ntuple(_ -> off, Val{nb}())
    scale = M^((2 - nb)/2)
    add!(c, as) = iszero(c) || push!(terms, BosonMonomial(ComplexF64(scale * c), as, ns))
    b(m) = m + (i - 1) * Nf

    if nb == 0
        add!(A[N, N], ())
    elseif nb == 1
        for m in 1:Nf
            add!(A[m, N], (L + b(m),))
            add!(A[N, m], (b(m),))
        end
    elseif nb == 2
        # The δ_mn A[N,N] is the depletion of the condensate, Sunny's own
        # `swt_hamiltonian_SUN!` writing exactly this coefficient.
        for m in 1:Nf, n in 1:Nf
            add!(A[m, n] - (m == n) * A[N, N], (L + b(m), b(n)))
        end
    elseif nb == 3
        # -(1/2) A[m,N] b†_m n̂ and its adjoint. Ordered so that the pair is
        # manifestly adjoint word by word, which is what makes the assembled
        # Hamiltonian Hermitian exactly rather than to within round-off.
        for m in 1:Nf, k in 1:Nf
            add!(-A[m, N] / 2, (L + b(m), L + b(k), b(k)))
            add!(-A[N, m] / 2, (L + b(k), b(k), b(m)))
        end
    end

    return terms
end

# Monomials of the K-boson term of the expansion in mode :SUN. Each onsite
# coupling, which `swt_data!` has already absorbed the Zeeman term into,
# contributes its K-boson word directly. Each pair coupling A_i ⊗ B_j
# contributes every split of the K bosons between its two sites; all K+1 of them
# carry the same power M^{2-K/2}, so all are kept. See `local_words` for `M`,
# which is 1 in every use outside the tests.
function sun_monomials(swt::SpinWaveTheory, ::Val{K}, M=1) where K
    (; sys) = swt
    @assert allequal(sys.Ns)
    Nf = nflavors(swt)
    L = nbands(swt)
    o = zero(Vec3)

    terms = BosonMonomial{K}[]
    for (i, int) in enumerate(sys.interactions_union)
        append!(terms, local_words(int.onsite, i, o, Val{K}(), Nf, L, M))

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break
            @assert i == bond.i
            n = Vec3(bond.n)

            for (A, B) in coupling.general.data, p in 0:K
                for ta in local_words(A, i, o, Val{p}(), Nf, L, M),
                    tb in local_words(B, bond.j, n, Val{K-p}(), Nf, L, M)
                    push!(terms, BosonMonomial(ta.c * tb.c, (ta.as..., tb.as...),
                                               (ta.ns..., tb.ns...)))
                end
            end
        end
    end

    # Merged because the cubic list is contracted once per loop wavevector, the
    # innermost loop of the whole calculation, and because the quartic list is
    # Wick decoupled term by term. Both are several times longer here than their
    # dipole-mode counterparts.
    return merge_monomials(terms)
end
