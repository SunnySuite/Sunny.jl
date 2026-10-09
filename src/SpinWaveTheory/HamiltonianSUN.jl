
# Bilinear-boson coefficients for a pair coupling A ⊗ B on units (i, j), used
# by the Ewald term below. The other couplings are compiled from the same
# expansion, see ExpansionSUN.jl.
#
# In SU(N) spin-wave theory each unit carries N-1 Holstein-Primakoff bosons b_m.
# Once operators are rotated so the condensate is basis state N, a local
# operator O expands (to the order that matters here) into three pieces:
#
#     O  =  O[N,N]                                 # mean ⟨O⟩
#        +  Σ_m O[m,N] b†_m + O[N,m] b_m           # linear: create / destroy a boson
#        +  Σ_mn (O[m,n] - δ_mn O[N,N]) b†_m b_n   # quadratic: hop between flavors
#
# The pair term A⊗B is kept only to boson-bilinear order, which happens two
# ways: quadratic×mean (stays on one unit) and linear×linear (couples the two
# units). That yields exactly five flavor-(m,n) coefficients:
#
#     c1  (A[m,n] - δ_mn⟨A⟩)⟨B⟩   quadratic A on i, mean of B    → b†_{m,i} b_{n,i}
#     c2  ⟨A⟩(B[m,n] - δ_mn⟨B⟩)   mean of A, quadratic B on j    → b†_{m,j} b_{n,j}
#     c3  A[m,N] B[N,n]           raise i, lower j (i→j hop)     → b†_{m,i} b_{n,j}
#     c4  A[N,m] B[n,N]           lower i, raise j (j→i hop)     → b_{m,i} b†_{n,j}
#     c5  A[m,N] B[n,N]           raise both (anomalous pairing) → b†_{m,i} b†_{n,j}
#
# `accum_pair_dense!` files these into the Nambu blocks (H11 = b†b, H22 = bb†
# mirror, H12 = b†b† pairing, H21 = H12†).
@inline function swt_pair_coefs(A, B, m, n, N)
    c1 = (A[m,n] - δ(m,n)*A[N,N]) * B[N,N]     # quadratic A on i, ⟨B⟩
    c2 = A[N,N] * (B[m,n] - δ(m,n)*B[N,N])     # ⟨A⟩, quadratic B on j
    c3 = A[m,N] * B[N,n]                       # i→j hop
    c4 = A[N,m] * B[n,N]                       # j→i hop
    c5 = A[m,N] * B[n,N]                       # anomalous (pairing)
    return (c1, c2, c3, c4, c5)
end

# Scatter a pair coupling A ⊗ B on units (i,j) into the dense Nambu blocks. The
# on-unit blocks (c1,c2) carry no wrapping phase, so they take the coupling at
# q=0, `J0`; the inter-unit blocks (c3,c4,c5) take the coupling at wavevector q,
# `Jq`. For Ewald these are J0 = A(q=0) and Jq = A(q), the Fourier-space
# kernel.
@inline function accum_pair_dense!(H11, H12, H21, H22, A, B, i, j, J0, Jq, N)
    for m in 1:N-1, n in 1:N-1
        (c1, c2, c3, c4, c5) = swt_pair_coefs(A, B, m, n, N)
        H11[m,i,n,i] += c1*J0;             H22[n,i,m,i] += c1*J0
        H11[m,j,n,j] += c2*J0;             H22[n,j,m,j] += c2*J0
        H11[m,i,n,j] += c3*Jq;             H22[n,j,m,i] += c3*conj(Jq)
        H11[n,j,m,i] += c4*conj(Jq);       H22[m,i,n,j] += c4*Jq
        H12[m,i,n,j] += c5*Jq;             H12[n,j,m,i] += c5*conj(Jq)
        H21[n,j,m,i] += conj(c5)*conj(Jq); H21[m,i,n,j] += conj(c5)*Jq
    end
end

# Adds the long-range dipole-dipole interaction to the dynamical matrix `H` in
# SU(N) mode. Everything else is compiled into the quadratic table. For an
# entangled system the Ewald data is that of its uncontracted system.
function accum_ewald_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; spin_ops) = data
    N = sys.Ns[1]
    Na = natoms(sys.crystal)
    L = (N-1) * Na
    usys = uncontracted_system(sys)

    blockdims = (N-1, Na, N-1, Na)
    H11 = reshape(view(H, 1:L, 1:L), blockdims)
    H12 = reshape(view(H, 1:L, L+1:2L), blockdims)
    H21 = reshape(view(H, L+1:2L, 1:L), blockdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L), blockdims)
    (; gs, ewald) = usys
    (; cache, A) = ewald
    Nbare = natoms(usys.crystal)

    # Interaction matrix for wavevector (0,0,0), equivalently:
    # real.(ewald_interaction_tensor(cache, zero(Vec3)))
    A0 = reshape(A, Nbare, Nbare)
    # Interaction matrix for wavevector q
    Aq = ewald_interaction_tensor(cache, q_reshaped)
    Aq = reshape(Aq, Nbare, Nbare)

    # For an entangled system, `contract_ewald!` has already stripped
    # intra-unit couplings in A0 and shifted them to an onsite coupling. The
    # same stripping must be done for the freshly computed Aq.
    if is_entangled(sys)
        for (; ai, aj, Adir) in intra_unit_dipole_terms(usys.crystal, get_entanglement(sys).units, cache.μ0_μB²)
            Aq[ai, aj] -= Adir
        end
    end

    # The Ewald matrices and g-tensors are indexed by bare crystal-atom `a`
    for i in 1:Na, j in 1:Na
        for ai in atoms_in_unit(sys, i), aj in atoms_in_unit(sys, j)
            # An ordered pair of magnetic moments contribute (μₐ A μ_b)/2 to the
            # energy, where μ = - g S. A symmetric contribution will appear for
            # the bond reversal (a, b) → (b, a).
            Jq = gs[ai]' * Aq[ai, aj] * gs[aj] / 2
            J0 = gs[ai]' * A0[ai, aj] * gs[aj] / 2

            for α in 1:3, β in 1:3
                Ai = spin_ops[α, ai]
                Bj = spin_ops[β, aj]
                accum_pair_dense!(H11, H12, H21, H22, Ai, Bj, i, j, J0[α,β], Jq[α,β], N)
            end
        end
    end
end
