
# Bilinear-boson coefficients for a pair coupling A ⊗ B on units (i, j)
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
# `Jq`. In the ordinary pair-coupling case J0 = J and Jq = J exp(2πi q·n), with
# J = 1 if the strength is already absorbed into A ⊗ B. The Ewald path reuses
# this same kernel with J0 = A(q=0) and Jq = A(q), the Fourier-space kernel.
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

# Matrix-free counterpart of `accum_pair_dense!`: accumulate y += (pair block) x
# over all wavevectors, keeping Holstein-Primakoff state in Nambu form
# X[q, flavor, unit, particle/hole]. Same weights as the dense builder, but the
# coupling at wavevector q is now an array `Jqs` over all q: in the ordinary case
# J0 = 1 and Jqs[q] is the periodic-wrapping phase; the entangled-bilinear path
# passes J0 = J and Jqs[q] = J·phase[q].
@inline function accum_pair_matvec!(Y, X, A, B, i, j, J0, Jqs, N)
    for m in 1:N-1, n in 1:N-1
        (c1, c2, c3, c4, c5) = swt_pair_coefs(A, B, m, n, N)
        @inbounds for q in axes(Y, 1)
            Jq = Jqs[q]
            Y[q, m, i, 1] += c1 * J0 * X[q, n, i, 1]
            Y[q, n, i, 2] += c1 * J0 * X[q, m, i, 2]

            Y[q, m, j, 1] += c2 * J0 * X[q, n, j, 1]
            Y[q, n, j, 2] += c2 * J0 * X[q, m, j, 2]

            Y[q, m, i, 1] += c3 * Jq * X[q, n, j, 1]
            Y[q, n, j, 2] += c3 * conj(Jq) * X[q, m, i, 2]

            Y[q, n, j, 1] += c4 * conj(Jq) * X[q, m, i, 1]
            Y[q, m, i, 2] += c4 * Jq * X[q, n, j, 2]

            Y[q, m, i, 1] += c5 * Jq * X[q, n, j, 2]
            Y[q, n, j, 1] += c5 * conj(Jq) * X[q, m, i, 2]
            Y[q, n, j, 2] += conj(c5) * conj(Jq) * X[q, m, i, 1]
            Y[q, m, i, 2] += conj(c5) * Jq * X[q, n, j, 1]
        end
    end
end

# Set the dynamical quadratic Hamiltonian matrix in SU(N) mode.
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped::Vec3)
    (; sys, data) = swt
    (; spin_ops) = data

    N = sys.Ns[1]
    Na = natoms(sys.crystal)
    L = (N-1) * Na

    # Clear the Hamiltonian
    @assert size(H) == (2L, 2L)
    H .= 0
    blockdims = (N-1, Na, N-1, Na)
    H11 = reshape(view(H, 1:L, 1:L), blockdims)
    H12 = reshape(view(H, 1:L, L+1:2L), blockdims)
    H21 = reshape(view(H, L+1:2L, 1:L), blockdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L), blockdims)

    for (i, int) in enumerate(sys.interactions_union)

        # Onsite coupling, including Zeeman. Note that op has already been
        # transformed according to the local frame of sublattice i.
        op = int.onsite
        for m in 1:N-1
            for n in 1:N-1
                c = op[m, n] - δ(m, n) * op[N, N]
                H11[m, i, n, i] += c
                H22[n, i, m, i] += c
            end
        end

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break

            @assert i == bond.i
            j = bond.j

            phase = cis(2π * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            # Set "general" pair interactions of the form Aᵢ⊗Bⱼ. Note that Aᵢ
            # and Bᵢ have already been transformed according to the local frames
            # of sublattice i and j, respectively.
            for (Ai, Bj) in coupling.general.data
                accum_pair_dense!(H11, H12, H21, H22, Ai, Bj, i, j, 1, phase, N)
            end
        end
    end

    # Long-range dipole-dipole interactions. In case of an entangled system,
    # this data would be stored in its uncontracted system.
    usys = uncontracted_system(sys)

    if !isnothing(usys.ewald)
        (; gs, ewald) = usys
        (; demag, μ0_μB², A) = ewald
        Nbare = natoms(usys.crystal)

        # Interaction matrix for wavevector (0,0,0). It could be recalculated as:
        # precompute_dipole_ewald(usys.crystal, (1,1,1), demag) * μ0_μB²
        A0 = reshape(A, Nbare, Nbare)
        # Interaction matrix for wavevector q
        Aq = precompute_dipole_ewald_at_wavevector(usys.crystal, (1,1,1), demag, q_reshaped) * μ0_μB²
        Aq = reshape(Aq, Nbare, Nbare)

        # For an entangled system, `contract_ewald!` has already stripped
        # intra-unit couplings in A0 and shifted them to an onsite coupling. The
        # same stripping must be done for the freshly computed Aq.
        if is_entangled(sys)
            for (; ai, aj, Adir) in intra_unit_dipole_terms(usys.crystal, get_entanglement(sys).units, μ0_μB²)
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

    # For an entangled system, inter-unit bilinear exchange is delegated to the
    # uncontracted clone (like Ewald).
    if is_entangled(sys)
        (; units) = get_entanglement(sys)
        for (ai, int) in enumerate(usys.interactions_union)
            i, _ = unit_and_part(units, ai)
            for pc in int.pair
                pc.isculled && break
                aj = pc.bond.j
                j, _ = unit_and_part(units, aj)
                J = pc.bilin isa Number ? Mat3(pc.bilin*I) : Mat3(pc.bilin)
                phase = cis(2π * dot(q_reshaped, contracted_bond(pc.bond, units).n))
                for α in 1:3, β in 1:3
                    iszero(J[α,β]) && continue
                    accum_pair_dense!(H11, H12, H21, H22, spin_ops[α,ai], spin_ops[β,aj], i, j, J[α,β], J[α,β]*phase, N)
                end
            end
        end
    end

    # H must be hermitian up to round-off errors
    @assert maxdiff(H, H') < 1e-12

    # Make H exactly hermitian
    hermitianpart!(H)

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i,i] += swt.regularization
    end
end


# Calculate y = H*x, where H is the quadratic Hamiltonian matrix (dynamical
# matrix). Note that x is assumed to be a 2D array with first index
# corresponding to q. 
function multiply_by_hamiltonian_SUN!(y::AbstractMatrix{ComplexF64}, x::AbstractMatrix{ComplexF64}, swt::SpinWaveTheory, qs_reshaped::Array{Vec3};
                                      phases=zeros(ComplexF64, size(qs_reshaped)))
    (; sys) = swt

    N = sys.Ns[1]
    Na = natoms(sys.crystal)
    L = (N-1) * Na

    Nq = length(qs_reshaped)
    @assert size(x) == size(y) == (Nq, 2L)
    X = reshape(x, (Nq, N-1, Na, 2))
    Y = reshape(y, (Nq, N-1, Na, 2))
    Y .= 0

    # All operators appearing in interactions have been pre-rotated to local
    # frame
    for (i, int) in enumerate(sys.interactions_union)

        # Onsite coupling, including Zeeman
        op = int.onsite
        for m in 1:N-1, n in 1:N-1
            c = op[m, n] - δ(m, n) * op[N, N]
            @inbounds for q in 1:Nq
                Y[q, m, i, 1] += c * X[q, n, i, 1]
                Y[q, n, i, 2] += c * X[q, m, i, 2]
            end
        end

        # General pair interactions
        for coupling in int.pair
            # Extract information common to bond
            (; isculled, bond) = coupling
            isculled && break

            @assert i == bond.i
            j = bond.j

            map!(phases, qs_reshaped) do q
                cis(2π*dot(q, bond.n))
            end

            # General pair interactions
            for (A, B) in coupling.general.data
                accum_pair_matvec!(Y, X, A, B, i, j, 1, phases, N)
            end
        end
    end

    # Inter-unit bilinear exchange delegated to the uncontracted clone. See the
    # corresponding block in `swt_hamiltonian_SUN!`.
    if is_entangled(sys)
        (; spin_ops) = swt.data
        (; units) = get_entanglement(sys)
        Jqs = zeros(ComplexF64, Nq)
        for (ai, int) in enumerate(uncontracted_system(sys).interactions_union)
            i, _ = unit_and_part(units, ai)
            for pc in int.pair
                pc.isculled && break
                aj = pc.bond.j
                j, _ = unit_and_part(units, aj)
                J = pc.bilin isa Number ? Mat3(pc.bilin*I) : Mat3(pc.bilin)
                map!(phases, qs_reshaped) do q
                    cis(2π * dot(q, contracted_bond(pc.bond, units).n))
                end
                for α in 1:3, β in 1:3
                    iszero(J[α,β]) && continue
                    @. Jqs = J[α,β] * phases
                    accum_pair_matvec!(Y, X, spin_ops[α,ai], spin_ops[β,aj], i, j, J[α,β], Jqs, N)
                end
            end
        end
    end

    if !isnothing(sys.ewald)
        error("Ewald not supported")
    end

    # Add small constant shift for positive-definiteness
    @inbounds @. Y += swt.regularization * X

    nothing
end
