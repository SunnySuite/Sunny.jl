
# The five bilinear block coefficients coupling flavors (m,n) for a pair of
# local-frame operators A (on unit i) and B (on unit j). The N-th slot of a
# rotated operator is the condensate direction Z, so A[m,N] = ⟨f_m|A|Z⟩,
# A[N,N] = ⟨A⟩, etc. This is the single source of truth shared by the dense
# builder and the matrix-free matvec, and by both the general-coupling and
# dipole-sector (Ewald / delegated bilinear) paths.
@inline function swt_pair_coefs(A, B, m, n, N)
    c1 = (A[m,n] - δ(m,n)*A[N,N]) * B[N,N]     # A-block on unit i, ⟨B⟩
    c2 = A[N,N] * (B[m,n] - δ(m,n)*B[N,N])     # ⟨A⟩, B-block on unit j
    c3 = A[m,N] * B[N,n]                        # i→j hop
    c4 = A[N,m] * B[n,N]                        # j→i hop
    c5 = A[m,N] * B[n,N]                        # particle-hole (anomalous)
    return (c1, c2, c3, c4, c5)
end

# Scatter the bilinear contribution of a pair (A on unit i, B on unit j) into the
# dense Hamiltonian blocks. `w_self` scales the intra-unit blocks (c1,c2);
# `w_hop` the inter-unit blocks (c3,c4,c5) — a periodic-wrapping phase for
# ordinary couplings, or a complex dipole coupling J[α,β] for the Ewald sector.
@inline function accum_pair_dense!(H11, H12, H21, H22, A, B, i, j, w_self, w_hop, N)
    for m in 1:N-1, n in 1:N-1
        (c1, c2, c3, c4, c5) = swt_pair_coefs(A, B, m, n, N)
        H11[m,i,n,i] += c1*w_self;  H22[n,i,m,i] += c1*w_self
        H11[m,j,n,j] += c2*w_self;  H22[n,j,m,j] += c2*w_self
        H11[m,i,n,j] += c3*w_hop;   H22[n,j,m,i] += c3*conj(w_hop)
        H11[n,j,m,i] += c4*conj(w_hop);  H22[m,i,n,j] += c4*w_hop
        H12[m,i,n,j] += c5*w_hop;   H12[n,j,m,i] += c5*conj(w_hop)
        H21[n,j,m,i] += conj(c5)*conj(w_hop);  H21[m,i,n,j] += conj(c5)*w_hop
    end
end

# Matrix-free counterpart of `accum_pair_dense!`: accumulate y += (pair block) x
# over all wavevectors, keeping Holstein-Primakoff state in Nambu form
# X[q, flavor, unit, particle/hole]. `w_self` scales intra-unit blocks; `w_hop`
# is a per-q inter-unit weight (typically the periodic-wrapping phase).
@inline function accum_pair_matvec!(Y, X, A, B, i, j, w_hop, N; w_self=1)
    for m in 1:N-1, n in 1:N-1
        (c1, c2, c3, c4, c5) = swt_pair_coefs(A, B, m, n, N)
        @inbounds for q in axes(Y, 1)
            ph = w_hop[q]
            Y[q, m, i, 1] += c1 * w_self * X[q, n, i, 1]
            Y[q, n, i, 2] += c1 * w_self * X[q, m, i, 2]

            Y[q, m, j, 1] += c2 * w_self * X[q, n, j, 1]
            Y[q, n, j, 2] += c2 * w_self * X[q, m, j, 2]

            Y[q, m, i, 1] += c3 * ph * X[q, n, j, 1]
            Y[q, n, j, 2] += c3 * conj(ph) * X[q, m, i, 2]

            Y[q, n, j, 1] += c4 * conj(ph) * X[q, m, i, 1]
            Y[q, m, i, 2] += c4 * ph * X[q, n, j, 2]

            Y[q, m, i, 1] += c5 * ph * X[q, n, j, 2]
            Y[q, n, j, 1] += c5 * conj(ph) * X[q, m, i, 2]
            Y[q, n, j, 2] += conj(c5) * conj(ph) * X[q, m, i, 1]
            Y[q, m, i, 2] += conj(c5) * ph * X[q, n, j, 1]
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

        # For an entangled system, A0 has already had its intra-unit
        # interactions stripped and shifted to an onsite coupling. The same
        # stripping must be done for Aq.
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
                J = gs[ai]' * Aq[ai, aj] * gs[aj] / 2
                J0 = gs[ai]' * A0[ai, aj] * gs[aj] / 2

                for α in 1:3, β in 1:3
                    Ai = spin_ops[α, ai]
                    Bj = spin_ops[β, aj]
                    accum_pair_dense!(H11, H12, H21, H22, Ai, Bj, i, j, J0[α,β], J[α,β], N)
                end
            end
        end
    end

    # Inter-unit bilinear exchange delegated to the uncontracted clone (like
    # Ewald). It lives as `bilin` on the clone's bare bonds rather than in the
    # unit's product-space tensor decomposition. Both factors are embedded spin
    # dipoles, so the contribution routes through the same kernel: w_self = J
    # (coupling at q=0), w_hop = J·phase (coupling at q), with the phase carried
    # by the corresponding contracted bond.
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
                accum_pair_matvec!(Y, X, A, B, i, j, phases, N)
            end
        end
    end

    # Inter-unit bilinear exchange delegated to the uncontracted clone. See the
    # corresponding block in `swt_hamiltonian_SUN!`.
    if is_entangled(sys)
        (; spin_ops) = swt.data
        (; units) = get_entanglement(sys)
        whop = zeros(ComplexF64, Nq)
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
                    @. whop = J[α,β] * phases
                    accum_pair_matvec!(Y, X, spin_ops[α,ai], spin_ops[β,aj], i, j, whop, N; w_self=J[α,β])
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
