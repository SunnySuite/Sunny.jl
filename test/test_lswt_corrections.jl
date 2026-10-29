# Tests of the 1/s corrections to linear spin wave theory, in two tiers. The
# derivation tier, in `test_lswt_corrections_extra.jl`, checks the expansion
# against exact diagonalization and each piece of the derivation against a
# private reimplementation: that is what to reach for when a result is in
# question, rather than what to run on every commit. They must still pass.
# `JuliaTestItems.toml` excludes that file from test discovery, both in VS Code
# and in `@run_package_tests`; delete the `exclude` line to run it.
#
# The default tier pins end-to-end output and checks the invariants that need no
# reference implementation.

# Models and truncated Fock spaces shared by both tiers. A `@testmodule` is
# evaluated once per `run_tests` call, rather than once per test item.
@testmodule CorrectionModels begin
    using Sunny, LinearAlgebra, SparseArrays

    # ---- One generic three-site cluster ----

    # Three sites give three bands, enough for off-diagonal elements of Σ̂; a
    # triclinic cell with generic positions leaves no site symmetry, so every
    # Stevens word is allowed; unequal spins exercise the per-site factors
    # σᵢ/σⱼ; and a generic field cants the moments out of collinearity, without
    # which Σ̂ would vanish identically. Every bond has zero offset, so the
    # Hamiltonian is 𝐪-independent and a single grid point integrates the
    # self-energy exactly. Ferromagnetic exchange is what keeps the Fock space
    # affordable, the anomalous mixing being ⟨n̂⟩ ≈ 0.004 here where a canted
    # antiferromagnet of any `s` would put it near 0.3.
    const cluster_cryst = Crystal(lattice_vectors(1, 1.1, 1.2, 80, 90, 100),
                                  [[0, 0, 0], [0.45, 0.05, 0.1], [0.1, 0.4, 0.05]], 1)
    const cluster_Js = [diagm([-0.6, -0.6, -1.3]), diagm([-1.4, -0.7, -0.7]),
                        -0.1*[1 0.3 -0.2; 0.25 1 0.15; -0.15 0.1 1]]
    const cluster_B = [0.48, -0.32, 0.8]

    # Onsite anisotropy for the cluster, with every Stevens coefficient carrying
    # a factor s^-k so that the classical energy landscape is held fixed as `s`
    # varies, making the rate at which a residual vanishes meaningful. The
    # order-6 word needs s ≥ 3 to be nonzero. `stevens_matrices` carries the
    # spin as a type parameter, so without `@nospecialize` this recompiles for
    # every `s` that any test uses.
    cluster_aniso(@nospecialize(O), s, i) = ((0.3*O[2, 0] + 0.15*O[2, 1])/s^2 + (0.02*O[4, 2] - 0.01*O[4, -3])/s^4 +
                                             0.004*O[6, i]/s^6)

    # Scalar biquadratic exchange on each bond of the cluster, carrying a factor
    # (s_i s_j)^-1 for the reason given above: (𝐒ᵢ⋅𝐒ⱼ)² is of order s⁴, and the
    # bilinear exchange of order s²
    const cluster_Ks = [0.15, -0.1, 0.05]

    # The `aniso` and `biquad` switches are `Bool`s rather than functions or
    # `nothing`, so that this and `cluster_errors` each compile one
    # specialization instead of one per call site.
    function cluster(ss; mode=:dipole, aniso=false, biquad=false)
        sys = System(cluster_cryst, [i => Moment(s=ss[i], g=1) for i in 1:3], mode)
        for (n, (i, j)) in enumerate([(1, 2), (2, 3), (1, 3)])
            set_exchange!(sys, cluster_Js[n], Bond(i, j, [0, 0, 0]); biquad = biquad * cluster_Ks[n] / (ss[i] * ss[j]))
        end
        set_field!(sys, cluster_B)
        aniso && for i in 1:3
            set_onsite_coupling!(sys, cluster_aniso(stevens_matrices(mode == :dipole ? ss[i] : Inf), ss[i], i), i)
        end
        # Deterministic, and converged to a torque of 1e-11 from this start
        polarize_spins!(sys, cluster_B)
        minimize_energy!(sys; jitter=0)
        return sys
    end

    # ---- One generic two-site SU(3) cluster ----

    # Mode :SUN expands in the number of boxes M of the symmetric SU(N)
    # representation rather than in s, and needs its own cluster: two sites
    # rather than three, because the M-box Hilbert space grows as
    # ((M+1)(M+2)/2)^Na. The couplings are generic — an anisotropic bilinear
    # exchange, a biquadratic term that mode :SUN handles with no special
    # treatment, an onsite anisotropy, and a field leaving no symmetry — and
    # every bond offset is zero, as in `cluster`.
    #
    # `pairscale` multiplies the pair coupling. The M-box classical energy is M
    # Σᵢ onsiteᵢ[N,N] + M² Σ A[N,N] B[N,N], and those differing powers make the
    # stationary reference state depend on M: the state to expand about at M
    # boxes is the one Sunny's minimizer finds for `pairscale = M`.
    const sun_cryst = Crystal(lattice_vectors(1, 1.1, 1.2, 80, 90, 100),
                              [[0, 0, 0], [0.45, 0.05, 0.1]], 1)

    function sun_cluster(; pairscale=1, field=cluster_B)
        sys = System(sun_cryst, [i => Moment(s=1, g=1) for i in 1:2], :SUN)
        set_pair_coupling!(sys, (Si, Sj) -> pairscale * (Si'*cluster_Js[1]*Sj + 0.3*(Si'*Sj)^2),
                           Bond(1, 2, [0, 0, 0]))
        set_field!(sys, field)
        S = spin_matrices(1)
        for i in 1:2
            set_onsite_coupling!(sys, 0.35*S[3]^2 + 0.2*(S[1]^2 - S[2]^2), i)
        end
        polarize_spins!(sys, cluster_B)
        minimize_energy!(sys; g_tol=1e-14, jitter=0)
        return sys
    end

    # Operators of the M-box symmetric representation of SU(N), on the basis
    # labeled by the occupations of the Nf = N-1 flavors that have left the
    # condensate, which holds the remaining M - Σn boxes. Returns that basis,
    # the generators E_{mn} = b†_m b_n, and the bosons themselves labeled as
    # `hamiltonian_monomials` labels them, flavor fastest with a ≤ L
    # annihilating. The generators are exact; a boson operator that would leave
    # the space is truncated, so residuals may only be read off blocks of low
    # occupation.
    function box_ops(Nf, M, Na)
        sb = [ns for ns in Iterators.product(ntuple(_ -> 0:M, Nf)...) if sum(ns) <= M]
        pos = Dict(ns => k for (k, ns) in enumerate(sb))
        occ(ns, m) = m <= Nf ? ns[m] : M - sum(ns)     # flavor N is the condensate
        ds = length(sb)
        op(O, i) = reduce(kron, (k == i ? sparse(O) : sparse(1.0I, ds, ds) for k in 1:Na))

        function E(m, n, i)
            A = spzeros(ComplexF64, ds, ds)
            for (k, ns) in enumerate(sb)
                if m == n
                    A[k, k] = occ(ns, m)
                    continue
                end
                iszero(occ(ns, n)) && continue
                ns′ = ntuple(f -> ns[f] + (f == m) - (f == n), Nf)
                all(>=(0), ns′) && sum(ns′) <= M || continue
                A[pos[ns′], k] = sqrt(occ(ns, n) * occ(ns′, m))
            end
            return op(A, i)
        end

        function bop(a)
            L = Nf * Na
            (i, m) = (div(mod1(a, L) - 1, Nf) + 1, mod1(mod1(a, L), Nf))
            A = spzeros(ComplexF64, ds, ds)
            for (k, ns) in enumerate(sb)
                iszero(ns[m]) && continue
                A[pos[ntuple(f -> ns[f] - (f == m), Nf)], k] = sqrt(ns[m])
            end
            return op(a <= L ? A : sparse(A'), i)
        end

        # Total boson number of each index of the Kronecker product
        nbs = [sum(j -> sum(sb[div(k-1, ds^(Na-j)) % ds + 1]), 1:Na) for k in 1:ds^Na]
        return (; dim=ds^Na, E, bop, nbs)
    end

    # ---- Exact diagonalization in a truncated boson Fock space ----

    # Boson operators on a truncated Fock space of `n` sites, labeled by the
    # Nambu index of a `BosonMonomial`: `a ≤ n` annihilates on site `a`, `a > n`
    # creates. Sparse, because the cluster above acts on 7³ states.
    function fock_ops(dims)
        n = length(dims)
        op(O, i) = reduce(kron, (k == i ? O : sparse(1.0I, dims[k], dims[k]) for k in 1:n))
        b(i) = spdiagm(1 => [√float(k) for k in 1:dims[i]-1])
        return a -> ComplexF64.(a <= n ? op(b(a), a) : op(b(a-n)', a-n))
    end

    # Normal-ordered quadratic Hamiltonian, read off from Sunny's own dynamical
    # matrix at 𝐪 = 0, together with that matrix. Any `terms2` are added to it,
    # which is how the anisotropy correction to H₂ enters.
    function fock_quadratic(swt, bop, dim; terms2=nothing)
        L = Sunny.nbands(swt)
        H = zeros(ComplexF64, 2L, 2L)
        Sunny.dynamical_matrix!(H, swt, zero(Sunny.Vec3))
        isnothing(terms2) || Sunny.accum_quadratic!(H, terms2, zero(Sunny.Vec3))
        H2 = spzeros(ComplexF64, dim, dim)
        for i in 1:L, j in 1:L
            H2 .+= ((H[i, j] + H[L+j, L+i])/2) * bop(L+i)*bop(j) +
                   (H[i, L+j]/2) * bop(L+i)*bop(L+j) + (H[L+i, j]/2) * bop(i)*bop(j)
        end
        return (H, H2)
    end

    expand(bop, terms, dim) = sum(t -> t.c * prod(bop, t.as), terms; init=spzeros(ComplexF64, dim, dim))

    # Bogoliubov transformations at each of the given wavevectors, each
    # diagonalized independently, and the vertex they contract to.
    function bogoliubov_matrices(swt, qs)
        L = Sunny.nbands(swt)
        H = zeros(ComplexF64, 2L, 2L)
        return map(qs) do q
            T = zeros(ComplexF64, 2L, 2L)
            Sunny.dynamical_matrix!(H, swt, q)
            Sunny.bogoliubov!(T, H)
            return T
        end
    end
    vertex(swt, terms, qs) = Sunny.vertex!(zeros(ComplexF64, ntuple(_ -> 2Sunny.nbands(swt), length(qs))),
                                           terms, qs, bogoliubov_matrices(swt, qs))

    # Compares the whole Holstein-Primakoff expansion, order by order in the
    # boson number, against exact diagonalization in the *spin* Hilbert space,
    # so that the truncation of the series is itself what is under test.
    function cluster_errors(ss; mode=:dipole, aniso=false, biquad=false)
        sys = cluster(ss; mode, aniso, biquad)
        # Regularization would otherwise leak into the quadratic coefficients
        swt = SpinWaveTheory(sys; measure=nothing, regularization=0)
        Ns = ntuple(i -> Int(2ss[i]+1), 3)
        dim = prod(Ns)
        bop = fock_ops(Ns)
        op(O, i) = reduce(kron, (k == i ? sparse(O) : sparse(1.0I, Ns[k], Ns[k]) for k in 1:3))
        S(a, i) = op(spin_matrices(ss[i])[a], i)

        # Exact cluster Hamiltonian, in the local frames that the boson
        # expansion uses and with the same Zeeman convention (+𝐁⋅𝐒). The
        # anisotropy is rotated by hand, so this is independent of the
        # implementation, which instead rotates the Stevens coefficients.
        Rs = swt.vacuum.frames
        Hex = spzeros(ComplexF64, dim, dim)
        for i in 1:3
            Bi = Rs[i]' * (sys.gs[1, 1, 1, i]' * sys.extfield[1, 1, 1, i])
            Hex .+= sum(a -> Bi[a] * S(a, i), 1:3)
            aniso || continue
            A = Hermitian(Matrix(cluster_aniso(stevens_matrices(ss[i]), ss[i], i)))
            Hex .+= op(Matrix(Sunny.rotate_operator(A, Rs[i])), i)
        end
        # Each bond appears twice, once culled. A scalar biquadratic exchange is
        # stored as the coefficient of Σ_q c_q 𝒪_2^q ⊗ 𝒪_2^q, with c_q the
        # `scalar_biquad_metric`, and in mode :dipole that coefficient carries
        # the renormalization λ₂ of each site, which is divided out here to
        # recover the operator. Mode :dipole has also absorbed the term
        # -(1/2)𝐒ᵢ⋅𝐒ⱼ of (𝐒ᵢ⋅𝐒ⱼ)² into the bilinear exchange stored, and its
        # constant into the scalar.
        for int in sys.interactions_union, c in int.pair
            c.isculled && continue
            @assert iszero(c.bond.n)
            (i, j) = (c.bond.i, c.bond.j)
            J = Rs[i]' * c.bilin * Rs[j]
            Hex .+= sum(J[a, b] * S(a, i) * S(b, j) for a in 1:3, b in 1:3) + c.scalar * I
            iszero(c.biquad) && continue
            λ(s) = mode == :dipole ? Sunny.rcs_factors(s)[2] : 1.0
            Q(i) = [op(Matrix(Sunny.rotate_operator(Hermitian(Matrix(stevens_matrices(ss[i])[2, q])), Rs[i])), i) for q in 2:-1:-2]
            (Qi, Qj) = (Q(i), Q(j))
            K = c.biquad / (λ(ss[i]) * λ(ss[j]))
            Hex .+= K * sum(Sunny.scalar_biquad_metric[q] * Qi[q] * Qj[q] for q in 1:5)
        end

        (; terms2, δE) = Sunny.anisotropy_correction(swt)
        H1 = expand(bop, Sunny.lswt_correction_monomials(swt, Val{1}()), dim)
        H3 = expand(bop, Sunny.cubic_monomials(swt), dim)
        H4 = expand(bop, Sunny.quartic_monomials(swt), dim)
        # Quadratic Hamiltonian of LSWT, plus its own correction at order 1/s
        (_, H2) = fock_quadratic(swt, bop, dim; terms2)

        # Flattened Kronecker index of a boson occupation triple, the states of
        # a given total boson number, and the boson numbers of an index
        idx(ms) = 1 + sum(k -> ms[k] * prod(Ns[k+1:3]), 1:3)
        states(n) = [idx(ms) for ms in Iterators.product(0:n, 0:n, 0:n) if sum(ms) == n]
        ms(a) = ntuple(k -> mod(div(a-1, prod(Ns[k+1:3])), Ns[k]), 3)
        (vac, n1, n2) = (idx((0, 0, 0)), states(1), states(2))
        E0 = real(Hex[vac, vac])
        R = Hex - E0*I - H1 - H2 - H3 - H4

        # With no anisotropy the first four blocks below are *exact*, not merely
        # asymptotic in `s`: a matrix element between one and two bosons is
        # purely cubic, since H₄ conserves boson number modulo two, and every
        # term of H₅ is a transverse operator times the classical Sᶻ of its
        # partner, hence proportional to the transverse field that vanishes at a
        # classical minimum.
        return (; classical = abs(E0 - (energy(sys) + 3δE)) / norm(Hex),
                  quadratic = norm(R[n1, n1]) / norm(Hex),
                  anomalous = norm(R[n2, [vac]]) / norm(Hex),
                  cubic = norm(R[n2, n1]) / norm(Hex),
                  # Unlike the blocks above this one does receive an H₆
                  # contribution, so it is asymptotic even with no anisotropy
                  quartic = norm(R[n2, n2]) / norm(H4[n2, n2]),
                  # Largest residual over every element four bosons can reach,
                  # which is the only available measure once anisotropy makes
                  # each word truncated rather than exact. A structural zero can
                  # never be the maximum, so only the stored entries are
                  # scanned; sweeping all dim² pairs instead costs more than
                  # every Sunny call here combined.
                  reachable = maximum((abs(v) for (a, b, v) in zip(findnz(R)...)
                                       if all(ms(a) .+ ms(b) .<= 4)); init=0.0) / norm(Hex),
                  coherent = max(abs(3δE), norm(H1)) / norm(Hex),
                  hermiticity = (norm(H3 - H3') + norm(H4 - H4')) / (norm(H3) + norm(H4)))
    end

    # The same idea in mode :SUN, whose expansion parameter is the number of
    # boxes M of the symmetric SU(N) representation. Sunny uses M = 1, where a
    # site cannot hold two bosons and the asymptotic blocks are unreachable, so
    # `hamiltonian_monomials` takes M as a test-only argument and this promotes
    # the very same local operators to the M-box representation exactly, via its
    # generators.
    function box_errors(M)
        # The couplings to expand are the unscaled ones, `local_monomials`
        # supplying the powers of M itself; all the scaling is for is the
        # reference state, so carry only the coherents over from the M-box
        # minimization.
        sysM = sun_cluster(; pairscale=M)
        sys0 = sun_cluster()
        sys0.coherents .= sysM.coherents
        swt = SpinWaveTheory(sys0; measure=nothing, regularization=0)
        # The couplings that `local_interactions` rotates into local frames and
        # absorbs the Zeeman term into are what `hamiltonian_monomials` reads.
        sys = swt.sys
        Nf = Sunny.nflavors(swt)
        N = Nf + 1
        (; dim, E, bop, nbs) = box_ops(Nf, M, Sunny.nsites(sys))
        prom(A, i) = sum(A[m, n] * E(m, n, i) for m in 1:N, n in 1:N if !iszero(A[m, n]))

        Hex = spzeros(ComplexF64, dim, dim)
        for i in 1:Sunny.nsites(sys)
            int = Sunny.local_interactions(sys, swt.vacuum.frames, i)
            iszero(int.onsite) || (Hex .+= prom(Matrix(int.onsite), i))
            for c in int.pair
                c.isculled && continue
                @assert iszero(c.bond.n)
                for (A, B) in c.general.data
                    Hex .+= prom(Matrix(A), c.bond.i) * prom(Matrix(B), c.bond.j)
                end
            end
        end

        Hs = [expand(bop, Sunny.hamiltonian_monomials(sys, swt.vacuum.frames, Val{K}(); M), dim) for K in 1:4]
        E0 = real(sum(t -> t.c, Sunny.hamiltonian_monomials(sys, swt.vacuum.frames, Val{0}(); M); init=0.0+0im))
        R = Hex - E0*I - sum(Hs)
        blk(n) = findall(==(n), nbs)
        scale = norm(Hex)

        # The classical, linear, quadratic, anomalous and cubic blocks are all
        # *exact* at the M-box stationary state, not merely asymptotic: the
        # first omitted word is the five-boson -(1/8)A[m,N] b†_m n̂², whose
        # coefficient summed over interactions is the gradient of that energy.
        # The quartic block, which the six-boson word does reach, is what must
        # fall off like 1/M.
        return (; classical = abs(E0 - real(Hex[blk(0)[1], blk(0)[1]])) / scale,
                  linear = norm(R[blk(1), blk(0)]) / scale,
                  quadratic = norm(R[blk(1), blk(1)]) / scale,
                  anomalous = norm(R[blk(2), blk(0)]) / scale,
                  cubic = norm(R[blk(2), blk(1)]) / scale,
                  quartic = norm(R[blk(2), blk(2)]) / norm(Hs[4][blk(2), blk(2)]),
                  hermiticity = (norm(Hs[3] - Hs[3]') + norm(Hs[4] - Hs[4]')) /
                                (norm(Hs[3]) + norm(Hs[4])))
    end

    # Exact retarded Green function of the cluster, in the quasi-particle
    # operators y = T⁻¹x = τ₃T†τ₃x, from which the self-energy follows as Σ̂ = ω
    # - diag(ε) - (Gτ₃)⁻¹, returned for each of the `ωs`. A complex frequency
    # keeps every denominator away from a pole, so the comparison is independent
    # of broadening.
    function cluster_self_energy(H2, Hpert, bop, T0, ε, λ, ωs)
        L = size(T0, 1) ÷ 2
        τ₃ = Diagonal([ones(L); -ones(L)])
        Y = [sum(a -> (τ₃ * T0' * τ₃)[m, a] * bop(a), 1:2L) for m in 1:2L]
        (Es, ψs) = eigen(Hermitian(Matrix(H2 + λ*Hpert)))
        ΔE = Es .- Es[1]
        us = [ψs' * (Y[m]' * ψs[:, 1]) for m in 1:2L]   # ⟨0|y_m|j⟩
        vs = [ψs' * (Y[m] * ψs[:, 1]) for m in 1:2L]    # ⟨0|y_m†|j⟩
        return map(ωs) do ω
            G = [sum(@. us[m]*conj(us[m′])/(ω - ΔE) - conj(vs[m])*vs[m′]/(ω + ΔE))
                 for m in 1:2L, m′ in 1:2L]
            return (ω*I - Diagonal(ε) - inv(G * τ₃)) / λ^2
        end
    end

    # ---- Square lattice ----

    const square_cryst = Sunny.square_crystal(; c=3)

    # Square-lattice antiferromagnet, optionally canted by a field. Sunny's
    # Zeeman coupling is +𝐁⋅𝐒, so the moments cant away from the field, with
    # cos θ = -B/8s. At B = 0 the structure is collinear Néel.
    function canted_square(s, B; mode=:dipole)
        sys = System(square_cryst, [1 => Moment(; s, g=1)], mode)
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        set_field!(sys, [0, 0, B])
        sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
        θ = acos(-B / 8s)
        set_dipole!(sys, [+sin(θ), 0, cos(θ)], (1, 1, 1, 1))
        set_dipole!(sys, [-sin(θ), 0, cos(θ)], (1, 1, 1, 2))
        @assert energy_per_site(sys) ≈ -2s^2 - B^2/16
        return sys
    end

    # Square-lattice antiferromagnet carrying everything at once: a generic
    # field that cants the two sublattices inequivalently, so that Σ̂ is
    # genuinely off-diagonal, and an anisotropy that is diagonal in the global
    # frame but in neither local frame, so that it contributes to every vertex.
    # Its bonds connect distinct cells, which is the one thing the zero-offset
    # `cluster` cannot reach. The ordered state is hard coded:
    # `minimize_energy!` reproduces it only to 1e-9, and an adaptive cubature
    # amplifies that into the last digits of a pinned intensity.
    function anisotropic_square()
        s = 2.0
        sys = System(square_cryst, [1 => Moment(; s, g=1)], :dipole)
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        O = stevens_matrices(s)
        set_onsite_coupling!(sys, (0.4*O[2, 0] + 0.1*O[4, 0] + 0.05*O[4, 4])/s^2, 1)
        set_field!(sys, [1.1, 0.3, 0.7])
        sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
        set_dipole!(sys, [0.8612292608182394, -1.772409186125766, 0.34183305464432445], (1, 1, 1, 1))
        set_dipole!(sys, [-1.082911766624903, 1.6076190866688223, -0.492811300482684], (1, 1, 1, 2))
        @assert energy_per_site(sys) ≈ -8.303053064914755 atol=1e-12
        return sys
    end

    # Easy-axis Néel order on the square lattice, whose gap makes every momentum
    # integral converge exponentially.
    function square_afm(; field)
        sys = System(square_cryst, [1 => Moment(s=1.0, g=1)], :dipole_uncorrected)
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> -0.5*S[3]^2, 1)
        set_field!(sys, field)
        sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
        set_dipole!(sys, [0, 0, +1], (1, 1, 1, 1))
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
        minimize_energy!(sys)
        return sys
    end

    # Static weight of the one-magnon bands, Σ_{n ≤ L} w[n, μ] conj(w[n, ν])
    # with w = T'ũ, contracted through the measure's combiner exactly as
    # `corrected_intensities` contracts the spectral function. With `metric`,
    # the hole bands n > L enter with a minus sign, giving the commutator w'Ĩw.
    # Passing `drop` subtracts |δA|², leaving the cross term 2Re(A conj(δA))
    # that is the correction of relative order 1/s, as the sum rule wants.
    function static_weights(swt, qs, δc; drop=false, metric=false)
        sys = swt.sys
        cryst = Sunny.orig_crystal(sys)
        L = Sunny.nbands(swt)
        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Sunny.num_observables(swt.measure))
        δu = zero(u)
        Ncells = Sunny.nsites(sys) / Sunny.natoms(cryst)
        return map(qs) do q
            q_global = cryst.recipvecs * q
            Sunny.excitations!(T, H, swt, q)
            Sunny.observable_vectors!(u, swt, q_global)
            fill!(δu, 0)
            isnothing(δc) || Sunny.accum_observable_corrections!(δu, swt, q_global, δc)
            w = T' * (u + δu)
            δw = T' * δu
            corr = map(swt.measure.corr_pairs) do (μ, ν)
                c = dot(view(w, 1:L, μ), view(w, 1:L, ν))
                metric && (c -= dot(view(w, L+1:2L, μ), view(w, L+1:2L, ν)))
                drop && (c -= dot(view(δw, 1:L, μ), view(δw, 1:L, ν)))
                c / Ncells
            end
            real(swt.measure.combiner(q_global, corr))
        end
    end

    # ---- Triangular lattice ----

    # Triangular-lattice antiferromagnet in a three-site cell, small enough for
    # the cubic self-energy to be affordable. Isotropic and in zero field, this
    # is the s = 1/2 model of arXiv:0901.4803 in its 120° state. Given an
    # easy-plane anisotropy and a tilted field it instead cants into a state
    # with three inequivalent sublattices, which is what makes the off-diagonal
    # elements of Σ̂ large. The state is built explicitly rather than minimized,
    # so that the chirality of the former and the degenerate direction of the
    # latter are fixed.
    const tri_cryst = Sunny.triangular_crystal(; a=1.0, c=10.0)
    const tri_ds = [[0.120272, -0.449974, -0.181818], [-0.465766, -0.002090, -0.181818],
                    [0.112161, 0.452064, -0.181818]]

    function triangular(; Δ=1.0, field=nothing, φ=0, g=2)
        sys = System(tri_cryst, [1 => Moment(s=1/2, g=g)], :dipole)
        set_exchange!(sys, diagm([1.0, 1.0, Δ]), Bond(1, 1, [1, 0, 0]))
        sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
        R = [cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1]
        if isnothing(field)
            Q = tri_cryst.recipvecs * [1/3, 1/3, 0]
            for site in eachsite(sys)
                θ = dot(Q, global_positions(sys)[site])
                set_dipole!(sys, [cos(θ), sin(θ), 0], site)
            end
            @assert energy_per_site(sys) ≈ -1.5 * (1/2)^2
        else
            set_field!(sys, R * field)
            for i in 1:3
                set_dipole!(sys, R * tri_ds[i], (1, 1, 1, i))
            end
            minimize_energy!(sys)
            @assert energy_per_site(sys) ≈ -0.51131313 atol=1e-8
        end
        return sys
    end
end


@testitem "1/s corrections on a lattice" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: canted_square, square_cryst, sun_cluster, cluster_B

    # Default accuracy for the momentum integrals, enough for the checks below
    # that compare against a reference value: the tightest of them, Oguchi's ζ,
    # comes out 1.4e-7 from its reference here, seven times inside the `atol`
    # asserted. Checks that instead compare two corrections computed from the
    # same quadrature, or extract a ratio from them, loosen it individually.
    tol = 1e-4

    # ---- Oguchi's Z_c ----

    # For a collinear structure the cubic vertex vanishes, so the mean field is
    # the entire O(1/s) shift of the dispersion, and it is a uniform rescaling
    # by Oguchi's Z_c = 1 + ζ/2s [Prog. Theor. Phys. 13, 148 (1960)]. The
    # uniformity in 𝐪 tests the four-boson coefficients, and ζ their overall
    # normalization.
    let s = 1/2
        swt = SpinWaveTheory(canted_square(s, 0); measure=nothing)
        qs = [[0.1, 0, 0], [0.3, 0.17, 0], [0.5, 0, 0], [0.13, -0.4, 0.22]]
        (; terms2, δE) = Sunny.hartree_fock_correction(swt; tol)
        Zc = dispersion(Sunny.swt_with_vacuum(swt; correction=terms2), qs) ./ dispersion(swt, qs)
        @test maximum(abs, Zc .- Zc[1]) < 1e-6
        @test 2s * (Zc[1] - 1) ≈ 0.1579474 atol=1e-6
        @test δE ≈ 0.01247369 atol=1e-7
    end

    # ---- What a collinear structure switches off ----

    # The ordered moment of the square-lattice antiferromagnet at s = 1/2 is the
    # textbook ⟨S⟩ = 0.3034, i.e. a boson density of 0.1966 (QMC gives 0.307).
    # Self-consistency is inert for this collinear antiferromagnet, where the
    # mean field merely rescales H₂ and so leaves the Bogoliubov transformation,
    # hence the mean fields themselves, unchanged. And every cubic monomial
    # carries a transverse component of the exchange in the local frame, so the
    # cubic vertex, the self-energy and the tadpole all vanish identically, the
    # last because the cubic monomials are proportional to a transverse
    # effective field that vanishes at a classical minimum.
    swt = SpinWaveTheory(canted_square(1/2, 0); measure=nothing)
    @test Sunny.boson_density(swt; tol=1e-4)[1] ≈ 0.19656 atol=1e-5
    @test norm(Sunny.corrected_magnetic_moments(swt; tol=1e-4)[1]) ≈ 0.303437 atol=1e-6
    # Para-unitarity of the Bogoliubov transform: ⟨bᵢb†ᵢ⟩ - ⟨b†ᵢbᵢ⟩ = 1
    L = Sunny.nbands(swt)
    o = Sunny.Vec3(0, 0, 0)
    g = Sunny.contractions(swt, [Sunny.BosonMonomial(1.0+0im, as, (o, o)) for as in ((L+1, 1), (1, L+1))], 1e-4)
    @test g(1, L+1, o) - g(L+1, 1, o) ≈ 1 atol=1e-10
    @test maximum(abs, Sunny.cubic_self_energy(swt, [[0.3, 0.1, 0]]; η=0.01, bz=(6, 6, 1))) < 1e-12
    grid32 = (32, 32, 1)
    Zcs = map((Sunny.hartree_fock_correction(swt; bz=grid32).terms2,
               Sunny.to_self_consistent_vacuum(swt; bz=grid32, scf_tol=1e-7).vacuum.correction)) do terms2
        dispersion(Sunny.swt_with_vacuum(swt; correction=terms2), [[0.3, 0.1, 0]])
    end
    @test Zcs[1] ≈ Zcs[2] atol=1e-7
    # The energy of a Gaussian state is stationary in the self-consistent
    # vacuum, against scaling its correction, and lies below that of the
    # harmonic vacuum, which is the harmonic energy plus ⟨H₄⟩ = -δE. At zero
    # field the mean field only rescales H₂ and the two vacua coincide, so the
    # structure is canted.
    let bz = (16, 16, 1), swt = SpinWaveTheory(canted_square(1/2, 1.0); measure=nothing)
        # Re-expanding about a copy of the frames already in `swt` must
        # reproduce `swt` itself: those frames already carry its own gauge
        # twist, so there is no hidden rotation left for the labels to pick up.
        let swt2 = Sunny.swt_with_vacuum(swt; frames=copy(swt.vacuum.frames))
            @test dispersion(swt2, [[0.3, 0.1, 0]]) ≈ dispersion(swt, [[0.3, 0.1, 0]])
            @test swt2.classical_energy ≈ swt.classical_energy
        end
        sc = Sunny.to_self_consistent_vacuum(swt; bz)
        E(λ) = Sunny.gaussian_energy_per_site(Sunny.swt_with_vacuum(swt; correction=[Sunny.BosonMonomial(λ*t.c, t.as, t.ns) for t in sc.vacuum.correction]); bz)
        @test abs(E(1 + 1e-3) - E(1 - 1e-3)) / 2e-3 < 1e-6
        @test E(1) ≈ Sunny.gaussian_energy_per_site(sc) atol=1e-14
        @test E(0) ≈ Sunny.corrected_energy_per_site(swt; bz) - Sunny.hartree_fock_correction(swt; bz).δE atol=1e-12
        @test E(1) < E(0)
        # The self-consistent mean field gaps the Goldstone mode at the ordering
        # wavevector, and a pin restores the harmonic one, which only the
        # regularization gaps
        gap(sw) = minimum(dispersion(sw, [[0, 0, 0]]))
        pinned = Sunny.to_self_consistent_vacuum(swt; bz, constraints=(Sunny.GoldstonePin([0, 0, 0]),))
        @test gap(sc) > 0.1
        @test gap(pinned) ≈ gap(swt) rtol=0.1
        # Relaxing the structure as well lowers the energy further
        st = Sunny.to_self_consistent_vacuum(swt; bz, relax_spins=true)
        @test Sunny.gaussian_energy_per_site(st) < E(1)
    end
    # A field along the moments shifts the magnons by ±B without changing the
    # Bogoliubov transformation. It leaves the Néel state energetically
    # unstable, but its positive-norm modes still define the same vacuum, so
    # the mean field is unchanged.
    let sysB = canted_square(1/2, 0)
        set_field!(sysB, [0.5, 0, 0])
        swtB = SpinWaveTheory(sysB; measure=nothing)
        @test_throws Sunny.InstabilityError dispersion(swtB, [[0.01, 0, 0]])
        hfB = Sunny.hartree_fock_correction(swtB; bz=grid32)
        @test dispersion(Sunny.swt_with_vacuum(swt; correction=hfB.terms2), [[0.3, 0.1, 0]]) ≈ Zcs[1] atol=1e-7
    end
    # The tadpole vanishes term by term, at 1e-34, so `tol` is irrelevant to it
    tad = Sunny.tadpole_correction(swt; tol=1e-3)
    @test all(t -> abs(t.c) < 1e-12, tad.terms2)
    @test abs(tad.δE) < 1e-12
    @test norm(tad.v) < 1e-12

    # `corrected_magnetic_moments` both depletes each moment (|μ| = 1 → 0.887)
    # and tilts it by the tadpole (0.37° away from the field here), and is
    # shaped and signed like `magnetic_moments`, i.e. μ = -g𝐒 indexed by
    # `Site`.
    let
        sys = canted_square(1, 3)
        swt = SpinWaveTheory(sys; measure=nothing)
        corrected = Sunny.corrected_magnetic_moments(swt; tol=1e-5)
        @test size(corrected) == size(magnetic_moments(sys)) == (1, 1, 1, 2)
        @test vec(corrected) ≈ [[-0.8244219, 0, 0.3256618], [0.8244219, 0, 0.3256618]] atol=1e-6
    end

    # Maximal-weight coherent states make this SU(N) system equivalent to the
    # dipole one, so both the boson densities and the corrected moments agree,
    # despite the SU(N) calculation counting N-1 bosons per site instead of one
    # and reading ⟨𝐒⟩ off a general N×N matrix rather than a Holstein-Primakoff
    # series. A loose `tol` suffices because the two modes integrate the same
    # function here, and error cancels in the comparison.
    let
        swt = SpinWaveTheory(canted_square(1, 3; mode=:SUN); measure=nothing)
        swt′ = SpinWaveTheory(canted_square(1, 3); measure=nothing)
        @test Sunny.boson_density(swt; tol=1e-3) ≈ Sunny.boson_density(swt′; tol=1e-3) rtol=1e-5
        @test Sunny.corrected_magnetic_moments(swt; tol=1e-3) ≈
              Sunny.corrected_magnetic_moments(swt′; tol=1e-3) rtol=1e-5
    end

    # ---- Modes :SUN and :dipole_uncorrected on the same s = 1/2 model ----

    # At s = 1/2 an SU(N) system has N = 2 and one boson per site, and its
    # expansion in the box number coincides term by term with the dipole
    # expansion in 1/s. So every correction must agree to machine precision,
    # which pins the :SUN vertices against the dipole ones that exact
    # diagonalization has certified above.
    #
    # Only gauge-invariant output may be compared. `bogoliubov!` fixes the phase
    # of each band independently, and the dipole `local_frames` puts a
    # deliberate rotation into each local frame, so the vertex tensors, the
    # Bogoliubov matrices and the dynamical matrices themselves all differ
    # between the two modes by phases.
    #
    # The tolerance below is loose because `hcubature` decides how finely to
    # subdivide by comparing an error estimate against `tol`, so the two modes
    # can take different refinement paths and leave errors of order `tol` that
    # fail to subtract. Even so, a scaled Hartree-Fock coefficient or a tilted
    # canting angle moves these outputs by far more than 1e-5, so real bugs are
    # still caught.
    let
        (s, B, tol) = (1/2, 0.6, 0.02)
        qs = [[0.23, 0.11, 0], [0.4, 0.3, 0]]
        rs = map((:dipole_uncorrected, :SUN)) do mode
            sys = canted_square(s, B; mode)
            swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
            hf = Sunny.hartree_fock_correction(swt; tol)
            tad = Sunny.tadpole_correction(swt; tol)
            return (; ε = dispersion(swt, qs),
                      E = Sunny.corrected_energy_per_site(swt; tol),
                      n = Sunny.boson_density(swt; tol),
                      δE = [hf.δE, tad.δE],
                      εc = dispersion(Sunny.swt_with_vacuum(swt; correction=[hf.terms2; tad.terms2]), qs),
                      Σ = Sunny.cubic_self_energy(swt, qs; η=0.05, bz=(8, 8, 1)),
                      I = Sunny.corrected_intensities(swt, qs, :one_loop; energies=range(0, 3, 61), η=0.1,
                                                      bz=(26, 26, 1)).data)
        end
        for k in keys(rs[1])
            @test maximum(abs, getfield(rs[1], k) .- getfield(rs[2], k)) < 1e-5
        end
        # Nontrivial: every quantity above is of order unity, and the
        # intensities in particular exercise the full Dyson resummation
        @test rs[1].Σ[1] != 0 && maximum(abs, rs[1].I) > 1
    end

    # ---- Collinearity switches off less in mode :SUN ----

    # Easy-axis Néel order at s = 1, where :SUN carries two flavors per site. The
    # cubic vertex vanishes identically in the dipole mode, every monomial carrying a
    # transverse exchange component, but not in :SUN: a single-ion level may decay
    # into two magnons of *different* flavors, a channel the dipole expansion has no
    # counterpart for. Here that leaves the vertex at 5% of H₂ and gives the upper
    # (single-ion) band a substantial width, while the dipole mode has neither.
    let
        function neel(mode)
            sys = System(square_cryst, [1 => Moment(s=1.0, g=1)], mode)
            set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
            set_onsite_coupling!(sys, S -> -0.5*S[3]^2, 1)
            sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
            set_dipole!(sys, [0, 0, +1], (1, 1, 1, 1))
            set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
            return sys
        end
        swt = SpinWaveTheory(neel(:dipole_uncorrected); measure=nothing)
        swt′ = SpinWaveTheory(neel(:SUN); measure=nothing)
        @test Sunny.cubic_vertex_vanishes(swt, Sunny.cubic_monomials(swt))
        @test !Sunny.cubic_vertex_vanishes(swt′, Sunny.cubic_monomials(swt′))

        # Bands 1-2 are the single-ion excitations at 8.0 and bands 3-4 the magnons.
        # Only the former decay: the two-magnon continuum they sit in starts at twice
        # the 4.39 magnon energy, just below. Damping of the magnons themselves is
        # pure broadening artifact, falling off linearly with η.
        q = [[0.3, 0.1, 0]]
        @test dispersion(swt′, q)[:] ≈ [8, 8, 4.387482, 4.387482] atol=1e-5
        Σs = map(((0.05, 24), (0.025, 48))) do (η, nk)
            Sunny.cubic_self_energy(swt′, q; η, bz=(nk, nk, 1))[:]
        end
        # Overdamped: the width exceeds a tenth of the energy, and is η-independent
        @test all(Σ -> -imag(Σ[1]) > 0.8 * 0.05, Σs)
        @test imag(Σs[1][1]) ≈ imag(Σs[2][1]) rtol=0.1
        @test all(Σ -> -imag(Σ[3]) < 0.01, Σs)
        # Im Σ ≤ 0 on shell: the decay channel damps
        @test all(Σ -> all(<=(1e-12), imag.(Σ)), Σs)
    end

    # ---- Goldstone modes survive every correction ----

    # Near a protected zero mode each correction separately diverges like 1/ε,
    # the Bogoliubov matrix doing so, which makes their cancellation stringent;
    # it is the product ε × shift that must vanish. Two symmetries are gated,
    # and they switch on different corrections.
    #
    # Rotation about the field axis leaves the canted structure a zero mode at
    # 𝐪 = 0, and there all three of mean field, tadpole and cubic self-energy
    # contribute. The Ward identity behind the cancellation holds point by
    # point in the loop wavevector, so it is exact on any grid provided the
    # mean fields are summed on the grid of the self-energy, as
    # `corrected_channels` does; only the regulator η is then left over.
    swt = SpinWaveTheory(canted_square(1, 3); measure=nothing)
    # Without a measure the one-loop pieces still build, and intensities refuse
    @test size(Sunny.OneLoop(swt; η=0.1, bz=(4, 4, 1)).δc) == (2Sunny.nbands(swt), 0)
    @test_throws "No observables" Sunny.corrected_intensities(swt, [[0, 0, 0]], :one_loop; energies=0:0.1:1, η=0.1)
    # A number `bz` scales the automatic loop grid linearly, with the
    # non-dispersing c axis left at one point and a floor of 4 elsewhere. For
    # static corrections a number is rejected in favor of `tol`.
    @test Sunny.auto_bzgrid(swt; η=0.1) == (57, 57, 1)
    @test Sunny.auto_bzgrid(swt; η=0.1, bz=2) == (113, 113, 1)
    @test Sunny.auto_bzgrid(swt; η=5.0) == (4, 4, 1)
    @test_throws "must be positive" Sunny.cubic_self_energy(swt, [[0, 0, 0]]; η=0.1, bz=0)
    @test_throws "use the keyword `tol`" Sunny.boson_density(swt; bz=1e-3)
    let bz = (8, 8, 1)
        t2 = [Sunny.hartree_fock_correction(swt; bz).terms2
              Sunny.tadpole_correction(swt; bz).terms2]
        δ = Sunny.static_self_energy(swt, [[0, 0, 0]], t2)[2]
        @test δ > 1e3
        @test abs(δ + real(Sunny.cubic_self_energy(swt, [[0, 0, 0]]; η=1e-8, bz)[2])) < 1e-8 * δ
    end
    # The on-shell bands take the cancelling static part at real ω = 0, so a
    # soft band does not move with η, even at ε ≈ η. At ε + iη it moved by 25%.
    let sys = canted_square(1/2, 1), swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
        soft(η) = minimum(Sunny.corrected_channels(swt, [[0.005, 0.005, 0]], :one_loop_on_shell;
                                                   energies=Float64[], η, bz=(24, 24, 1)).bands.disp)
        @test soft(0.05) ≈ soft(0.02) rtol=1e-3
    end

    # The second symmetry is broken only by an onsite anisotropy. Every term
    # below is invariant under rotation about ẑ, so the in-plane moment of this
    # easy-plane ferromagnet leaves a zero mode at 𝐪 = 0, while the cubic and
    # linear vertices vanish identically by that same symmetry. That leaves
    # `anisotropy_correction` to cancel the mean field on its own -- which it
    # can only do because `hamiltonian_monomials` keeps exactly one order in
    # 1/s per word. Since a perturbation of the quadratic form opens a gap like
    # its square root, expanding the anisotropy exactly instead, thereby
    # injecting a partial set of O(1/s²) terms, would gap the mode at O(1/s).
    for mode in (:dipole, :dipole_uncorrected)
        s = 2.0
        sys = System(Sunny.square_crystal(; c=3), [1 => Moment(; s, g=1)], mode)
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        O = stevens_matrices(mode == :dipole ? s : Inf)
        set_onsite_coupling!(sys, (0.5*O[2, 0] - 0.2*O[4, 0]/s^2) / (3s^2), 1)
        set_dipole!(sys, [1, 0, 0], (1, 1, 1, 1))
        sw = SpinWaveTheory(sys; measure=nothing)
        @test maximum(abs(t.c) for t in Sunny.cubic_monomials(sw)) < 1e-12
        @test maximum(abs(t.c) for t in Sunny.lswt_correction_monomials(sw, Val{1}()); init=0.0) < 1e-12
        ε = dispersion(sw, [[0, 0, 0]])[1]
        mf = Sunny.hartree_fock_correction(sw; tol).terms2
        resid(t2) = Sunny.static_self_energy(sw, [[0, 0, 0]], t2)[1] * ε
        @test abs(resid([mf; Sunny.anisotropy_correction(sw).terms2])) < 1e-7 < abs(resid(mf))
    end
end

@testitem "1/s corrections to intensities" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: triangular, square_afm, square_cryst, canted_square, static_weights

    sys = triangular()
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)

    # ---- Harmonic results ----

    # The energy per site is -0.5388 J, and the magnon energy at the M point of
    # the original lattice is 2Js, the lowest of the three folded bands.
    @test Sunny.corrected_energy_per_site(swt; tol=1e-4) ≈ -0.53881 atol=1e-4
    q = [[1/2, 0, 0]]
    @test dispersion(swt, q)[:] ≈ [√2.5, √2.5, 1] atol=1e-6

    # ---- The cubic self-energy on a lattice ----

    # Being a symmetric energy minimum, the 120° structure cannot be tilted by
    # zero-point fluctuations. Unlike the collinear case the cubic monomials are
    # individually nonzero, so their cancellation here tests their relative
    # phases.
    tad = Sunny.tadpole_correction(swt; tol=1e-3)
    @test maximum(t -> abs(t.c), tad.terms2) < 1e-6
    @test abs(tad.δE) < 1e-12

    terms2 = [Sunny.hartree_fock_correction(swt; tol=1e-3).terms2; tad.terms2]
    δ = Sunny.static_self_energy(swt, q, terms2)[:]
    Σs = map(nk -> Sunny.cubic_self_energy(swt, q; η=0.02, bz=(nk, nk, 1))[:], (24, 48))

    # The self-energy converges like 1/nk, so a Richardson step gives the O(1/s)
    # magnon energy at the M point. It falls 27% below the harmonic value, most
    # of that coming from the cubic self-energy rather than the mean-field
    # shift.
    εs = [1 + δ[3] + real(Σ[3]) for Σ in Σs]
    @test 2εs[2] - εs[1] ≈ 0.7316 atol=2e-3

    # Because the harmonic dispersion vanishes at Γ and K, the lower edge of the
    # two-magnon continuum touches the one-magnon branch at every wavevector,
    # and a magnon acquires a width only where the branch lies strictly inside
    # the continuum. There the width survives η → 0, as it does for the top of
    # the band here, where 2Γ/ε extrapolates to about 0.2, of the order of the
    # maximum ~0.3 that the reference reports. The M-point magnon instead sits
    # on the boundary, and its apparent width is entirely the Lorentzian tail of
    # the regularization, falling off like η.
    Σ = Sunny.cubic_self_energy(swt, q; η=0.01, bz=(48, 48, 1))[:]
    @test imag(Σ[3]) ≈ imag(Σs[2][3]) / 2 rtol=0.01
    @test imag(Σ[1]) / imag(Σs[2][1]) > 0.85

    # ---- The corrected spectral function ----

    # The result is a retarded response, so its frequency integral is the
    # commutator w'Ĩw, for any self-energy and on any grid, and it is positive
    # at ω > 0 where the auxiliary model is stable. Both hold at s = 1/2. The
    # spectrum sums its mean fields on its loop grid, so the reference is built
    # from the same grid.
    bz = (12, 12, 1)
    swt2 = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
    qs2 = [[0.476, 0, 0], [0.375, 0.125, 0]]

    δc = Sunny.observable_corrections(swt2; v=Sunny.tadpole_correction(swt2; bz).v, bz)
    refs = static_weights(swt2, qs2, δc)
    comms = static_weights(swt2, qs2, δc; metric=true)

    # A Lorentzian tail needs range rather than resolution, so the window is wide and
    # the step is a fraction of η.
    η = 0.06
    energies = range(-20, 24, 1501)
    chans = Sunny.corrected_channels(swt2, qs2, :one_loop; energies, η, bz)
    (; transverse) = chans
    pos = energies .> 0
    @test !any(chans.breakdown)
    @test all(≥(0), (transverse + chans.cross + chans.direct)[pos, :])
    # The commutator is small for a trace measure, so it is compared on the scale of
    # the static weight. The residual is the truncated Lorentzian tail.
    @test maximum(abs, vec(sum(transverse; dims=1)) * step(energies) - comms) < 2e-3 * maximum(refs)
    # About 10% of the transverse weight at ω > 0 moves to the continuum and to
    # the negative-frequency poles
    @test 0.8 < minimum(vec(sum(transverse[pos, :]; dims=1)) * step(energies) ./ refs) < 0.95

    # The public `SelfEnergy` is the propagator of `:one_loop`, whose spectral
    # matrix it must reproduce
    let z = energies[800] + im*η, Σ = Sunny.SelfEnergy(Sunny.OneLoop(swt2; η, bz), qs2[1])
        Ĩ = Diagonal([ones(L); -ones(L)])
        G = inv(z * Ĩ - Diagonal(abs.(Σ.ε)) - Σ(z))
        c = Sunny.corrected_channels(swt2, qs2[1:1], :one_loop; energies=[real(z)], η, bz, spectral=true)
        @test ((G' - G) / (2π*im))[1:L, 1:L] ≈ c.specfunc[:, :, 1, 1]
    end

    # The reproduction schemes propagate the resolvent of a stable quadratic
    # model here, and assemble the magnon term and the continuum separately, so
    # both are positive. The magnon term of `:one_loop_projected` keeps the
    # amplitudes and metric, so it satisfies the commutator sum rule;
    # `:one_loop_on_shell` keeps a unit-weight Lorentzian per pole.
    for scheme in (:one_loop_projected, :one_loop_on_shell)
        c = Sunny.corrected_channels(swt2, qs2, scheme; energies, η, bz, spectral=true)
        @test all(≥(0), (c.transverse + c.cross + c.direct)[pos, :])
        if scheme == :one_loop_projected
            @test maximum(abs, vec(sum(c.transverse; dims=1)) * step(energies) - comms) < 2e-3 * maximum(refs)
        else
            # Every pole decays forward in time, by construction
            @test all(≥(0), c.bands.widths)
            # The residual is the Lorentzian tails outside the window
            nb = size(c.specfunc, 1)
            weights = [real(tr(c.specfunc[:, :, iω, iq])) for iω in eachindex(energies), iq in eachindex(qs2)]
            @test vec(sum(weights; dims=1)) * step(energies) ≈ fill(nb, length(qs2)) rtol=1e-2
        end
    end

    # The interference is invisible to a trace measure integrated over the zone:
    # 𝐒·𝐒 is a scalar, so its expansion in bosons has no term linking an odd
    # number of them to an even one, and the cancellation is exact for every ω
    # once Σ_𝐪 restores momentum conservation. It is the only check here that
    # constrains the *relative phase* of the two routes to a pair, a phase that
    # each channel alone is free of. The window must reach below ω = 0, or a
    # pair at x ≈ 0 contributes just half of its Lorentzian.
    qs4 = vec([[i, j, 0] ./ 3 for i in 0:2, j in 0:2])
    chans4 = Sunny.corrected_channels(swt2, qs4, :one_loop; energies=range(-4, 12, 161), η=0.3,
                                      bz=(3, 3, 1))
    @test abs(sum(chans4.cross)) < 1e-3 * sum(chans4.direct)

    # ---- Gauge invariance, where Σ̂ is genuinely off-diagonal ----

    # With three inequivalent sublattices the off-diagonal elements of Σ̂ are a
    # third of the diagonal, whereas for the 120° structure above they vanish
    # identically: its branches sit in momentum sectors that the cubic vertex
    # cannot connect, as do those of the umbrella that a field along z produces.
    # Only the exact-diagonalization cluster of the derivation tier constrains
    # them otherwise. They are also the one part of Σ̂ sensitive to the per-band
    # phase that `bogoliubov!` fixes independently, since a gauge conjugation Σ̂
    # → DΣ̂D† with D diagonal and unitary leaves the eigenvalues of the Dyson
    # denominator alone, and with them the poles, and leaves τ₃Σ̂ Hermitian;
    # only the contraction against the observable amplitudes T'u sees it.
    # Rotation about z is an exact symmetry of this model, fixing both diagm([1,
    # 1, Δ]) and the trace structure factor while moving every local frame, so
    # the corrected intensities must be invariant under it. Breaking that gauge
    # violates this by 5%. The loop grid is coarse because an invariance holds
    # grid by grid and needs no converged integral.
    #
    # The rotated side is evaluated at -𝐪 and negative frequency, which the
    # retarded response relates to +𝐪 by χ''(𝐪, ω) = -χ''(-𝐪, -ω). The structure is
    # non-reciprocal, S(𝐪, ω) and S(-𝐪, ω) differing by 50%, so this also pins
    # the hole block of the Dyson equation to the energies and eigenvectors of
    # -𝐪 rather than to mirror images of those at 𝐪. Taking the bare hole
    # energies from ε(𝐪) violates it by 20%, a mistake the rotation alone and
    # every reciprocal model are blind to.
    swts = map(φ -> let sys3 = triangular(; Δ=0.6, field=[0.7, 0, 1.2], φ, g=1)
                        SpinWaveTheory(sys3; measure=ssf_trace(sys3; apply_g=false))
                    end, (0.0, 0.9))
    qs3 = [[0.23, 0.11, 0], [0.37, 0.09, 0]]
    Σ3 = Sunny.cubic_self_energy_matrix(swts[1], qs3[1], dispersion(swts[1], qs3[1:1])[1:1] .+ 0.05im;
                                        η=0.05, bz=(6, 6, 1))[1:L, 1:L, 1]
    @test maximum(abs, Σ3 - Diagonal(diag(Σ3))) > 0.2 * maximum(abs, diag(Σ3))
    energies = range(0.2, 2.0, 25)
    for scheme in (:one_loop, :two_magnon_ladder, :one_loop_projected, :one_loop_on_shell)
        spectral(swt, qs, energies) = let c = Sunny.corrected_channels(swt, qs, scheme; energies, η=0.15, bz=(4, 4, 1))
            c.transverse + c.cross + c.direct
        end
        data = spectral(swts[1], qs3, energies)
        mirror = spectral(swts[2], -qs3, -reverse(energies))
        @test data ≈ -reverse(mirror; dims=1) rtol=1e-6
        @test Sunny.corrected_intensities(swts[1], qs3, scheme; energies, η=0.15, bz=(4, 4, 1)).data ≈ data
    end

    # ---- Quantum sum rule on the square lattice ----

    # A field-polarized ferromagnet conserves Sᶻ. Its Bogoliubov transformation
    # is trivial, so the two-magnon channel must carry no weight at all, and
    # nothing else in the 1/s expansion is nonzero for this state either: every
    # mean field vanishes in the empty vacuum, and a collinear structure has no
    # cubic vertex. So the corrected intensities must reduce to those of linear
    # spin wave theory, which is itself exact here, the polarized state and its
    # one-magnon excitations being eigenstates. Being a retarded response, the
    # result also carries the mirror pole of -𝐪 at negative frequency, whose
    # tail reaches ω > 0.
    let
        sys = System(square_cryst, [1 => Moment(s=1, g=2)], :dipole)
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        set_field!(sys, [0, 0, 0.5])
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 1))
        @assert energy_per_site(sys) ≈ -2 - 1.0
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        qs = [[0.3, 0.2, 0]]
        (energies, η) = (range(0, 10, 101), 0.2)
        chans = Sunny.corrected_channels(swt, qs, :one_loop; energies, η, bz=(4, 4, 1))
        @test maximum(abs, chans.direct) < 1e-25
        kernel = lorentzian(fwhm=2η)
        mirror = reverse(intensities(swt, -qs; energies=-reverse(energies), kernel).data)
        @test chans.transverse + chans.cross + chans.direct ≈
              intensities(swt, qs; energies, kernel).data - mirror atol=1e-12

        # Expanding instead about a vacuum whose quadratic Hamiltonian is
        # shifted by a uniform b†b term leaves its eigenvectors alone, and with
        # no loop to dress, the counterterm restores the harmonic propagator
        # exactly.
        o = Sunny.Vec3(0, 0, 0)
        L = Sunny.nbands(swt)
        shift = [Sunny.BosonMonomial(0.3+0im, (L+1, 1), (o, o))]
        c = Sunny.corrected_channels(Sunny.swt_with_vacuum(swt; correction=shift), qs, :one_loop; energies, η, bz=(4, 4, 1))
        @test c.transverse + c.cross + c.direct ≈ chans.transverse + chans.cross + chans.direct atol=1e-12
        @test c.disp ≈ chans.disp .+ 0.3
    end

    # The fast path of `corrected_channels`. A collinear structure has no cubic
    # vertex, but its monomials cancel on merging rather than being absent, so
    # the gate must be on magnitude: the coefficients of the 26 monomials here
    # sum to 1e-12 of the quadratic Hamiltonian, against 5e-2 once a field cants
    # the same model. Its observable consequence is that `cross` comes out
    # identically zero rather than at the 1e-24 of the round-off it replaces,
    # the rest of the result being untouched — checked to be bit-identical,
    # which the two tests below cannot see.
    let
        # Largest cubic coefficient relative to the quadratic Hamiltonian, the
        # dimensionless quantity the gate thresholds
        function vertex_scale(swt)
            L = Sunny.nbands(swt)
            H = zeros(ComplexF64, 2L, 2L)
            Sunny.dynamical_matrix!(H, swt, Sunny.Vec3(0, 0, 0))
            terms3 = Sunny.cubic_monomials(swt)
            @test length(terms3) == 26
            return maximum(abs(t.c) for t in terms3) / norm(H)
        end

        sys = square_afm(; field=[0.0, 0, 0])
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        @test vertex_scale(swt) < 1e-11
        canted = SpinWaveTheory(square_afm(; field=[1.5, 0, 0]); measure=nothing)
        @test vertex_scale(canted) > 1e-2

        chans = Sunny.corrected_channels(swt, [[0.3, 0.2, 0]], :one_loop; energies=range(0, 12, 121),
                                         η=0.2, bz=(8, 8, 1))
        @test iszero(chans.cross) && !iszero(chans.direct)
    end

    # The bare ladder lacks the Adler zero of a Goldstone mode, so on a fine
    # enough grid (56² or more) it binds soft pairs below ω = 0 at 𝐪 = 0 and at
    # the ordering wavevector of the Néel state. Check that a breakdown is
    # detected at these wavevectors, but not at their neighbors.
    let
        sys = canted_square(1/2, 0)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        qs = [[0, 0, 0], [1/2, 1/2, 0], [0.05, 0, 0], [1/2, 0.45, 0]]
        chans = Sunny.corrected_channels(swt, qs, :two_magnon_ladder; energies=[1.0], η=0.1, bz=(64, 64, 1))
        @test chans.breakdown == [true, true, false, false]

        # A breakdown is marked by NaN in every component of a measure that
        # returns a vector
        swt = SpinWaveTheory(sys; measure=ssf_custom((q, ssf) -> diag(ssf), sys))
        res = @test_logs (:warn, r"will not resolve") Sunny.corrected_intensities(swt, qs[1:1], :two_magnon_ladder;
                                                                                  energies=[-0.05, 0.05], η=0.1, bz=(64, 64, 1))
        @test all(x -> all(isnan, x), res.data)
    end
end


@testitem "1/s corrections, end to end" setup=[CorrectionModels] begin
    using .CorrectionModels: cluster, anisotropic_square

    # Two models pinned through the public entry point, so that any change to
    # any part of the 1/s machinery moves a number here. Between them they carry
    # inequivalent sublattices, an off-diagonal Σ̂, a canted structure with a
    # nonzero tadpole, anisotropy words at two Stevens orders, a biquadratic
    # exchange, unequal spins, nonzero bond offsets, a Stevens observable, and a
    # triclinic cell with no site symmetry — there is no invariant being
    # checked, only reproducibility. The two Stevens orders come from
    # `anisotropic_square`, since the k > 2 words of `cluster_aniso` vanish
    # identically at the small spins that keep the cluster cheap.
    #
    # The mean fields are summed on the loop grid passed below, so no adaptive
    # cubature enters, and the pair measures are binned at η/10, tight enough
    # here that the binning is invisible at the `rtol` asserted.
    # The ordered states are deterministic, either from an explicit
    # `polarize_spins!` start or hard coded, so a rebuild reproduces both sums
    # bit-for-bit.
    #
    # Re-pinned when the mean fields moved onto the loop grid. The cluster is
    # 𝐪-independent, so one point already integrated them exactly and only
    # `:one_loop_on_shell` moved there, by its change of definition. On the
    # square lattice every scheme moved by about 1%, the error of the mean
    # fields on an 8×8 grid; that is the price of summing them where the
    # self-energy is summed, and it is small against η. `:one_loop_on_shell` on
    # the cluster was re-pinned again, by 3e-6, when its static part moved to
    # ω = 0.
    # Each `scheme` is pinned, as (sum, maximum) of the intensities.
    function pinned(swt, qs, energies, η, dims, scheme)
        res = Sunny.corrected_intensities(swt, qs, scheme; energies, η, bz=dims)
        return [sum(res.data), maximum(res.data)]
    end

    let sys = cluster((1.0, 3/2, 1.0); aniso=true, biquad=true)
        # Zero bond offsets make every vertex 𝐪-independent, so a single loop
        # wavevector integrates the self-energy exactly. The quadrupole is
        # mostly a pair creator, so its words of every rank enter; the schemes
        # differ only downstream of those words, so one suffices.
        quadrupole = Sunny.ssf_multipoles(sys; max_k=2) do q, C
            real(C[2, 1, 2, 1])
        end
        refs = [ssf_trace(sys; apply_g=false) => (; one_loop = [134.67891417203998, 8.349860176352514],
                                                    one_loop_projected = [134.1001430057365, 8.330843191193743],
                                                    one_loop_on_shell = [134.03117398331528, 8.269158004258141]),
                quadrupole => (; one_loop = [16.22226535928757, 0.35843088787581123])]
        for (measure, refs) in refs
            swt = SpinWaveTheory(sys; measure)
            for (scheme, ref) in pairs(refs)
                @test pinned(swt, [[0.21, -0.33, 0.12], [0.5, 0, 0]], range(0.1, 6.0, 120), 0.1, (1, 1, 1), scheme) ≈ ref rtol=1e-6
            end
        end
    end

    let sys = anisotropic_square()
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        refs = (; one_loop = [79.02848991053109, 3.8296352863687693],
                  one_loop_projected = [74.58444179829294, 3.146640845289903],
                  one_loop_on_shell = [74.7926152607717, 3.1485933083002178])
        for (scheme, ref) in pairs(refs)
            @test pinned(swt, [[0.3, 0.1, 0], [0.37, 0.22, 0]], range(0.2, 14.0, 200), 0.2, (8, 8, 1), scheme) ≈ ref rtol=1e-6
        end
    end
end
