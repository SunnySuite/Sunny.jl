# The derivation tier of the 1/s corrections tests, see the top of
# `test_lswt_corrections.jl`, which also defines the `CorrectionModels` used here.

@testitem "1/s corrections against exact diagonalization" setup=[CorrectionModels] begin
    using LinearAlgebra, SparseArrays
    using .CorrectionModels: cluster_errors, box_errors, anisotropic_square, bogoliubov_matrices, vertex

    # ---- The boson expansion against the exact cluster Hamiltonian ----

    # `cluster_errors` compares the whole Holstein-Primakoff expansion, order by
    # order in the boson number, against exact diagonalization in the *spin*
    # Hilbert space, so that the truncation of the series is itself what is
    # under test.

    # Exchange and field only. Two `s` triples: one certifies every exact block,
    # the second supplies the 1/s scaling of the leftover H₆ piece, which would
    # instead approach unity were the quartic term wrong at O(s⁰).
    rs = map(((1.0, 3/2, 1.0), (2.0, 3.0, 2.0))) do ss
        err = cluster_errors(ss)
        @test err.classical < 1e-12
        @test err.quadratic < 1e-8
        @test err.anomalous < 1e-12
        @test err.cubic < 1e-8
        @test err.hermiticity < 1e-12
        return err.quartic
    end
    @test 0.4 < rs[2] / rs[1] < 0.6

    # Now `cluster_aniso` on every site. Here the words are truncated rather
    # than exact — one order in 1/s is kept per word — so what is left over is
    # the next order. Requiring it to vanish faster than 1/s, the size of the
    # retained correction itself, pins every retained coefficient: an error at
    # the order kept would leave a residual of the same size as the correction,
    # and duplicating a term that LSWT already holds would leave one of order
    # unity. Mode :dipole renormalizes the stored Stevens coefficients and
    # :dipole_uncorrected does not, so both are checked.
    for mode in (:dipole, :dipole_uncorrected)
        # Every spin tuple here is `NTuple{3,Float64}`, matching the
        # exchange-only case above; mixing in an integer tuple would specialize
        # `cluster_errors` a second time, which costs more in compilation than
        # the whole block does in arithmetic.
        # A biquadratic exchange is the product of two quadrupoles, whose words
        # are truncated likewise, so it is covered by the same checks.
        e = cluster_errors((3.0, 3.0, 3.0); mode, aniso=true, biquad=true)
        @test e.hermiticity < 1e-12
        @test e.reachable < 0.02
        # The Stevens coefficients are renormalized in mode :dipole so that the
        # classical energy function is exact in a spin coherent state, to all
        # orders in 1/s. That makes the corrections to the energy and to the
        # linear term vanish identically, and leaves the anomalous words as the
        # only correction to LSWT's H₂.
        @test mode == :dipole ? e.coherent < 1e-12 : e.coherent > 1e-3
        @test mode == :dipole ? max(e.classical, e.quadratic) < 1e-12 : e.quadratic > 1e-5
        # The sub-leading words of the quadrupoles enter the anomalous block,
        # which the biquadratic exchange alone isolates. With them its residual
        # falls about sevenfold from s = 2 to 4, and without them only threefold.
        (r2, r4) = map(s -> cluster_errors((s, s, s); mode, biquad=true).anomalous, (2.0, 4.0))
        @test r4 < 0.2 * r2
        # The 1/s falloff of that leftover is checked in one mode only, the word
        # truncation being the same in both.
        mode == :dipole || continue
        @test cluster_errors((4.0, 4.0, 4.0); mode, aniso=true, biquad=true).reachable < 0.5 * e.reachable
    end


    # ---- The :SUN boson expansion against the M-box Hamiltonian ----

    # The same idea in mode :SUN, whose expansion parameter is the number of
    # boxes M of the symmetric SU(N) representation. Sunny uses M = 1, where a
    # site cannot hold two bosons and the asymptotic blocks are unreachable, so
    # `hamiltonian_monomials` takes M as a test-only argument and this promotes
    # the very same local operators to the M-box representation exactly, via its
    # generators.

    # Two values of M give the scaling of the quartic residual, which would
    # instead approach a constant were the quartic term wrong at leading order.
    quartics = map((4, 8)) do M
        err = box_errors(M)
        @test err.classical < 1e-12
        @test err.linear < 1e-12
        @test err.quadratic < 1e-12
        @test err.anomalous < 1e-12
        @test err.cubic < 1e-12
        @test err.hermiticity < 1e-12
        return err.quartic
    end
    @test 0.4 < quartics[2] / quartics[1] < 0.6


    # ---- A single ion, where both modes are exactly solvable ----

    # Mode :dipole is a symplectic restriction of :SUN: the coherent states it
    # explores are the SU(2) orbit of the maximal-weight state inside the full
    # projective space, and RCS renormalizes the Stevens coefficients so that
    # the classical energy is exactly ⟨n̂|Ĥ|n̂⟩ on that orbit. A single ion with
    # no couplings makes both statements checkable against the N × N spectrum.
    #
    # Stevens operators grow like sᵏ, so the higher orders are divided by sᵏ to
    # keep the easy axis dominant and the ground state at m = s for every s.
    ion_aniso(s) = let O = stevens_matrices(s)
        -O[2,0]/s^2 + 0.08O[4,0]/s^4 - 0.03O[6,0]/s^6
    end

    function single_ion(s, mode, B)
        cryst = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0, 0, 0]], 1)
        sys = System(cryst, [1 => Moment(; s, g=1)], mode)
        set_onsite_coupling!(sys, ion_aniso(s), 1)
        set_field!(sys, B)
        polarize_spins!(sys, [0, 0, 1])
        minimize_energy!(sys; g_tol=1e-15, jitter=0)
        return sys
    end

    # The exact Hamiltonian of that ion as an N × N matrix. The Zeeman term is
    # read back from `sys.extfield` so that the unit conversion of `set_field!`
    # needs no duplicating here.
    function ion_hamiltonian(s, sys)
        S = spin_matrices(s)
        return Matrix(ion_aniso(s)) + sum(sys.extfield[1][d] * Matrix(S[d]) for d in 1:3)
    end

    # Stevens words of order k vanish for 2s < k, so: s = 1 keeps k = 2 alone, s
    # = 2 adds k = 4, and s = 3 adds k = 6.
    for s in (1, 2, 3)
        # Along -ẑ, the Zeeman energy being +𝐁⋅(g𝐒), so that m = s is the minimum
        B = [0, 0, -0.3]
        # Axial, so the exact levels are already labeled by m = s, s-1, …, -s
        H = ion_hamiltonian(s, single_ion(s, :dipole, B))
        @test norm(H - Diagonal(diag(H))) < 1e-12
        lv = real(diag(H))
        gaps = sort(lv)[2:end] .- minimum(lv)

        # In :SUN the condensate is an arbitrary N-vector, so the quadratic form
        # is the exactly projected Hamiltonian: all N-1 bands come out exact,
        # and the expansion terminates. The cubic words are pure round-off and
        # there are no quartic ones at all, so every 1/M correction vanishes
        # identically.
        swt = SpinWaveTheory(single_ion(s, :SUN, B); measure=nothing, regularization=0)
        @test sort(dispersion(swt, [[0, 0, 0]])[:]) ≈ gaps
        @test maximum(abs, [t.c for t in Sunny.cubic_monomials(swt)]; init=0.0) < 1e-12
        @test isempty(Sunny.quartic_monomials(swt))
        # Exactly zero at any tolerance, there being no term left to integrate
        @test iszero(Sunny.hartree_fock_correction(swt; tol=0.1).δE)
        @test iszero(Sunny.boson_density(swt; tol=0.1))

        # Mode :dipole keeps one band, and RCS makes it the exact gap to the
        # level one unit of magnetization down — *not* the smallest gap, which
        # for an easy axis is the nearly degenerate m = -s partner. Only m = s-1
        # survives at linear order in u = sin²(θ/2), since |⟨m|θ⟩|² ∝ u^(s-m);
        # the 2s of the expansion cancels the 1/s relating curvature to
        # frequency, so the identity is s-independent and holds at every Stevens
        # order at once.
        swt = SpinWaveTheory(single_ion(s, :dipole, B); measure=nothing, regularization=0)
        @test dispersion(swt, [[0, 0, 0]])[1] ≈ lv[2] - lv[1]
    end

    # That identity is not generic: it needs the axial symmetry that pins n̂.
    # Adding a transverse field tilts the moment, and then the classical minimum
    # of ⟨n̂|Ĥ|n̂⟩ is no longer where the quantum gradient vanishes — a coherent
    # state not being an eigenstate — leaving ⟨s-1|Ĥ|s⟩ ≠ 0 in the local frame
    # and the one-band gap merely approximate. Mode :SUN, whose variational
    # family is the whole projective space, stays exact.
    let s = 2, B = [0.4, 0, -0.3]
        sys = single_ion(s, :SUN, B)
        swt = SpinWaveTheory(sys; measure=nothing, regularization=0)
        lv = sort(real(eigvals(Hermitian(ion_hamiltonian(s, sys)))))
        @test sort(dispersion(swt, [[0, 0, 0]])[:]) ≈ lv[2:end] .- lv[1]
    end


    # ---- Symmetries of the vertex at every slot count ----

    # Invariances that need no reference tensor, checked on the one model whose
    # bonds connect distinct cells, so that the Fourier phases the zero-offset
    # cluster above cannot reach are exercised. The contraction against an
    # explicit reference tensor is in the derivation tier.
    swt = SpinWaveTheory(anisotropic_square(); measure=nothing)
    L = Sunny.nbands(swt)
    q(x, y) = Sunny.Vec3(x, y, 0)
    cases = ((Sunny.cubic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.54, -0.22))),
             (Sunny.quartic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.22, 0.35), q(-0.32, -0.57))))
    for (terms, qs) in cases
        K = length(qs)
        Ts = bogoliubov_matrices(swt, qs)
        buf() = zeros(ComplexF64, ntuple(_ -> 2L, K))
        U = Sunny.vertex!(buf(), terms, qs, Ts)

        # Permuting the (momentum, band) slots together leaves the vertex
        # invariant
        for p in Sunny.slot_permutations(K)
            Up = Sunny.vertex!(buf(), terms, ntuple(t -> qs[p[t]], K), ntuple(t -> Ts[p[t]], K))
            @test permutedims(Up, invperm(collect(p))) ≈ U
        end

        # Hermiticity of the underlying operator implies U(-𝐤; n̄) = conj(U(𝐤;
        # n)), where n̄ = n + L exchanges creation and annihilation. Only the
        # magnitudes are compared, because `bogoliubov!` fixes the phase of each
        # band independently at 𝐤 and -𝐤.
        Um = vertex(swt, terms, ntuple(t -> -qs[t], K))
        bar(n) = mod1(n + L, 2L)
        @test abs.(U) ≈ [abs(Um[CartesianIndex(map(bar, Tuple(i)))]) for i in CartesianIndices(U)]
    end
end


@testitem "1/s corrections: sum rules (derivation)" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: square_afm, fock_ops, fock_quadratic, expand, static_weights

    # Weights of the three channels into which the quantum sum rule decomposes,
    # all per site and in units where a trace measure is used, so that no local
    # frame projection survives.
    #
    # Because Sᶻ = s - b†b is exact in the local frame, the two-magnon spectrum
    # saturates the longitudinal sum rule ⟨(δSᶻ)²⟩ = n(1+n) + |Δ|², where n =
    # ⟨b†b⟩ and Δ = ⟨bb⟩ follow from Wick's theorem. The transverse channel
    # obeys an identity sharper still: because the truncated S⁺ = σ(b - b†bb/4s)
    # makes S⁻S⁺ = n̂(2s+1-n̂) exact, and because completeness turns a sum of
    # one-magnon weights over bands and wavevectors into the static expectation
    # value ⟨Â†Â⟩, the one-magnon bands must carry ⟨(Sˣ)² + (Sʸ)²⟩ = s + 2s⟨n̂⟩
    # - ⟨n̂²⟩. Linear spin wave theory produces only the first two terms; the
    # -⟨n̂²⟩ is supplied entirely by `observable_corrections`, so this fixes
    # both the sign and the magnitude of that correction. `tol` is set by the
    # two identities below, which come out 3e-8 from exact here and 4e-7 at
    # `tol=1e-5`. The sum rule itself is limited instead by the energy integral,
    # whose 3.7e-4 does not move with `tol` at all.
    function channel_weights(sys; nq=16, tol=1e-6)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        L = Sunny.nbands(swt)

        # Onsite ⟨b†b⟩ and ⟨bb⟩, from which ⟨n̂²⟩ = ⟨n̂⟩² + ⟨n̂⟩(1+⟨n̂⟩) + |Δ|²
        o = Sunny.Vec3(0, 0, 0)
        g = Sunny.contractions(swt, [Sunny.BosonMonomial(1.0+0im, as, (o, o)) for i in 1:L for as in ((L+i, i), (i, i))],
                               tol)
        ss = [swt.sys.κs[i] for i in 1:L]
        n = [real(g(L+i, i, o)) for i in 1:L]
        n2 = [n[i]^2 + n[i] * (1 + n[i]) + abs2(g(i, i, o)) for i in 1:L]

        # A uniform grid offset by half a step cancels the phase factors of all
        # correlations at distances below `nq`, leaving only the onsite ones
        # above.
        qs = vec([[(a - 0.5)/nq, (b - 0.5)/nq, 0] for a in 1:nq, b in 1:nq])
        δc = Sunny.observable_corrections(swt; tol)
        harm = sum(static_weights(swt, qs, nothing)) / length(qs)
        transverse = sum(static_weights(swt, qs, δc; drop=true)) / length(qs)

        # Elastic weight of the ordered moment, shortened by zero-point
        # fluctuations. In dipole mode the shortening is the boson density `n`
        # already gathered above.
        elastic = sum(i -> (ss[i] - n[i])^2, 1:L) / L

        # Two-magnon continuum, integrated over energy. Its wavevector average
        # converges quickly enough to use a coarser grid, which matters because
        # each point requires its own momentum-space integral. For a one-atom
        # chemical cell this average is the longitudinal weight per site.
        qs2 = vec([[(a - 0.5)/3, (b - 0.5)/3, 0] for a in 1:3, b in 1:3])
        energies = range(-2, 16, 181)
        direct = Sunny.corrected_channels(swt, qs2, :one_loop; energies, η=0.2, bz=(12, 12, 1)).direct
        longitudinal = sum(direct) * step(energies) / length(qs2)

        return (; harm, transverse, elastic, longitudinal,
                harm_ref = sum(@. ss + 2ss*n) / L,
                transverse_ref = sum(@. ss + 2ss*n - n2) / L,
                casimir = sum(@. ss * (ss + 1)) / L)
    end

    # Field causes canting, such that the onsite ⟨bb⟩ is nonzero.
    let w = channel_weights(square_afm(; field=[1.5, 0, 0]))
        @test abs(w.harm / w.harm_ref - 1) < 1e-7
        @test abs(w.transverse / w.transverse_ref - 1) < 1e-7

        # The quantum sum rule, and the point of the whole exercise. Because
        # 𝐒⋅𝐒 is a Casimir, the elastic weight of the ordered moment, the
        # one-magnon bands and the two-magnon continuum must together carry
        # exactly s(s+1), and each is produced by a different part of this
        # module. Linear spin wave theory saturates the rule only through O(s):
        # using its uncorrected one-magnon weights instead overshoots by ⟨n̂²⟩,
        # which is the entire O(s⁰) content of the rule, and some 3% of s(s+1)
        # here. The residual error is that of the energy integral above.
        @test abs((w.elastic + w.transverse + w.longitudinal) / w.casimir - 1) < 1e-3
        @test (w.elastic + w.harm + w.longitudinal) / w.casimir - 1 > 0.02
    end

    # The sum rule, order by order. The weight of a correlation integrated over
    # ω > 0 equals the static weight ⟨Q†Q⟩ - |⟨Q⟩|², an identity that holds
    # separately at each order of the 1/s expansion. Scale every 1/s piece by a
    # bookkeeping parameter ε: the cubic vertex by ε; the static self-energy,
    # the amplitude corrections δc and δpair by ε²; and, because a K-boson word
    # of an observable is smaller than its one-boson word by ε^(K-1), its pair
    # amplitude by ε. Then the ε² coefficient of the integrated corrected
    # spectrum must equal that of the strict static weight, the ground state of
    # H₂ + εH₃ + ε²(H₄ + δH₂ + H₁) with Q = Σ_K ε^(K-1) W_K, computed in a
    # truncated Fock space. That is exact, not asymptotic, so the gate needs a
    # single s and is limited only by Fock truncation and finite differences.
    #
    # The cluster is a canted antiferromagnetic dimer, so that H₃ ≠ 0 and every
    # channel is live: a uniform field cants the moments, and a staggered
    # transverse field and easy axis pin the canting, keeping the ground state
    # gapped. Without the pin it tunnels between its two chiralities, with a gap
    # that closes exponentially in s, and is not the state that the expansion
    # describes. The observables span ranks 1 to 3, among them a quadrupole
    # whose spectrum is mostly two-magnon and an octupole whose one-magnon
    # amplitude is suppressed, cases where the leading channel is not the
    # one-magnon one.
    let s = 4.0, nmax = 10
        cryst = Crystal(lattice_vectors(1, 1.1, 10, 90, 90, 100), [[0, 0, 0], [0.4, 0.3, 0.1]], 1)
        sys = System(cryst, [1 => Moment(; s, g=1), 2 => Moment(; s, g=1)], :dipole)
        set_exchange!(sys, 1/s, Bond(1, 2, [0, 0, 0]))
        O = stevens_matrices(s)
        set_onsite_coupling!(sys, -0.1 * O[2, 2] / s, 1)
        set_onsite_coupling!(sys, -0.1 * O[2, 2] / s, 2)
        sys = to_inhomogeneous(sys)
        set_field_at!(sys, [-0.25, 0, 0.8], (1, 1, 1, 1))
        set_field_at!(sys, [0.25, 0, 0.8], (1, 1, 1, 2))
        set_dipole!(sys, [0.9, 0, -0.4], (1, 1, 1, 1))
        set_dipole!(sys, [-0.9, 0, -0.4], (1, 1, 1, 2))
        minimize_energy!(sys; maxiters=10_000)

        # Gauss-Legendre nodes, by Golub-Welsch
        function gausslegendre(n)
            (x, V) = eigen(SymTridiagonal(zeros(n), [k / √(4k^2 - 1) for k in 1:n-1]))
            return (x, 2 .* V[1, :].^2)
        end
        (xs, ws) = gausslegendre(100)
        # Taylor coefficients of f at 0, through second order
        function taylor2(f; h=0.02)
            (fm2, fm1, f0, fp1, fp2) = f.((-2h, -h, 0.0, h, 2h))
            return (; f0, c1 = (fm2 - 8fm1 + 8fp1 - fp2) / 12h, c2 = (-fm2 + 16fm1 - 30f0 + 16fp1 - fp2) / 24h^2)
        end
        scaled(ts, x) = [Sunny.BosonMonomial(x * t.c, t.as, t.ns) for t in ts]

        # Weight at ω > 0 of the :one_loop response at q = 0, a finite sum of
        # poles χ(z) = Σ rₙ/(z - ωₙ). It is (Σ rₙ + Σ rₙ sign ωₙ)/2, the second
        # sum being -(1/π) ∫ χ(iy) dy, so that no broadening enters.
        function corrected_weight(swt, ε)
            L = Sunny.nbands(swt)
            ol = Sunny.OneLoop(swt; η=1e-3, bz=(1, 1, 1))
            ol = Sunny.OneLoop(ol.swt, ol.η, ol.loop_grid, ol.bin_width, scaled(ol.terms3, ε), ol.terms4,
                               scaled(ol.terms2, ε^2), ε^2 * ol.δc, map(ts -> scaled(ts, ε^2), ol.δpair))
            se = Sunny.SelfEnergy(ol, [0, 0, 0])
            # The row of the bath for the pair amplitude, by ε
            d = 2L + 1
            for ρ in (se.decay, se.source), v in values(ρ.bins), n in 1:ρ.dim
                v[n ≤ d ? Sunny.packed_index(n, d) : Sunny.packed_index(d, n)] *= (n == d ? ε^2 : ε)
            end
            p = 1:2L
            E = Diagonal(abs.(se.ε))
            Ĩ = Diagonal([ones(L); -ones(L)])
            w = se.w[:, 1]
            function χ(zs)
                (; Σ, K) = Sunny.dyson_model(se, :one_loop, zs)
                return map(eachindex(zs)) do i
                    Kz = K[:, :, i]
                    G = inv(zs[i] * Ĩ - E - Σ - Kz[p, p])
                    w' * G * w + (Kz[d:d, p] * G * w)[1] + (w' * G * Kz[p, d:d])[1] + Kz[d, d] + (Kz[d:d, p] * G * Kz[p, d:d])[1]
                end
            end
            θ = (π/2) .* xs
            signed = -real(sum(ws .* (π/2) ./ cos.(θ).^2 .* χ(im .* tan.(θ)))) / π
            total = real(im * 1e8 * only(χ([1e8im])))
            return (total + signed) / 2
        end

        function strict_weight(swt, ε)
            L = Sunny.nbands(swt)
            Ns = (nmax+1, nmax+1)
            dim = prod(Ns)
            bop = fock_ops(Ns)
            (_, H2) = fock_quadratic(swt, bop, dim)
            H3 = expand(bop, Sunny.cubic_monomials(swt), dim)
            H4 = expand(bop, [Sunny.quartic_monomials(swt); Sunny.lswt_correction_monomials(swt, Val{2}())], dim) +
                 expand(bop, Sunny.lswt_correction_monomials(swt, Val{1}()), dim)
            Q = sum(K -> ε^(K-1) * expand(bop, reduce(vcat, (Sunny.observable_monomials(swt, 1, i, Val{K}()) for i in 1:L)), dim), 1:4)
            ψ = eigen(Hermitian(Matrix(H2 + ε * H3 + ε^2 * H4))).vectors[:, 1]
            return real(ψ' * Q' * Q * ψ) - abs2(ψ' * Q * ψ)
        end

        for A in (spin_matrices(s)[1], spin_matrices(s)[3], O[2, 1], O[2, 2], O[3, -1])
            st = Sunny.StevensExpansion(Sunny.matrix_to_stevens_coefficients(Hermitian(ComplexF64.(Matrix(A)))))
            measure = Sunny.MeasureSpec(fill(st, 1, 1, 1, 1, 2, 1), [(1, 1)], (q, corr) -> real(only(corr)),
                                        fill(one(FormFactor), 1, 2, 1))
            swt = SpinWaveTheory(sys; measure, regularization=0)
            a = taylor2(ε -> corrected_weight(swt, ε))
            b = taylor2(ε -> strict_weight(swt, ε))
            @test a.f0 ≈ b.f0 rtol=1e-7
            @test abs(a.c1) < 1e-8 * b.f0
            @test a.c2 ≈ b.c2 rtol=1e-6
        end

        # For an observable with no one-magnon part, the ε² term above is its
        # leading order, and its first correction is the ε⁴ term. Here the
        # corrections do NOT satisfy the gate: the pair amplitude is taken in
        # the harmonic vacuum, missing the shift of ⟨W₂†W₂⟩ that the vacuum
        # takes at first order. The quadrupole 𝒪₂² of a collinear Néel dimer
        # isolates it, since there H₃ = 0 and 𝒪₂² creates only pairs.
        #
        # Measured, the ε⁴ coefficient is -4.0668 here and -2.5739 strict. The
        # strict value splits into the words of the observable, -4.0668, which
        # the code has exactly, and the vacuum shift, +1.4929, which it lacks.
        # That shift is in turn a Gaussian part, +1.714, from the mean field of
        # H₄ acting on the two lines of the pair, and a four-magnon admixture,
        # -0.2215, from H₄ coupling the pair to two further magnons.
        #
        # Planned, and deferred: the four-magnon admixture needs a double loop
        # over momenta, and it must come with the Gaussian part. Taking the
        # Gaussian part alone, by building the bath from the magnons of H₂ plus
        # the static mean field, closes it on this dimer exactly, but breaks the
        # Casimir sum rule of the canted square lattice (1.6%, against < 1e-3)
        # and leaves the dressed H₂ indefinite near Γ. Partial ε⁴ terms appear
        # to spoil identities that hold at the order kept, so the fix is all or
        # nothing. In a canted state the second order of H₃ enters too.
        #
        # The same order also takes the tilt of the pair amplitude, one
        # displaced leg of the three-boson words of the observable, by analogy
        # with the tilt of `observable_corrections`; see
        # `observable_pair_corrections`. Its amplitude is O(ε³), so it reaches
        # the weight only at ε⁴, and it vanishes here, where there is no
        # tadpole. Gating it needs a pure pair creator about a canted state,
        # e.g. 𝒪₂² in the local frame of each site of the canted dimer above.
        # Once all three pieces are in, this becomes `@test`.
        sys = System(cryst, [1 => Moment(; s, g=1), 2 => Moment(; s, g=1)], :dipole)
        set_exchange!(sys, 1/s, Bond(1, 2, [0, 0, 0]))
        set_onsite_coupling!(sys, -0.15 * O[2, 0] / s, 1)
        set_onsite_coupling!(sys, -0.15 * O[2, 0] / s, 2)
        sys = to_inhomogeneous(sys)
        set_field_at!(sys, [0, 0, -0.3], (1, 1, 1, 1))
        set_field_at!(sys, [0, 0, 0.3], (1, 1, 1, 2))
        set_dipole!(sys, [0, 0, 1], (1, 1, 1, 1))
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
        st = Sunny.StevensExpansion(Sunny.matrix_to_stevens_coefficients(Hermitian(ComplexF64.(Matrix(O[2, 2])))))
        measure = Sunny.MeasureSpec(fill(st, 1, 1, 1, 1, 2, 1), [(1, 1)], (q, corr) -> real(only(corr)),
                                    fill(one(FormFactor), 1, 2, 1))
        swt = SpinWaveTheory(sys; measure, regularization=0)
        # Fourth-order Taylor coefficient of f at 0, from a symmetric stencil
        function taylor4(f; h=0.05)
            v = f.((-3h, -2h, -h, 0.0, h, 2h, 3h))
            return (; c2 = (2v[1] - 27v[2] + 270v[3] - 490v[4] + 270v[5] - 27v[6] + 2v[7]) / 360h^2,
                      c4 = (-v[1] + 12v[2] - 39v[3] + 56v[4] - 39v[5] + 12v[6] - v[7]) / 144h^4)
        end
        a = taylor4(ε -> corrected_weight(swt, ε))
        b = taylor4(ε -> strict_weight(swt, ε))
        @test a.c2 ≈ b.c2 rtol=1e-7
        @test_broken a.c4 ≈ b.c4 rtol=1e-4
    end
end


@testitem "1/s corrections: boson expansion (derivation)" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: cluster, fock_ops, fock_quadratic, expand, cluster_self_energy, vertex

    # ---- Mean fields, self-energy and pair amplitudes, one Fock space ----

    sys = cluster((1.0, 3/2, 1.0))
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)
    nmax = 6
    dim = (nmax+1)^L
    bop = fock_ops(ntuple(_ -> nmax+1, L))
    (H, H2) = fock_quadratic(swt, bop, dim)
    ψ = eigen(Hermitian(Matrix(H2))).vectors[:, 1]
    expect(O) = dot(ψ, O, ψ)

    # The moments cant by up to 20°, which is what makes Σ̂ nonvanishing below
    ds = [normalize(sys.dipoles[1, 1, 1, i]) for i in 1:L]
    @test maximum(norm(ds[i] × ds[j]) for i in 1:L, j in 1:L) ≈ 0.35 atol=0.02

    # Truncation is harmless only if the vacuum has no weight in the top sector
    ψn = reshape(ψ, ntuple(_ -> nmax+1, L))
    @test maximum(i -> norm(selectdim(ψn, i, nmax+1)), 1:L) < 1e-4

    # Every Nambu correlation, including the anomalous ones
    o = Sunny.Vec3(0, 0, 0)
    tol = 1e-8
    g = Sunny.contractions(swt, [Sunny.BosonMonomial(1.0+0im, (a, a′), (o, o)) for a in 1:2L for a′ in 1:2L], tol)
    @test [g(a, a′, o) for a in 1:2L, a′ in 1:2L] ≈ [expect(bop(a)*bop(a′)) for a in 1:2L, a′ in 1:2L] atol=1e-7

    # One contraction of H₄ gives the mean-field form, two the constant
    terms4 = Sunny.quartic_monomials(swt)
    g = Sunny.contractions(swt, terms4, tol)
    terms2 = Sunny.wick_reduce(terms4, g, nothing, Val{2}())
    δE = -Sunny.wick_expectation(terms4, g, nothing; nfactors=2)
    H4 = expand(bop, terms4, dim)
    H4mf = expand(bop, terms2, dim)

    # Wick's theorem is exact in a Gaussian state, so the constant subtracted by
    # the decoupling is precisely -⟨H₄⟩. This pins the three pairings, their
    # sign, and the fact that the constant is removed once rather than twice.
    @test -δE ≈ expect(H4) atol=1e-7

    # Defining property of the decoupling: the mean-field operator reproduces
    # the response of H₄ to every quadratic perturbation.
    Qs = [bop(a)*bop(a′) for a in 1:2L, a′ in 1:2L]
    @test [expect(H4*Q - Q*H4) for Q in Qs] ≈ [expect(H4mf*Q - Q*H4mf) for Q in Qs] atol=1e-6

    # The analogous decoupling of H₃, which Wick-contracts down to the linear
    # operator that tadpole relaxation must cancel. Only linear perturbations
    # test anything here, since a Gaussian state gives ⟨[H₃, Q]⟩ = 0 for
    # quadratic Q; a commutator of two linear operators is a c-number, so this
    # pins the three cubic pairings exactly rather than to integration accuracy.
    terms3 = Sunny.cubic_monomials(swt)
    g = Sunny.contractions(swt, terms3, tol)
    ℓ = Sunny.nambu_vector(Sunny.wick_reduce(terms3, g, nothing, Val{1}()), L, 1e-6)
    H3raw = expand(bop, terms3, dim)
    H3mf = sum(a -> ℓ[a] * bop(a), 1:2L)
    @test [expect(H3raw*bop(a) - bop(a)*H3raw) for a in 1:2L] ≈
          [expect(H3mf*bop(a) - bop(a)*H3mf) for a in 1:2L] atol=1e-8
    # The linear part must be removed from the perturbation below: normal
    # ordering H₃ in the quasi-particle basis leaves that piece behind, and it
    # is the tadpole correction rather than a loop.
    H3 = H3raw - H3mf

    τ₃ = Diagonal([ones(L); -ones(L)])
    T0 = zeros(ComplexF64, 2L, 2L)
    ε = copy(Sunny.bogoliubov!(T0, H))

    # The whole Nambu matrix of the cubic self-energy, against the exact
    # retarded Green function. This fixes the off-diagonal elements, which mix
    # bands and correct spectral weights, and the anomalous blocks, which admix
    # three magnons into the ground state. Nothing cheaper constrains the
    # off-diagonal ones: a collinear structure has Σ̂ = 0 outright, and the
    # branches of a spiral or an umbrella live in momentum sectors that the
    # cubic vertex cannot mix, so Σ̂ is diagonal there too. They are also the
    # one part of Σ̂ sensitive to the per-band phase that `bogoliubov!` fixes
    # independently; leaving it free corrupts them by 70%, two orders of
    # magnitude above this tolerance.
    ωs = [0.5 + 0.3im, -1.1 + 0.25im]
    Σed = cluster_self_energy(H2, H3, bop, T0, ε, 0.06, ωs)
    Σm = Sunny.cubic_self_energy_matrix(swt, [0, 0, 0], ωs; η=1e-10, bz=(1, 1, 1))
    @test collect(eachslice(Σm; dims=3)) ≈ Σed atol=3e-4
    @test maximum(abs, Σm[:, :, 1] - Diagonal(diag(Σm[:, :, 1]))) > 0.2 * maximum(abs, diag(Σm[:, :, 1]))

    # Below the three-magnon threshold, unitarity requires τ₃Σ̂ to be Hermitian,
    # which is what makes the Dyson equation preserve spectral weight. The
    # broadening η is what breaks it, by an amount η ∂Σ/∂ω.
    Σh = Sunny.cubic_self_energy_matrix(swt, [0, 0, 0], [0.5 + 1e-10im]; η=1e-10, bz=(1, 1, 1))[:, :, 1]
    @test τ₃ * Σh ≈ (τ₃ * Σh)' atol=1e-9

    # Second-order perturbation theory in H₃, evaluated exactly in the truncated
    # Fock space, against the on-shell `cubic_self_energy`, which takes no
    # frequencies. The vacuum shift constrains the source channel alone, and
    # the level shifts the diagonal of Σ̂ at the on-shell frequency. Three
    # magnons are created and destroyed in the former, and the unrestricted sum
    # over their bands supplies 3! orderings, cancelling one of the two factors
    # of 3! that relate the symmetrized vertex to Γ₂. Levels are identified by
    # proximity to the harmonic band energy, the Fock spectrum interleaving
    # two-magnon states among the one-magnon ones, so a band permutation would
    # be caught here. The residual is the O(h²) error of the differencing.
    function levels(λ)
        E = eigen(Hermitian(Matrix(H2 + λ*H3))).values
        ΔE = E .- E[1]
        return [E[1]; [ΔE[argmin(abs.(ΔE .- ε[n]))] for n in 1:L]]
    end
    shifts = ((levels(0.01) + levels(-0.01))/2 - levels(0)) / 0.01^2
    Σ = Sunny.cubic_self_energy(swt, [[0, 0, 0]]; η=1e-10, bz=(1, 1, 1))
    @test real(vec(Σ)) ≈ shifts[2:L+1] atol=1e-3
    @test maximum(abs, imag(Σ)) < 1e-8
    U3 = vertex(swt, terms3, ntuple(_ -> zero(Sunny.Vec3), 3))
    @test shifts[1] ≈ -6 * sum(abs2(U3[L+n1, L+n2, L+n3]) / (ε[n1] + ε[n2] + ε[n3])
                           for n1 in 1:L, n2 in 1:L, n3 in 1:L) atol=1e-6

    # Both amplitudes for the observable to create a pair of magnons, against
    # exact matrix elements in the same Fock space. This is what pins their
    # relative phase, on which the interference of `corrected_channels` depends
    # and nothing else does: each squared amplitude is invariant under a phase
    # on either one separately, so the two diagonal blocks of the measure — the
    # self-energy and the `direct` channel — are blind to an error here. So are
    # the sum rule and the rotation-invariance checks elsewhere, both being
    # homogeneous in the interference. A sign error was found this way,
    # `pair_amplitude` having returned the amplitude of +b†b where Sᶻ = s - b†b
    # carries a minus.
    #
    # The comparison needs a two-magnon state and the operator that reaches it.
    # For a normalized pair |ab⟩ = y_a†y_b†|0⟩/√(1+δ_ab), each of Sunny's
    # amplitudes is conj(⟨ab|·|·⟩) up to the combinatorial factor √(2/(1+δ_ab))
    # by which the symmetrized state differs from Sunny's sum over ordered
    # pairs. The readout is the full S^{μν} tensor, so that every component is
    # reachable and nothing is protected by the scalar sum rule that a trace
    # measure enjoys.
    swt2m = SpinWaveTheory(sys; measure=ssf_custom((q, ssf) -> ssf, sys; apply_g=false))
    Nobs = Sunny.num_observables(swt2m.measure)
    u2 = zeros(ComplexF64, 2L, Nobs)
    q2m = [0.23, -0.41, 0.17]
    qr2 = Sunny.to_supercell_rlu(sys, q2m)
    qg2 = Sunny.orig_crystal(sys).recipvecs * q2m
    Sunny.observable_vectors!(u2, swt2m, qg2)
    words2m = Sunny.observable_pair_monomials(swt2m, qg2)

    # Observable A_ν(q) as a linear form in the Nambu vector. The index swap is
    # forced, not chosen: `intensities_bands` forms Avec[μ] = dot(u[:,μ],
    # T[:,n]) for the left amplitude ⟨0|A†|n⟩, so A = Σ_a u[ā] x_a. The
    # `contractions` block above fixes the state, so nothing here is free,
    # and the transverse amplitude must come out equal to Sunny's own w = T†u.
    Aodd = [sum(a -> u2[Sunny.nambu_conj(a, L), ν] * bop(a), 1:2L) for ν in 1:Nobs]
    # The even part, whose words `observable_pair_monomials` supplies. They
    # carry the conjugated Fourier phase, being the amplitude to create a pair,
    # so they are conjugated back to describe the same operator as `Aodd`.
    Aeven = [sum(w -> conj(w.c) * bop(w.as[1]) * bop(w.as[2]), words2m[ν]) for ν in 1:Nobs]

    Y = [sum(a -> (τ₃ * T0' * τ₃)[m, a] * bop(a), 1:2L) for m in 1:2L]
    ψ1 = [Y[m]' * ψ for m in 1:L]
    @test [dot(ψ1[n], Aodd[ν] * ψ) for n in 1:L, ν in 1:Nobs] ≈ (T0' * u2)[1:L, :] atol=1e-7

    U3m = zeros(ComplexF64, 2L, 2L, 2L)
    T3m = [conj(T0[Sunny.nambu_conj(a, L), Sunny.nambu_conj(b, L)]) for a in 1:2L, b in 1:2L]
    Sunny.vertex!(U3m, terms3, ntuple(_ -> zero(Sunny.Vec3), 3), (T0, T0, T3m), similar(U3m))
    w2m = T0' * u2

    # Each amplitude separately, then both at once, which is the quantity the
    # interference actually depends on. A phase convention shared by both would
    # cancel from each separately and survive in the coherent sum, so that sum
    # is compared against first-order perturbation theory for ⟨ab|A_ν|0⟩: the
    # even part reaches the pair directly, the odd part through a one-magnon
    # intermediate.
    #
    # The same matrix elements, accumulated into the measure that
    # `corrected_channels` actually consumes, certify it in its own storage
    # convention rather than as rederived here. Linear bin splitting preserves
    # the zeroth and first moments of the measure exactly, so those two moments
    # pin it without replicating the binning. The 12 block is the interference,
    # which no published calculation constrains.
    ρ2m = Sunny.pair_measures(swt2m, terms3, qr2, Sunny.LoopGrid([zero(Sunny.Vec3)], [1.0], 1);
                              bin_width=1e-3, words2=words2m).decay
    rows = [1:L; 2L+1:2L+Nobs]
    bin(v) = [Sunny.bin_entry(v, n, n′) for n in rows, n′ in rows]
    M0 = sum(bin(v) for v in values(ρ2m.bins))
    M1 = sum(j * ρ2m.Δ * bin(v) for (j, v) in ρ2m.bins)
    (E0, E1) = (zeros(ComplexF64, L+Nobs, L+Nobs), zeros(ComplexF64, L+Nobs, L+Nobs))
    for a in 1:L, b in 1:L
        ψpair = (Y[a]' * (Y[b]' * ψ)) / sqrt(1 + (a == b))
        c = sqrt(2 / (1 + (a == b)))
        x = ε[a] + ε[b]
        # Magnon-mediated route: the pair is reached from one magnon through H₃.
        # Fixes the √18 and the external leg.
        med = [dot(ψpair, H3 * ψ1[m]) for m in 1:L]
        @test med ≈ c * [conj(√18 * U3m[a, b, Sunny.nambu_conj(m, L)]) for m in 1:L] atol=1e-6
        # Direct route, including the sign of Sᶻ = s - b†b
        dir = [dot(ψpair, Aeven[ν] * ψ) for ν in 1:Nobs]
        @test dir ≈ c * [conj(Sunny.pair_amplitude(words2m[ν], T0, T0, a, b)) for ν in 1:Nobs] atol=1e-6
        # Both routes at once
        mediated(ν) = sum(m -> conj(w2m[m, ν]) * √18 * U3m[a, b, Sunny.nambu_conj(m, L)] / (x - ε[m]), 1:L)
        both = (dir + [sum(m -> med[m] * dot(ψ1[m], Aodd[ν] * ψ) / (x - ε[m]), 1:L) for ν in 1:Nobs]) ./ c
        @test both ≈ [conj(Sunny.pair_amplitude(words2m[ν], T0, T0, a, b) + mediated(ν)) for ν in 1:Nobs] atol=1e-6
        y = conj([med; dir]) ./ c
        E0 .+= y * y'
        E1 .+= x .* (y * y')
    end
    @test M0 ≈ E0 atol=1e-6
    @test M1 ≈ E1 atol=1e-5

    # The measures `corrected_channels` resums, against the exact retarded
    # response χ^{μν}(z) = Σ_j ⟨0|A_μ†|j⟩⟨j|A_ν|0⟩/(z - E_j) -
    # ⟨0|A_ν|j⟩⟨j|A_μ†|0⟩/(z + E_j) of the cluster to H₂ + λH₃, with observables
    # A_odd + λA_even. Its O(λ²) part must be w'G₀K_mm G₀w + K_dm G₀w + w'G₀K_md
    # + K_dd, which exercises both channels, all 2L Nambu legs, and the relative
    # orientation of every block. Frequencies of both signs reach the decay and
    # source channels resonantly in turn.
    (; decay, source) = Sunny.pair_measures(swt2m, terms3, qr2, Sunny.LoopGrid([zero(Sunny.Vec3)], [1.0], 1);
                                            bin_width=1e-4, words2=words2m)
    function χexact(λ, z)
        (Es, ψs) = eigen(Hermitian(Matrix(H2 + λ*H3)))
        A = [Matrix(Aodd[ν] + λ*Aeven[ν]) for ν in 1:Nobs]
        r = [ψs' * (A[ν] * ψs[:, 1]) for ν in 1:Nobs]
        l = [ψs' * (A[ν]' * ψs[:, 1]) for ν in 1:Nobs]
        ΔE = Es .- Es[1]
        return [sum(@. conj(r[μ]) * r[ν] / (z - ΔE) - conj(l[ν]) * l[μ] / (z + ΔE)) for μ in 1:Nobs, ν in 1:Nobs]
    end
    for z in (0.7 + 0.3im, -1.1 + 0.25im)
        K = Sunny.cauchy_transform((decay, source), [z])[:, :, 1]
        (m, d) = (1:2L, 2L+1:2L+Nobs)
        G0 = inv(z*τ₃ - Diagonal(abs.(ε)))
        χ2 = w2m' * G0 * K[m, m] * G0 * w2m + K[d, m] * G0 * w2m + w2m' * G0 * K[m, d] + K[d, d]
        ref = (χexact(0.01, z) + χexact(-0.01, z) - 2χexact(0.0, z)) / (2 * 0.01^2)
        @test χ2 ≈ ref rtol=1e-3
    end

    # ---- Onsite anisotropy in closed form ----

    # For the words of at most two bosons the correction has a closed form: the
    # leading word that LSWT already holds, times ℓ = -binomial(k, 2)/2s, which is the
    # leading deviation of `rcs_factors` from unity and so vanishes in mode :dipole,
    # plus a further 1/4s on the anomalous A₂ alone. That last piece is present in
    # both modes because LSWT reads A₂ off the classical energy, in effect using a
    # boson coherent state, whose amplitude on two spin deviations exceeds a spin
    # coherent state's by 1/√(1 - 1/2s). Unlike the tests above this is exact rather
    # than asymptotic, so it pins both scalars including their signs. One Stevens
    # order at a time is required, since ℓ depends on k.
    function closed_form_errors(mode, k, s)
        sys = System(Crystal(lattice_vectors(1, 1.3, 1.7, 88, 95, 100), [[0, 0, 0]], 1),
                     [1 => Moment(; s, g=1)], mode)
        O = stevens_matrices(mode == :dipole ? s : Inf)
        set_onsite_coupling!(sys, (0.3*O[k, 0] + 0.2*O[k, 1] - 0.1*O[k, -2])/s^k, 1)
        polarize_spins!(sys, [0.3, 0.5, 0.8])
        minimize_energy!(sys; jitter=0)
        swt = SpinWaveTheory(sys; measure=nothing, regularization=0)

        # LSWT's own coefficients: A₁ of the word b†b and A₂ of the anomalous b†b†
        H = zeros(ComplexF64, 2, 2)
        Sunny.dynamical_matrix!(H, swt, zero(Sunny.Vec3))
        (; terms2, δE) = Sunny.anisotropy_correction(swt)
        coef(as) = sum(t.c for t in terms2 if t.as == as; init=0.0im)

        ℓ = mode == :dipole ? 0 : -binomial(k, 2)/2s
        E = energy_per_site(sys)
        # Errors are measured against a common scale, because in mode :dipole the
        # first two comparisons are between exact zeros
        rel(a, b) = abs(a - b) / max(abs(E), abs(H[1, 1]), abs(H[1, 2]))
        return (; energy = rel(δE, ℓ*E),
                  diagonal = rel(coef((2, 1)), ℓ*H[1, 1]),
                  anomalous = rel(2coef((2, 2)), (ℓ + 1/4s)*H[1, 2]))
    end

    for mode in (:dipole, :dipole_uncorrected), k in (2, 4, 6)
        err = closed_form_errors(mode, k, 4.5)
        @test err.energy < 1e-12
        @test err.diagonal < 1e-12
        @test err.anomalous < 1e-12
    end
end


@testitem "1/s corrections: vertex contraction (derivation)" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: anisotropic_square, bogoliubov_matrices

    # The contraction itself, against a reference that builds the whole symmetrized
    # Nambu tensor explicitly and transforms every slot at once. `vertex!` instead
    # walks the monomial list, which is why the two share nothing but the definition.
    # Its permutation and conjugation invariances need no reference tensor and are
    # checked in the default tier.
    swt = SpinWaveTheory(anisotropic_square(); measure=nothing)
    L = Sunny.nbands(swt)
    function reference_vertex(terms, qs, Ts)
        K = length(qs)
        W = zeros(ComplexF64, ntuple(_ -> 2L, K))
        perms = Sunny.slot_permutations(K)
        for t in terms, p in perms
            c = t.c / length(perms)
            for slot in 1:K
                c *= cis(2π * dot(qs[slot], t.ns[p[slot]]))
            end
            W[CartesianIndex(ntuple(slot -> t.as[p[slot]], K))] += c
        end
        U = zero(W)
        for n in CartesianIndices(U), a in CartesianIndices(W)
            U[n] += W[a] * prod(slot -> Ts[slot][a[slot], n[slot]], 1:K)
        end
        return U
    end

    q(x, y) = Sunny.Vec3(x, y, 0)
    cases = ((Sunny.cubic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.54, -0.22))),
             (Sunny.quartic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.22, 0.35), q(-0.32, -0.57))))
    for (terms, qs) in cases
        K = length(qs)
        Ts = bogoliubov_matrices(swt, qs)
        U = Sunny.vertex!(zeros(ComplexF64, ntuple(_ -> 2L, K)), terms, qs, Ts)
        @test U ≈ reference_vertex(terms, qs, Ts)
    end
end


@testitem "1/s corrections: mean fields and tadpole (derivation)" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: canted_square

    tol = 1e-6

    # ---- Two independent routes to the tadpole, canted ----

    # Thermodynamic consistency of the corrected moments, and the observable
    # correction, both on the canted structure.
    #
    # The correction to the uniform magnetization can be read from
    # `corrected_magnetic_moments`, which sums the tadpole tilt and the boson
    # depletion, or obtained as -∂/∂B of the zero-point energy. Those routes share
    # no machinery, and the second is s-independent. The identity is exact at this
    # order rather than asymptotic in 1/s, so the agreement below is limited only by
    # the quadrature and by differencing the energy; both sides being the moment
    # μ = -g𝐒 rather than the spin dipole is what makes it μ = -∂E/∂B.
    #
    # The tilt also corrects the amplitude for creating one magnon, at relative order
    # 1/s. Here `observable_corrections` displaces the boson as b → b + v within the
    # untilted frame, whereas re-expressing the transverse spin components about the
    # tilted axis rotates the observable vectors. The rotation is exact in the tilt
    # angle while the displacement is linear in it, so their difference must be smaller
    # than the correction itself by one more power of the tilt. Since the local frames
    # make `v` complex, this is sensitive to the conjugations and to the Nambu labeling.
    function canted_routes(s)
        # Differencing the zero-point energy loses accuracy, so 1e-6 is as tight as the
        # magnetization can be read anyway; the tilt comparison is a ratio of two
        # corrections from the same quadrature, so it is looser still
        B = 3s
        sys = canted_square(s, B)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        tad = Sunny.tadpole_correction(swt; tol=1e-6)
        corrected = Sunny.corrected_magnetic_moments(swt; tol=1e-6)
        δmz = sum(i -> (corrected[1, 1, 1, i] - magnetic_moments(sys)[1, 1, 1, i])[3], 1:2) / 2
        # Only the 1/s part of the energy, since `δmz` is likewise only the correction to
        # the magnetization; the classical -∂E/∂B would add the classical moment itself
        function zp(B′)
            sys′ = canted_square(s, B′)
            swt′ = SpinWaveTheory(sys′; measure=nothing)
            return Sunny.corrected_energy_per_site(swt′; tol=1e-6) - energy_per_site(sys′)
        end
        δzp = -(zp(B + 1e-4) - zp(B - 1e-4)) / 2e-4

        δc = Sunny.observable_corrections(swt; v=tad.v, tol=1e-4) -
             Sunny.observable_corrections(swt; tol=1e-4)
        (err, mag) = (0.0, 0.0)
        for i in 1:2
            σ = √(2swt.sys.κs[i])
            # The tilted axis in the local frame, where LSWT's is ẑ
            tilt = Sunny.Vec3(σ*real(tad.v[i]), σ*imag(tad.v[i]), swt.sys.κs[i])
            Rot = Sunny.rotation_between_vectors([0, 0, 1], tilt)
            for μ in 1:3
                O = Sunny.local_observable(swt.sys, swt.vacuum.frames, swt.measure, μ, i)
                O′ = Rot' * O
                δcp = ((O′[1] - O[1]) - im*(O′[2] - O[2])) / 2
                δcm = ((O′[1] - O[1]) + im*(O′[2] - O[2])) / 2
                err = max(err, abs(σ*δcp - δc[i, μ]), abs(σ*δcm - δc[2+i, μ]))
                mag = max(mag, abs(δc[i, μ]), abs(δc[2+i, μ]))
            end
        end
        return (δmz, δzp, err / mag)
    end
    rs = map(canted_routes, (1, 2))
    # Sign flipped relative to the spin dipole, since μ = -g𝐒 with g = 1 here
    @test rs[1][1] ≈ -0.0493385 atol=1e-6
    @test rs[1][2] ≈ rs[2][2] atol=1e-7
    for r in rs
        @test r[1] ≈ r[2] atol=1e-7
    end
    @test rs[1][3] < 0.005
    @test rs[1][3] / rs[2][3] ≈ 2 rtol=0.02

    # ---- The tadpole sourced by an anisotropy alone ----

    # An onsite anisotropy sources the tadpole on its own, because the quantum
    # correction δE₀ that `anisotropy_correction` makes to the classical energy depends
    # on the direction of the moment. The linear monomial that `tadpole_correction` adds
    # to its source ℓ must therefore be the gradient of δE₀. Displacing the boson by v
    # tilts the moment to (Sˣ, Sʸ) = σ(Re v, Im v) in the local frame, and a linear term
    # c b† + h.c. contributes 2 Re(c v̄), so the gradient is 2(Re c, Im c)/σ. The
    # identity relates the words of one operator in a rotated frame, so it holds
    # separately at every order in 1/s; both sides are read off at the same order,
    # making it exact rather than asymptotic.
    function anisotropic_site(n)
        cryst = Crystal(lattice_vectors(1, 1.1, 1.3, 88, 92, 95), [[0, 0, 0]], 1)
        sys = System(cryst, [1 => Moment(s=3, g=1)], :dipole_uncorrected)
        O = stevens_matrices(Inf)
        set_onsite_coupling!(sys, 0.3*O[2, 0] - 0.2*O[2, -1] + 0.02*O[4, 2], 1)
        set_dipole!(sys, n, (1, 1, 1, 1))
        return SpinWaveTheory(sys; measure=nothing)
    end
    swt = anisotropic_site([0.3, 0.5, 0.8])
    Rloc = swt.vacuum.frames[1]
    δE₀(n) = Sunny.anisotropy_correction(anisotropic_site(n)).δE
    # Tilting n by t along a transverse axis of the local frame moves that component of
    # the dipole by s t, to first order
    grad = [(δE₀(normalize(Rloc[:, 3] + 1e-5*Rloc[:, k])) - δE₀(normalize(Rloc[:, 3] - 1e-5*Rloc[:, k]))) / (2e-5 * 3)
            for k in 1:2]
    terms1 = Sunny.lswt_correction_monomials(swt, Val{1}())
    cb = only(t.c for t in terms1 if t.as == (2,))  # coefficient of b†
    @test only(t.c for t in terms1 if t.as == (1,)) ≈ conj(cb)
    @test grad ≈ 2 * [real(cb), imag(cb)] / √6 rtol=1e-6

    # ---- The mean-field operator, and its two consumers ----

    # On the canted structure, whose bonds connect distinct cells, `accum_quadratic!`
    # must agree with the phase convention of `vertex!`, which the exact-diagonalization
    # tests check independently. Contracting a quadratic monomial list at momenta
    # (𝐪, -𝐪) gives U₂ with Σ U₂[n₁,n₂] y_𝐪[n₁] y_{-𝐪}[n₂], whereas `accum_quadratic!`
    # produces H with (1/2) x†_𝐪 H x_𝐪 = (1/2) y†_𝐪 T†HT y_𝐪. Using y_𝐪[m]† = y_{-𝐪}[m̄],
    # the two agree once the latter is symmetrized over its slots. The Bogoliubov matrix
    # at -𝐪 is built from the one at 𝐪, via T_{-𝐪}[a,m] = conj(T_𝐪[ā,m̄]), because
    # `bogoliubov!` fixes the phase of each band independently and only a consistent
    # pair of matrices can be compared.
    swt = SpinWaveTheory(canted_square(1, 3); measure=nothing)
    L = Sunny.nbands(swt)
    terms2 = Sunny.hartree_fock_correction(swt; tol).terms2
    bar(m) = mod1(m + L, 2L)
    function mean_field_matrix(q)
        H = zeros(ComplexF64, 2L, 2L)
        Sunny.accum_quadratic!(H, terms2, q)
        return H
    end
    let q = Sunny.Vec3(0.23, -0.41, 0.17)
        Hq = zeros(ComplexF64, 2L, 2L)
        T = zeros(ComplexF64, 2L, 2L)
        Sunny.dynamical_matrix!(Hq, swt, q)
        Sunny.bogoliubov!(T, Hq)
        Tm = [conj(T[bar(a), bar(m)]) for a in 1:2L, m in 1:2L]
        U = Sunny.vertex!(zeros(ComplexF64, 2L, 2L), terms2, (q, -q), (T, Tm))
        A = T' * mean_field_matrix(q) * T / 2
        B = Tm' * mean_field_matrix(-q) * Tm / 2
        @test U ≈ [(A[bar(n2), n1] + B[bar(n1), n2])/2 for n1 in 1:2L, n2 in 1:2L]
    end

    # Iterating the mean fields to self-consistency must reach a fixed point that does
    # not depend on where the iteration starts. Both paths share the same
    # quadrature, so its accuracy is irrelevant here and a coarse grid keeps the
    # iterations cheap.
    ress = map((swt, Sunny.swt_with_vacuum(swt; correction=terms2))) do start
        sc = Sunny.to_self_consistent_vacuum(start; bz=(8, 8, 1))
        return (Sunny.gaussian_energy_per_site(sc), dispersion(sc, [[0.3, 0.1, 0]]))
    end
    @test ress[1][1] ≈ ress[2][1] atol=1e-8
    @test ress[1][2] ≈ ress[2][2] atol=1e-7

    # `static_self_energy` is the term linear in the correction of the shift
    # that `dispersion(swt_with_vacuum(swt; correction=terms2), q)` obtains by
    # rediagonalizing.
    scaled(λ) = [Sunny.BosonMonomial(λ*t.c, t.as, t.ns) for t in terms2]
    q = [[0.3, 0.1, 0]]
    δdisp = (dispersion(Sunny.swt_with_vacuum(swt; correction=scaled(1e-4)), q) - dispersion(Sunny.swt_with_vacuum(swt; correction=scaled(-1e-4)), q)) / 2e-4
    @test δdisp ≈ Sunny.static_self_energy(swt, q, terms2) atol=1e-6
end


@testitem "1/s corrections: quadrature (derivation)" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: triangular, square_afm

    sys = triangular()
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)
    q = [[1/2, 0, 0]]

    # The loop grid must keep both internal lines, 𝐩 and 𝐪-𝐩, off the zone centre,
    # where the cubic vertex diverges. Offsetting by half a step does that only for 𝐩:
    # 𝐪 = [0, 1/4, 0] is [1/4, 1/4, 0] in supercell RLU, so at nk = 26 the
    # reflected grid hits the zone centre exactly, and one point out of 26² then
    # dominates the integral. Nothing about this 𝐪 is singular, so grids on either side
    # of it must agree; a fixed half-step offset instead gave +1.07 - 0.23im for the
    # first band, wrong even in sign, and -2.82 for the second.
    Σgrid = [Sunny.cubic_self_energy(swt, [[0, 1/4, 0]]; η=0.02, bz=(nk, nk, 1))[:] for nk in (24, 26)]
    @test all(Σ -> isapprox(Σ, [-0.655 - 0.086im, -0.996 - 0.040im, -0.996 - 0.040im]; atol=0.012), Σgrid)

    # Binning the measure in the pair energy is a choice of quadrature, so it
    # must reproduce the frequency loop it replaces, both channels live. The
    # error is second order in the bin width relative to the imaginary part of
    # the frequencies. The reference is written out term by term from the
    # formula at the head of SelfEnergy.jl, sharing nothing but
    # `foreach_cubic_line`. The binned K is Σ̂ up to the metric and a transpose,
    # being oriented to contract as w' G w.
    #
    # The grid must be built at 𝐪, as every caller builds it: the multiplicity
    # trick of `loop_wavevectors` visits one point of each pair {𝐩, 𝐪-𝐩} and
    # doubles it, which is exact only on a grid closed under that involution.
    # Omitting 𝐪 gives a grid closed under 𝐩 ↦ -𝐩 instead, and the doubling
    # then lands on the wrong partner: the missing (b, a) contribution is what
    # symmetrizes the band block, so the two sides came out transposed and
    # disagreed by 7% at every bin width.
    terms3 = Sunny.cubic_monomials(swt)
    ωs = range(-2, 2, 6) .+ im*0.06
    k = Sunny.to_supercell_rlu(sys, q[1])
    grid = Sunny.loop_wavevectors((12, 12, 1), k)
    τ₃ = Diagonal([ones(L); -ones(L)])
    Σloop = let Σ = zeros(ComplexF64, 2L, 2L, length(ωs))
        Sunny.foreach_cubic_line(swt, terms3, k, grid, 2L) do a, _b, w, u, x, _p, _T1, _T2
            for m′ in 1:2L, m in 1:2L, iω in eachindex(ωs)
                Σ[m, m′, iω] += (a > L ? -18w : 18w) * τ₃[m, m] * conj(u[m]) * u[m′] / (ωs[iω] - x)
            end
        end
        Σ ./ grid.npts
    end
    Σbin = let (; decay, source) = Sunny.pair_measures(swt, terms3, k, grid; bin_width=0.06/16, words2=[])
        K = Sunny.cauchy_transform((decay, source), ωs)
        stack(τ₃ * transpose(K[:, :, iω]) for iω in eachindex(ωs))
    end
    @test maximum(abs, Σbin - Σloop) / maximum(abs, Σloop) < 1e-3

    # Reference for the bare two-magnon continuum K_dd, which is the whole
    # `direct` channel under `:one_loop_on_shell`. The reference is the same sum
    # over pairs of magnons, both channels, but broadening each pair
    # individually instead of binning its energy first, and contracting the
    # three observable amplitudes β by `contract` rather than through a
    # `measure`. The binning is therefore the only approximation separating the
    # two, so this doubles as the gate on `bin_width`, which
    # `corrected_channels` derives from `η` rather than exposing. Too slow for a
    # converged grid, but an identity holds grid by grid.
    function direct_unbinned(contract, swt, qs, energies, η, grid)
        (; sys) = swt
        L = Sunny.nbands(swt)
        Ncells = Sunny.nsites(sys) / Sunny.natoms(Sunny.orig_crystal(sys))
        ref = zeros(length(energies), length(qs))
        for (iq, q) in enumerate(qs)
            q_supercell = Sunny.to_supercell_rlu(sys, Sunny.Vec3(q))
            q_global = Sunny.orig_crystal(sys).recipvecs * Sunny.Vec3(q)
            words2 = Sunny.observable_pair_monomials(swt, q_global)
            lg = Sunny.loop_wavevectors(grid, q_supercell)
            Sunny.foreach_magnon_pair(swt, q_supercell, lg) do _p, w, T1, T2, ε1, ε2
                # Decay pairs at x > 0 and source pairs, of negative mass, at x < 0
                for b in 1:2L, a in 1:2L
                    (a > L) == (b > L) || continue
                    β = ntuple(μ -> Sunny.pair_amplitude(words2[μ], T1, T2, a, b), 3)
                    x = ε1[a] + ε2[b]
                    for (iω, ω) in enumerate(energies)
                        ref[iω, iq] += (a > L ? -w : w) * contract(β) * (η/π) / ((ω - x)^2 + η^2)
                    end
                end
            end
            view(ref, :, iq) ./= Ncells * lg.npts
        end
        return ref
    end

    # Orientation of the μν pair in the channels of `corrected_channels`. Every
    # combiner used elsewhere (`ssf_trace`, `ssf_perp`) puts zero weight on
    # off-diagonal `corr_pairs` and on imaginary parts, so a μν transpose —
    # which for a Hermitian S is a complex conjugation — is invisible to all of
    # them, as are the rotation and sum-rule checks, the latter being
    # homogeneous in the interference. This measure sees it, Im Sˣʸ being
    # genuinely nonzero on the 120° structure, so `direct` is fixed in sign and
    # not just in magnitude. A transpose would deviate by 2 rather than by the
    # 3e-4 of the binning.
    let
        swt5 = SpinWaveTheory(sys; measure=ssf_custom((q, ssf) -> imag(ssf[1, 2]), sys; apply_g=false))
        (η, grid) = (0.15, (8, 8, 1))
        qs = [[0.3, 0.2, 0], [1/6, 1/6, 0]]
        energies = range(0, 4, 61)
        direct = Sunny.corrected_channels(swt5, qs, :one_loop_on_shell; energies, η, bz=grid).direct
        ref = direct_unbinned(β -> imag(β[1] * conj(β[2])), swt5, qs, energies, η, grid)
        scale = maximum(abs, ref)
        @test scale > 1e-2   # the measure is not trivially zero
        @test maximum(abs, direct - ref) < 1e-3 * scale
    end

    # The sum rules elsewhere constrain the total weight of the continuum, which
    # converges exponentially here, but not its distribution in energy. This
    # gates the shape, and with it the two discretizations that produce it.
    # Measured against a 64×64 grid, the error falls as 3.5e-2, 4.5e-3, 2.8e-4
    # at 8, 16 and 32 points per direction, so the first figure below is a
    # convergence rate; the second isolates the binning of the pair energy,
    # which is an order of magnitude smaller than the grid error it is paired
    # with and so never the limiting approximation.
    let
        sys = square_afm(; field=[1.5, 0, 0])
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        qs = [[0.3, 0.2, 0]]
        (energies, η) = (range(-2, 16, 181), 0.2)
        direct(grid) = Sunny.corrected_channels(swt, qs, :one_loop_on_shell; energies, η, bz=grid).direct
        ref = direct((32, 32, 1))
        scale = maximum(abs, ref)
        @test 1e-3 < maximum(abs, direct((16, 16, 1)) - ref) / scale < 1e-2
        unbinned = direct_unbinned(β -> sum(abs2, β), swt, qs, energies, η, (32, 32, 1))
        @test maximum(abs, unbinned - ref) < 3e-4 * scale
    end
end
