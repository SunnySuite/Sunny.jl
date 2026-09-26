# Tests of the 1/s corrections to linear spin wave theory, in two tiers. The
# items marked `skip=true` each check one piece of the derivation against a
# private reimplementation: that is what to reach for when a result is in
# question, rather than what to run on every commit. They must still pass, so
# flip the keyword to `false` to run one, and flip it back. A literal
# `skip=true` is honored before the test item's module is created, so a skipped
# item costs nothing at all.
#
# The default tier pins end-to-end output, checks the invariants that need no
# reference implementation, and certifies the whole boson expansion against
# exact diagonalization, so that a question of correctness remains decidable
# without it.

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

    # The `aniso` switch is a `Bool` rather than a function or `nothing`, so
    # that this and `cluster_errors` each compile one specialization instead of
    # one per call site.
    function cluster(ss; mode=:dipole, aniso=false)
        sys = System(cluster_cryst, [i => Moment(s=ss[i], g=1) for i in 1:3], mode)
        for (n, (i, j)) in enumerate([(1, 2), (2, 3), (1, 3)])
            set_exchange!(sys, cluster_Js[n], Bond(i, j, [0, 0, 0]))
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
    # representation rather than in s, and needs its own cluster: two sites rather
    # than three, because the M-box Hilbert space grows as ((M+1)(M+2)/2)^Na. The
    # couplings are generic — an anisotropic bilinear exchange, a biquadratic term
    # that mode :SUN handles with no special treatment, an onsite anisotropy, and a
    # field leaving no symmetry — and every bond offset is zero, as in `cluster`.
    #
    # `pairscale` multiplies the pair coupling. The M-box classical energy is
    # M Σᵢ onsiteᵢ[N,N] + M² Σ A[N,N] B[N,N], and those differing powers make the
    # stationary reference state depend on M: the state to expand about at M boxes
    # is the one Sunny's minimizer finds for `pairscale = M`.
    const sun_cryst = Crystal(lattice_vectors(1, 1.1, 1.2, 80, 90, 100),
                              [[0, 0, 0], [0.45, 0.05, 0.1]], 1)

    function sun_cluster(; pairscale=1)
        sys = System(sun_cryst, [i => Moment(s=1, g=1) for i in 1:2], :SUN)
        set_pair_coupling!(sys, (Si, Sj) -> pairscale * (Si'*cluster_Js[1]*Sj + 0.3*(Si'*Sj)^2),
                           Bond(1, 2, [0, 0, 0]))
        set_field!(sys, cluster_B)
        S = spin_matrices(1)
        for i in 1:2
            set_onsite_coupling!(sys, 0.35*S[3]^2 + 0.2*(S[1]^2 - S[2]^2), i)
        end
        polarize_spins!(sys, cluster_B)
        minimize_energy!(sys; g_tol=1e-14, jitter=0)
        return sys
    end

    # Operators of the M-box symmetric representation of SU(N), on the basis labeled
    # by the occupations of the Nf = N-1 flavors that have left the condensate, which
    # holds the remaining M - Σn boxes. Returns that basis, the generators E_{mn} =
    # b†_m b_n, and the bosons themselves labeled as `sun_monomials` labels them,
    # flavor fastest with a ≤ L annihilating. The generators are exact; a boson
    # operator that would leave the space is truncated, so residuals may only be read
    # off blocks of low occupation.
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


@testitem "LSWT correction to classical energy" begin
    J = 1
    s = 1
    δE_afm1_ref = 0.488056/(2s) * (-2*J*s^2)

    # The results are taken from Phys. Rev. B 102, 220405(R) (2020) for the AFM1
    # phase on the FCC lattice
    function correction(mode)
        a = 1
        latvecs = lattice_vectors(a, a, a, 90, 90, 90)
        positions = [[0, 0, 0]]
        fcc = Crystal(latvecs, positions, 225)
        sys_afm1 = System(fcc, [1 => Moment(; s, g=1)], mode)
        set_exchange!(sys_afm1, J, Bond(1, 2, [0, 0, 0]))
        set_dipole!(sys_afm1, (0, 0,  1), position_to_site(sys_afm1, (0, 0, 0)))
        set_dipole!(sys_afm1, (0, 0, -1), position_to_site(sys_afm1, (1/2, 1/2, 0)))
        set_dipole!(sys_afm1, (0, 0, -1), position_to_site(sys_afm1, (1/2, 0, 1/2)))
        set_dipole!(sys_afm1, (0, 0,  1), position_to_site(sys_afm1, (0, 1/2, 1/2)))
        swt_afm1 = SpinWaveTheory(sys_afm1; measure=nothing)
        # A fixed budget rather than a `tol`, so that both modes integrate the
        # same function to the same accuracy; `tol` alone stalls in :SUN mode,
        # whose Nambu space is larger and whose norm the relative test is
        # measured against. `corrected_energy_per_site` reports an absolute
        # energy, so the classical part is subtracted off to compare against the
        # published correction. Using `sys_afm1` rather than the clone inside
        # `swt_afm1`, whose exchange has been rotated.
        δE_afm1 = Sunny.corrected_energy_per_site(swt_afm1; maxevals=2000) -
                  energy_per_site(sys_afm1)
        return isapprox(δE_afm1_ref, δE_afm1; atol=1e-3)
    end

    for mode in (:dipole, :SUN)
        @test correction(mode)
    end

    # The onsite coupling contributes a constant at this same order, which
    # `corrected_energy_per_site` now includes. It vanishes identically in
    # `:dipole` mode, where `rcs_factors` leaves the classical energy exact, so
    # only `:dipole_uncorrected` sees it; the reference above is therefore
    # unaffected, having no anisotropy at all. The two modes describe the same
    # model with different truncations of it, so their corrected energies need
    # not agree, only the correction must be present in one and absent in the
    # other.
    function easy_plane(mode)
        cryst = Sunny.square_crystal(; c=3)
        sys = System(cryst, [1 => Moment(s=2, g=1)], mode)
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> 0.3*S[3]^2 + 0.05*S[3]^4, 1)
        sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
        set_dipole!(sys, [+1, 0, 0], (1, 1, 1, 1))
        set_dipole!(sys, [-1, 0, 0], (1, 1, 1, 2))
        minimize_energy!(sys)
        return SpinWaveTheory(sys; measure=nothing)
    end
    @test iszero(Sunny.anisotropy_correction(easy_plane(:dipole)).δE)
    @test Sunny.anisotropy_correction(easy_plane(:dipole_uncorrected)).δE ≈ 0.0542857 rtol=1e-5
end


@testitem "LSWT correction to the ordered moments (s maximized)" begin
    # Test example 1: The magnetization is maximized to `s`. Reference result
    # comes from Phys. Rev. B 79, 144416 (2009) Eq. (45) for the 120° order on
    # the triangular lattice.
    J = 1
    s = 1/2
    a = 1
    δS_ref = -0.261302

    function δS_triangular(mode)
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 120)
        cryst = Crystal(latvecs, [[0, 0, 0]])
        sys = System(cryst, [1 => Moment(s=s, g=2)], mode)
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        polarize_spins!(sys, [0, 1, 0])
        sys = repeat_periodically_as_spiral(sys, (3, 3, 1); k=[2/3, -1/3, 0], axis=[0, 0, 1])
        swt = SpinWaveTheory(sys; measure=nothing)
        # Calculate first 3 digits for faster testing. At s = 1/2 there is one
        # boson per site in either mode, so the density is the dipole
        # shortening.
        δS = -Sunny.boson_density(swt; tol=1e-3)[1]

        return isapprox(δS_ref, δS, atol=1e-3)
    end

    for mode in (:dipole, :SUN)
        @test δS_triangular(mode)
    end
end


@testitem "LSWT correction to the ordered moments (s not maximized)" begin
    using LinearAlgebra
    # Test example 2: The magnetization is smaller than `s` due to easy-plane
    # single-ion anisotropy The results are derived in the Supplemental
    # Information (Note 12) of https://doi.org/10.1038/s41467-021-25591-7.
    a = b = 8.3193
    c = 5.3348
    lat_vecs = lattice_vectors(a, b, c, 90, 90, 90)
    types = ["Fe"]
    positions = [[0, 0, 0]]
    cryst = Crystal(lat_vecs, positions, 113; types)

    s = 1
    J₁  = 0.266
    J₁′ = 0.1J₁
    Δ = Δ′ = 1/3
    D = 1.42
    gab, gcc = 2.18, 1.93
    g = diagm([gab, gab, gcc])

    # Easy-plane anisotropy shrinks the classical dipole in :SUN mode, but cannot in
    # :dipole mode, where |𝐒| is pinned to s. The two therefore have different
    # classical states and different boson densities.
    for (mode, d_ref, n_ref) in [(:dipole, 1.0, 0.113839), (:SUN, 0.772087, 0.042901)]
        sys = System(cryst, [1 => Moment(; s, g)], mode; dims=(1, 1, 2), seed=0)
        set_exchange!(sys, diagm([J₁, J₁, J₁*Δ]),  Bond(1, 2, [0, 0, 0]))
        set_exchange!(sys, diagm([J₁′, J₁′, J₁′*Δ′]), Bond(1, 1, [0, 0, 1]))
        set_onsite_coupling!(sys, S -> D*S[3]^2, 1)

        randomize_spins!(sys)
        minimize_energy!(sys; maxiters=1000)
        swt = SpinWaveTheory(sys; measure=nothing)

        @test norm(sys.dipoles[1]) ≈ d_ref atol=1e-6
        @test Sunny.boson_density(swt; tol=1e-5)[1] ≈ n_ref atol=1e-6
    end
end


@testitem "1/s corrections against exact diagonalization" setup=[CorrectionModels] begin
    using LinearAlgebra, SparseArrays
    using .CorrectionModels: cluster, cluster_aniso, fock_ops, fock_quadratic, expand,
                             anisotropic_square, sun_cluster, box_ops

    # ---- The boson expansion against the exact cluster Hamiltonian ----

    # Compares the whole Holstein-Primakoff expansion, order by order in the boson
    # number, against exact diagonalization in the *spin* Hilbert space, so that the
    # truncation of the series is itself what is under test.
    function cluster_errors(ss; mode=:dipole, aniso=false)
        sys = cluster(ss; mode, aniso)
        # Regularization would otherwise leak into the quadratic coefficients
        swt = SpinWaveTheory(sys; measure=nothing, regularization=0)
        Ns = ntuple(i -> Int(2ss[i]+1), 3)
        dim = prod(Ns)
        bop = fock_ops(Ns)
        op(O, i) = reduce(kron, (k == i ? sparse(O) : sparse(1.0I, Ns[k], Ns[k]) for k in 1:3))
        S(a, i) = op(spin_matrices(ss[i])[a], i)

        # Exact cluster Hamiltonian, in the local frames that the boson expansion
        # uses and with the same Zeeman convention (+𝐁⋅𝐒). The anisotropy is rotated
        # by hand, so this is independent of the implementation, which instead works
        # from the Stevens coefficients that `swt_data!` stored.
        Rs = swt.data.local_rotations
        Hex = spzeros(ComplexF64, dim, dim)
        for i in 1:3
            Bi = Rs[i]' * (sys.gs[1, 1, 1, i]' * sys.extfield[1, 1, 1, i])
            Hex .+= sum(a -> Bi[a] * S(a, i), 1:3)
            aniso || continue
            A = Hermitian(Matrix(cluster_aniso(stevens_matrices(ss[i]), ss[i], i)))
            Hex .+= op(Matrix(Sunny.rotate_operator(A, Rs[i])), i)
        end
        # Each bond appears twice, once culled
        for int in sys.interactions_union, c in int.pair
            c.isculled && continue
            @assert iszero(c.bond.n)
            (i, j) = (c.bond.i, c.bond.j)
            J = Rs[i]' * c.bilin * Rs[j]
            Hex .+= sum(J[a, b] * S(a, i) * S(b, j) for a in 1:3, b in 1:3)
        end

        (; terms2, δE) = Sunny.anisotropy_correction(swt)
        H1 = expand(bop, Sunny.anisotropy_monomials(swt, Val{1}()), dim)
        H3 = expand(bop, Sunny.cubic_monomials(swt), dim)
        H4 = expand(bop, Sunny.quartic_monomials(swt), dim)
        # Quadratic Hamiltonian of LSWT, plus its own correction at order 1/s
        (_, H2) = fock_quadratic(swt, bop, dim; terms2)

        # Flattened Kronecker index of a boson occupation triple, the states of a
        # given total boson number, and the boson numbers of an index
        idx(ms) = 1 + sum(k -> ms[k] * prod(Ns[k+1:3]), 1:3)
        states(n) = [idx(ms) for ms in Iterators.product(0:n, 0:n, 0:n) if sum(ms) == n]
        ms(a) = ntuple(k -> mod(div(a-1, prod(Ns[k+1:3])), Ns[k]), 3)
        (vac, n1, n2) = (idx((0, 0, 0)), states(1), states(2))
        E0 = real(Hex[vac, vac])
        R = Hex - E0*I - H1 - H2 - H3 - H4

        # With no anisotropy the first four blocks below are *exact*, not merely
        # asymptotic in `s`: a matrix element between one and two bosons is purely
        # cubic, since H₄ conserves boson number modulo two, and every term of H₅ is
        # a transverse operator times the classical Sᶻ of its partner, hence
        # proportional to the transverse field that vanishes at a classical minimum.
        return (; classical = abs(E0 - (energy(sys) + 3δE)) / norm(Hex),
                  quadratic = norm(R[n1, n1]) / norm(Hex),
                  anomalous = norm(R[n2, [vac]]) / norm(Hex),
                  cubic = norm(R[n2, n1]) / norm(Hex),
                  # Unlike the blocks above this one does receive an H₆ contribution,
                  # so it is asymptotic even with no anisotropy
                  quartic = norm(R[n2, n2]) / norm(H4[n2, n2]),
                  # Largest residual over every element four bosons can reach, which
                  # is the only available measure once anisotropy makes each word
                  # truncated rather than exact. A structural zero can never be the
                  # maximum, so only the stored entries are scanned; sweeping all
                  # dim² pairs instead costs more than every Sunny call here combined.
                  reachable = maximum((abs(v) for (a, b, v) in zip(findnz(R)...)
                                       if all(ms(a) .+ ms(b) .<= 4)); init=0.0) / norm(Hex),
                  coherent = max(abs(3δE), norm(H1)) / norm(Hex),
                  hermiticity = (norm(H3 - H3') + norm(H4 - H4')) / (norm(H3) + norm(H4)))
    end

    # Exchange and field only. Two `s` triples: one certifies every exact block, the
    # second supplies the 1/s scaling of the leftover H₆ piece, which would instead
    # approach unity were the quartic term wrong at O(s⁰).
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

    # Now `cluster_aniso` on every site. Here the words are truncated rather than
    # exact — one order in 1/s is kept per word — so what is left over is the next
    # order. Requiring it to vanish faster than 1/s, the size of the retained
    # correction itself, pins every retained coefficient: an error at the order kept
    # would leave a residual of the same size as the correction, and duplicating a
    # term that LSWT already holds would leave one of order unity. Mode :dipole
    # renormalizes the stored Stevens coefficients and :dipole_uncorrected does not,
    # so both are checked.
    for mode in (:dipole, :dipole_uncorrected)
        # Every spin tuple here is `NTuple{3,Float64}`, matching the exchange-only
        # case above; mixing in an integer tuple would specialize `cluster_errors` a
        # second time, which costs more in compilation than the whole block does in
        # arithmetic.
        es = [cluster_errors(f .* (3.0, 3.0, 3.0); mode, aniso=true) for f in (1.0, 4/3)]
        @test es[1].hermiticity < 1e-12
        @test es[1].reachable < 0.02
        @test es[2].reachable < 0.5 * es[1].reachable
        # The Stevens coefficients are renormalized in mode :dipole so that the
        # classical energy function is exact in a spin coherent state, to all orders
        # in 1/s. That makes the corrections to the energy and to the linear term
        # vanish identically, and leaves the anomalous coefficient A₂ as the only
        # correction to LSWT's H₂.
        @test mode == :dipole ? es[1].coherent < 1e-12 : es[1].coherent > 1e-3
    end


    # ---- The :SUN boson expansion against the M-box Hamiltonian ----

    # The same idea in mode :SUN, whose expansion parameter is the number of boxes M
    # of the symmetric SU(N) representation. Sunny uses M = 1, where a site cannot
    # hold two bosons and the asymptotic blocks are unreachable, so `sun_monomials`
    # takes M as a test-only argument and this promotes the very same local operators
    # to the M-box representation exactly, via its generators.
    function box_errors(M)
        # The couplings to expand are the unscaled ones, `local_words` supplying the
        # powers of M itself; all the scaling is for is the reference state, so carry
        # only the coherents over from the M-box minimization.
        sysM = sun_cluster(; pairscale=M)
        sys0 = sun_cluster()
        sys0.coherents .= sysM.coherents
        swt = SpinWaveTheory(sys0; measure=nothing, regularization=0)
        # `swt.sys` is the private clone that `swt_data!` rotated into local frames
        # and absorbed the Zeeman term into; that, not `sys0`, is what
        # `sun_monomials` reads.
        sys = swt.sys
        Nf = Sunny.nflavors(swt)
        N = Nf + 1
        (; dim, E, bop, nbs) = box_ops(Nf, M, Sunny.nsites(sys))
        prom(A, i) = sum(A[m, n] * E(m, n, i) for m in 1:N, n in 1:N if !iszero(A[m, n]))

        Hex = spzeros(ComplexF64, dim, dim)
        for (i, int) in enumerate(sys.interactions_union)
            iszero(int.onsite) || (Hex .+= prom(Matrix(int.onsite), i))
            for c in int.pair
                c.isculled && continue
                @assert iszero(c.bond.n)
                for (A, B) in c.general.data
                    Hex .+= prom(Matrix(A), c.bond.i) * prom(Matrix(B), c.bond.j)
                end
            end
        end

        Hs = [expand(bop, Sunny.sun_monomials(swt, Val{K}(), M), dim) for K in 1:4]
        E0 = real(sum(t -> t.c, Sunny.sun_monomials(swt, Val{0}(), M); init=0.0+0im))
        R = Hex - E0*I - sum(Hs)
        blk(n) = findall(==(n), nbs)
        scale = norm(Hex)

        # The classical, linear, quadratic, anomalous and cubic blocks are all *exact*
        # at the M-box stationary state, not merely asymptotic: the first omitted word
        # is the five-boson -(1/8)A[m,N] b†_m n̂², whose coefficient summed over
        # interactions is the gradient of that energy. The quartic block, which the
        # six-boson word does reach, is what must fall off like 1/M.
        return (; classical = abs(E0 - real(Hex[blk(0)[1], blk(0)[1]])) / scale,
                  linear = norm(R[blk(1), blk(0)]) / scale,
                  quadratic = norm(R[blk(1), blk(1)]) / scale,
                  anomalous = norm(R[blk(2), blk(0)]) / scale,
                  cubic = norm(R[blk(2), blk(1)]) / scale,
                  quartic = norm(R[blk(2), blk(2)]) / norm(Hs[4][blk(2), blk(2)]),
                  hermiticity = (norm(Hs[3] - Hs[3]') + norm(Hs[4] - Hs[4]')) /
                                (norm(Hs[3]) + norm(Hs[4])))
    end

    # Two values of M give the scaling of the quartic residual, which would instead
    # approach a constant were the quartic term wrong at leading order.
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
    # projective space, and RCS renormalizes the Stevens coefficients so that the
    # classical energy is exactly ⟨n̂|Ĥ|n̂⟩ on that orbit. A single ion with no
    # couplings makes both statements checkable against the N × N spectrum.
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

    # The exact Hamiltonian of that ion as an N × N matrix. The Zeeman term is read
    # back from `sys.extfield` so that the unit conversion of `set_field!` needs no
    # duplicating here.
    function ion_hamiltonian(s, sys)
        S = spin_matrices(s)
        return Matrix(ion_aniso(s)) + sum(sys.extfield[1][d] * Matrix(S[d]) for d in 1:3)
    end

    for s in (1, 3/2, 2, 5/2, 3)
        # Along -ẑ, the Zeeman energy being +𝐁⋅(g𝐒), so that m = s is the minimum
        B = [0, 0, -0.3]
        # Axial, so the exact levels are already labeled by m = s, s-1, …, -s
        H = ion_hamiltonian(s, single_ion(s, :dipole, B))
        @test norm(H - Diagonal(diag(H))) < 1e-12
        lv = real(diag(H))
        gaps = sort(lv)[2:end] .- minimum(lv)

        # In :SUN the condensate is an arbitrary N-vector, so the quadratic form is
        # the exactly projected Hamiltonian: all N-1 bands come out exact, and the
        # expansion terminates. The cubic words are pure round-off and there are no
        # quartic ones at all, so every 1/M correction vanishes identically.
        swt = SpinWaveTheory(single_ion(s, :SUN, B); measure=nothing, regularization=0)
        @test sort(dispersion(swt, [[0, 0, 0]])[:]) ≈ gaps
        @test maximum(abs, [t.c for t in Sunny.sun_monomials(swt, Val{3}())]; init=0.0) < 1e-12
        @test isempty(Sunny.sun_monomials(swt, Val{4}()))
        # Exactly zero at any tolerance, there being no term left to integrate
        @test iszero(Sunny.hartree_fock_correction(swt; maxiters=1, tol=0.1).δE)
        @test iszero(Sunny.boson_density(swt; tol=0.1))

        # Mode :dipole keeps one band, and RCS makes it the exact gap to the level
        # one unit of magnetization down — *not* the smallest gap, which for an easy
        # axis is the nearly degenerate m = -s partner. Only m = s-1 survives at
        # linear order in u = sin²(θ/2), since |⟨m|θ⟩|² ∝ u^(s-m); the 2s of the
        # expansion cancels the 1/s relating curvature to frequency, so the identity
        # is s-independent and holds at every Stevens order at once.
        swt = SpinWaveTheory(single_ion(s, :dipole, B); measure=nothing, regularization=0)
        @test dispersion(swt, [[0, 0, 0]])[1] ≈ lv[2] - lv[1]
    end

    # That identity is not generic: it needs the axial symmetry that pins n̂. Adding
    # a transverse field tilts the moment, and then the classical minimum of
    # ⟨n̂|Ĥ|n̂⟩ is no longer where the quantum gradient vanishes — a coherent state
    # not being an eigenstate — leaving ⟨s-1|Ĥ|s⟩ ≠ 0 in the local frame and the
    # one-band gap merely approximate. Mode :SUN, whose variational family is the
    # whole projective space, stays exact.
    for s in (1, 2)
        B = [0.4, 0, -0.3]
        sys = single_ion(s, :SUN, B)
        swt = SpinWaveTheory(sys; measure=nothing, regularization=0)
        lv = sort(real(eigvals(Hermitian(ion_hamiltonian(s, sys)))))
        @test sort(dispersion(swt, [[0, 0, 0]])[:]) ≈ lv[2:end] .- lv[1]
    end


    # ---- Symmetries of the vertex at every slot count ----

    # Invariances that need no reference tensor, checked on the one model whose bonds
    # connect distinct cells, so that the Fourier phases the zero-offset cluster above
    # cannot reach are exercised. The contraction against an explicit reference tensor
    # is in the derivation tier.
    swt = SpinWaveTheory(anisotropic_square(); measure=nothing)
    L = Sunny.nbands(swt)
    q(x, y) = Sunny.Vec3(x, y, 0)
    cases = ((Sunny.cubic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.54, -0.22))),
             (Sunny.quartic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.22, 0.35), q(-0.32, -0.57))))
    for (terms, qs) in cases
        K = length(qs)
        Ts = Sunny.bogoliubov_matrices(swt, qs)
        buf() = zeros(ComplexF64, ntuple(_ -> 2L, K))
        U = Sunny.vertex!(buf(), terms, qs, Ts)

        # Permuting the (momentum, band) slots together leaves the vertex invariant
        for p in Sunny.slot_permutations(K)
            Up = Sunny.vertex!(buf(), terms, ntuple(t -> qs[p[t]], K), ntuple(t -> Ts[p[t]], K))
            @test permutedims(Up, invperm(collect(p))) ≈ U
        end

        # Hermiticity of the underlying operator implies U(-𝐤; n̄) = conj(U(𝐤; n)),
        # where n̄ = n + L exchanges creation and annihilation. Only the magnitudes
        # are compared, because `bogoliubov!` fixes the phase of each band
        # independently at 𝐤 and -𝐤.
        Um = Sunny.vertex(swt, terms, ntuple(t -> -qs[t], K))
        bar(n) = mod1(n + L, 2L)
        @test abs.(U) ≈ [abs(Um[CartesianIndex(map(bar, Tuple(i)))]) for i in CartesianIndices(U)]
    end
end


@testitem "1/s corrections on a lattice" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: canted_square, square_cryst

    # Default accuracy for the momentum integrals, enough for the checks below that
    # compare against a reference value: the tightest of them, Oguchi's ζ, comes out
    # 2e-8 from its reference here, thirty times inside the `atol` asserted, and
    # tightening to 1e-6 does not move any error below while it doubles the cost.
    # Checks that instead compare two corrections computed from the same quadrature,
    # or extract a ratio from them, loosen it individually.
    tol = 1e-5

    # ---- Oguchi's Z_c, on a supercell and on a multi-atom basis ----

    # For a collinear structure the cubic vertex vanishes, so the mean field is the
    # entire O(1/s) shift of the dispersion, and it is a uniform rescaling by Oguchi's
    # Z_c = 1 + ζ/2s [Prog. Theor. Phys. 13, 148 (1960)]. Both the uniformity in 𝐪 and
    # the 1/s scaling are strong tests of the four-boson coefficients; the numerical
    # values of ζ test their overall normalization. Two values of `s` suffice for the
    # scaling, and ⟨H₄⟩ being of order s⁰ makes the energy correction s-independent.
    # The two lattices differ in how the sublattices arise: a reshaped supercell for
    # the square, and the two atoms of the chemical cell for the honeycomb, whose
    # ζ = 1 - ⟨√(1-|γ_𝐪|²)⟩ with γ_𝐪 = (1 + e^{iq₁} + e^{iq₂})/3 was obtained by
    # direct quadrature.
    function neel_honeycomb(s)
        sys = System(Sunny.hexagonal_crystal(; c=3), [1 => Moment(; s, g=1)], :dipole)
        set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
        set_dipole!(sys, [0, 0, +1], (1, 1, 1, 1))
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
        @assert energy_per_site(sys) ≈ -3s^2/2
        return sys
    end

    qs = [[0.1, 0, 0], [0.3, 0.17, 0], [0.5, 0, 0], [0.13, -0.4, 0.22]]
    for (model, ζ, δE₀, ss) in ((s -> canted_square(s, 0), 0.1579474, 0.01247369, (1/2, 2)),
                                (neel_honeycomb, 0.20984170, 0.01651258, (1/2,)))
        for s in ss
            swt = SpinWaveTheory(model(s); measure=nothing)
            (; terms2, δE) = Sunny.hartree_fock_correction(swt; tol)
            Zc = Sunny.corrected_dispersion(swt, qs, terms2) ./ dispersion(swt, qs)
            @test maximum(abs, Zc .- Zc[1]) < 1e-6
            @test 2s * (Zc[1] - 1) ≈ ζ atol=1e-6
            @test δE ≈ δE₀ atol=1e-7
        end
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
    gs = Sunny.nambu_correlations(swt, [(L+1, 1, (0, 0, 0)), (1, L+1, (0, 0, 0))],
                                  Sunny.BosonMonomial{2}[]; tol=1e-4)
    @test gs[2] - gs[1] ≈ 1 atol=1e-10
    @test maximum(abs, Sunny.cubic_self_energy(swt, [[0.3, 0.1, 0]]; η=0.01, grid=(6, 6, 1))) < 1e-12
    Zcs = map(maxiters -> Sunny.corrected_dispersion(swt, [[0.3, 0.1, 0]],
                  Sunny.hartree_fock_correction(swt; maxiters, scf_tol=1e-9, tol=1e-4).terms2), (1, 100))
    @test Zcs[1] ≈ Zcs[2] atol=1e-7
    # The tadpole vanishes term by term, at 1e-34, so `tol` is irrelevant to it
    tad = Sunny.tadpole_correction(swt; tol=1e-3)
    @test all(t -> abs(t.c) < 1e-12, tad.terms2)
    @test abs(tad.δE) < 1e-12
    @test tad.dipoles ≈ [[1, 0, 0], [-1, 0, 0]] / 2

    # `corrected_magnetic_moments` both shortens each moment (|μ| = 1 → 0.886) and
    # tilts it by the tadpole (0.42° away from the field here), and is shaped and
    # signed like `magnetic_moments`, i.e. μ = -g𝐒 indexed by `Site`.
    let
        sys = canted_square(1, 3)
        swt = SpinWaveTheory(sys; measure=nothing)
        corrected = Sunny.corrected_magnetic_moments(swt; tol=1e-5)
        @test size(corrected) == size(magnetic_moments(sys)) == (1, 1, 1, 2)
        @test vec(corrected) ≈ [[-0.8240908, 0, 0.3264180], [0.8240908, 0, 0.3264180]] atol=1e-6
    end

    # Maximal-weight coherent states make this SU(N) system equivalent to the dipole
    # one, so the boson densities agree despite counting N-1 bosons per site instead of
    # one. The moments must refuse in :SUN mode, where the correction to ⟨𝐒⟩ is not
    # purely radial.
    #
    # FIXME -- once `corrected_magnetic_moments` handles :SUN mode, assert the corrected
    # moments here instead of the refusal, including the transverse part.
    let
        swt = SpinWaveTheory(canted_square(1, 3; mode=:SUN); measure=nothing)
        swt′ = SpinWaveTheory(canted_square(1, 3); measure=nothing)
        @test Sunny.boson_density(swt; tol=1e-6) ≈ Sunny.boson_density(swt′; tol=1e-6) rtol=1e-5
        @test_throws ErrorException Sunny.corrected_magnetic_moments(swt; tol=1e-3)
    end

    # ---- Modes :SUN and :dipole_uncorrected on the same s = 1/2 model ----

    # At s = 1/2 an SU(N) system has N = 2 and one boson per site, and its expansion
    # in the box number coincides term by term with the dipole expansion in 1/s. So
    # every correction must agree to machine precision, which pins the :SUN vertices
    # against the dipole ones that exact diagonalization has certified above.
    #
    # Only gauge-invariant output may be compared. `bogoliubov!` fixes the phase of
    # each band independently, and the dipole `swt_data!` puts a deliberate rotation
    # into each local frame, so the vertex tensors, the Bogoliubov matrices and the
    # dynamical matrices themselves all differ between the two modes by phases.
    #
    # The quadratures here must be *tight*, even though what is asserted is a
    # difference between two modes running the same integrand. Discretization error
    # does not cancel between them: the two modes reach the same Hamiltonian only to
    # round-off, and `hcubature` subdivides by comparing an error estimate against
    # `tol`, so a 1e-14 disagreement can send them down different refinement paths
    # and leave errors of order `tol` that do not subtract. Those land in the
    # Hartree-Fock coefficients, hence in Σstat, hence — divided by a near-degenerate
    # Dyson denominator — in `εc` and the intensities, at 1e-7 for `tol` = 1e-4.
    # Freezing Σstat across such a perturbation drops the discrepancy to 1e-11,
    # confirming the mean-field cubature as the sole amplifier. At 1e-6 every
    # refinement path is converged well past the threshold below. Pinning
    # `loop_grid` keeps the wavevector loop at the size `tol` = 0.02 would have
    # chosen, so tightening costs 5x rather than 900x.
    let
        (s, B) = (1/2, 0.6)
        qs = [[0.23, 0.11, 0], [0.4, 0.3, 0]]
        rs = map((:dipole_uncorrected, :SUN)) do mode
            sys = canted_square(s, B; mode)
            swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
            hf = Sunny.hartree_fock_correction(swt; maxiters=1, tol=1e-6)
            tad = Sunny.tadpole_correction(swt; tol=1e-6)
            return (; ε = dispersion(swt, qs),
                      E = Sunny.corrected_energy_per_site(swt; tol=1e-6),
                      n = Sunny.boson_density(swt; tol=1e-6),
                      δE = [hf.δE, tad.δE],
                      εc = Sunny.corrected_dispersion(swt, qs, [hf.terms2; tad.terms2]),
                      Σ = Sunny.cubic_self_energy(swt, qs; η=0.05, grid=(8, 8, 1)),
                      I = Sunny.corrected_intensities(swt, qs; energies=range(0, 3, 61),
                                                      η=0.1, tol=1e-6, loop_grid=(26, 26, 1)).data)
        end
        for k in keys(rs[1])
            @test maximum(abs, getfield(rs[1], k) .- getfield(rs[2], k)) < 1e-11
        end
        # Nontrivial: every quantity above is of order unity, and the intensities in
        # particular exercise the full Dyson resummation
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
            Sunny.cubic_self_energy(swt′, q; η, grid=(nk, nk, 1))[:]
        end
        # Overdamped: the width exceeds a tenth of the energy, and is η-independent
        @test all(Σ -> -imag(Σ[1]) > 0.8 * 0.05, Σs)
        @test imag(Σs[1][1]) ≈ imag(Σs[2][1]) rtol=0.1
        @test all(Σ -> -imag(Σ[3]) < 0.01, Σs)
        # Im Σ ≤ 0 in the particle block, as the Dyson resummation requires
        @test all(Σ -> all(<=(1e-12), imag.(Σ)), Σs)
    end

    # ---- Goldstone modes survive every correction ----

    # Near a protected zero mode each correction separately diverges like 1/ε, the
    # Bogoliubov matrix doing so, which makes their cancellation stringent; it is the
    # product ε × shift that must vanish. Two symmetries are gated, and they switch on
    # different corrections.
    #
    # Rotation about the field axis leaves the canted structure a zero mode at 𝐪 = 0,
    # and there all three of mean field, tadpole and cubic self-energy contribute. What
    # remains is the discretization error of the self-energy integral, falling off
    # like 1/nk.
    # What is checked is the cancellation between the static shift and the self-energy,
    # both of which the quadrature affects together, so a loose `tol` leaves the ratio
    # below unchanged in its first three digits
    swt = SpinWaveTheory(canted_square(1, 3); measure=nothing)
    t2 = [Sunny.hartree_fock_correction(swt; tol=1e-4).terms2
          Sunny.tadpole_correction(swt; tol=1e-4).terms2]
    δ = Sunny.static_self_energy(swt, [[0, 0, 0]], t2)[2]
    @test δ > 1e3
    rs = map(nk -> (δ + real(Sunny.cubic_self_energy(swt, [[0, 0, 0]]; η=0.005, grid=(nk, nk, 1))[2])) / δ, (16, 32))
    @test rs[1] ≈ 2 * rs[2] rtol=0.01
    @test rs[2] < 0.02

    # The second symmetry is broken only by an onsite anisotropy. Every term below is
    # invariant under rotation about ẑ, so the in-plane moment of this easy-plane
    # ferromagnet leaves a zero mode at 𝐪 = 0, while the cubic and linear vertices
    # vanish identically by that same symmetry. That leaves `anisotropy_correction` to
    # cancel the mean field on its own -- which it can only do because `anisotropy_words`
    # keeps exactly one order in 1/s per word. Since a perturbation of the quadratic
    # form opens a gap like its square root, expanding the anisotropy exactly instead,
    # thereby injecting a partial set of O(1/s²) terms, would gap the mode at O(1/s).
    for mode in (:dipole, :dipole_uncorrected)
        s = 2.0
        sys = System(Sunny.square_crystal(; c=3), [1 => Moment(; s, g=1)], mode)
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        O = stevens_matrices(mode == :dipole ? s : Inf)
        set_onsite_coupling!(sys, (0.5*O[2, 0] - 0.2*O[4, 0]/s^2) / (3s^2), 1)
        set_dipole!(sys, [1, 0, 0], (1, 1, 1, 1))
        sw = SpinWaveTheory(sys; measure=nothing)
        @test maximum(abs(t.c) for t in Sunny.cubic_monomials(sw)) < 1e-12
        @test maximum(abs(t.c) for t in Sunny.anisotropy_monomials(sw, Val{1}()); init=0.0) < 1e-12
        ε = dispersion(sw, [[0, 0, 0]])[1]
        mf = Sunny.hartree_fock_correction(sw; tol).terms2
        resid(t2) = Sunny.static_self_energy(sw, [[0, 0, 0]], t2)[1] * ε
        @test abs(resid([mf; Sunny.anisotropy_correction(sw).terms2])) < 1e-7 < abs(resid(mf))
    end
end

@testitem "1/s corrections to intensities" setup=[CorrectionModels] begin
    using LinearAlgebra
    using .CorrectionModels: triangular, square_afm, square_cryst

    sys = triangular()
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)

    # ---- Harmonic results, and the folding convention ----

    # The energy per site is -0.5388 J, and the magnon energy at the M point of the
    # original lattice is 2Js, the lowest of the three folded bands. Both Γ and K fold
    # onto 𝐪 = 0, where the harmonic dispersion vanishes, so all three bands are
    # gapless there. Over the whole zone the dispersion is Eq. (11) of arXiv:1306.1231,
    # written in the reciprocal lattice units of the original one-site cell; the
    # three-site cell folds 𝐪 together with 𝐪 ± 𝐊, so each wavevector gates all three
    # bands at once, and with them the folding convention.
    @test Sunny.corrected_energy_per_site(swt; tol=1e-4) ≈ -0.53881 atol=1e-4
    q = [[1/2, 0, 0]]
    @test dispersion(swt, q)[:] ≈ [√2.5, √2.5, 1] atol=1e-6
    γ(q) = (cos(2π*q[1]) + cos(2π*q[2]) + cos(2π*(q[1] + q[2]))) / 3
    ε11(q) = 3 * (1/2) * sqrt(max(0, (1 - γ(q)) * (1 + 2γ(q))))
    K = [1/3, 1/3, 0]
    @test all([[h, k, 0] for h in 0.1:0.4:0.9, k in 0.1:0.4:0.9]) do q
        isapprox(sort(dispersion(swt, [q])[:]), sort([ε11(q + n*K) for n in -1:1]); atol=1e-6)
    end

    # ---- The cubic self-energy on a lattice ----

    # Being a symmetric energy minimum, the 120° structure cannot be tilted by
    # zero-point fluctuations. Unlike the collinear case the cubic monomials are
    # individually nonzero, so their cancellation here tests their relative phases.
    tad = Sunny.tadpole_correction(swt; tol=1e-3)
    @test maximum(t -> abs(t.c), tad.terms2) < 1e-6
    @test abs(tad.δE) < 1e-12

    terms2 = [Sunny.hartree_fock_correction(swt; tol=1e-3).terms2; tad.terms2]
    δ = Sunny.static_self_energy(swt, q, terms2)[:]
    Σs = map(nk -> Sunny.cubic_self_energy(swt, q; η=0.02, grid=(nk, nk, 1))[:], (24, 48))

    # The self-energy converges like 1/nk, so a Richardson step gives the O(1/s) magnon
    # energy at the M point. It falls 27% below the harmonic value, most of that coming
    # from the cubic self-energy rather than the mean-field shift.
    εs = [1 + δ[3] + real(Σ[3]) for Σ in Σs]
    @test 2εs[2] - εs[1] ≈ 0.7316 atol=2e-3

    # Because the harmonic dispersion vanishes at Γ and K, the lower edge of the
    # two-magnon continuum touches the one-magnon branch at every wavevector, and a
    # magnon acquires a width only where the branch lies strictly inside the continuum.
    # There the width survives η → 0, as it does for the top of the band here, where
    # 2Γ/ε extrapolates to about 0.2, of the order of the maximum ~0.3 that the
    # reference reports. The M-point magnon instead sits on the boundary, and its
    # apparent width is entirely the Lorentzian tail of the regularization, falling off
    # like η.
    Σ = Sunny.cubic_self_energy(swt, q; η=0.01, grid=(48, 48, 1))[:]
    @test imag(Σ[3]) ≈ imag(Σs[2][3]) / 2 rtol=0.01
    @test imag(Σ[1]) / imag(Σs[2][1]) > 0.85

    # ---- The corrected spectral function ----

    # The corrected structure factor is a spectral function in its own right, not
    # merely one to the order worked to. Because the Dyson equation is solved in the
    # particle block, with the source channel of the self-energy frozen on shell, its
    # denominator has imaginary part at least the regulator η, so the intensity is
    # non-negative and no taller than a resolution-limited peak of the same weight; and
    # because that denominator grows as ωI, the transverse weight of each 𝐪 is exactly
    # the static weight of the corrected observables. All three properties are violated
    # at s = 1/2 by inverting the full Nambu denominator instead: near 𝐪 = [0.476, 0, 0]
    # a mirror pole is pushed up through ω = 0, giving intensities of -2.8 and +3.6
    # against a bound of 2.3, and at 𝐪 = [0.375, 0.125, 0] the weight comes out 27% low.
    # The momentum-space integrals need no great accuracy here: the identities hold for
    # any self-energy and any observable amplitudes. The options below are those that
    # `tol` selects, so the reference weight is built from the same mean fields as the
    # spectrum.
    opts = (; tol=0.01, maxevals=100_000)
    swt2 = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
    qs2 = [[0.476, 0, 0], [0.375, 0.125, 0]]

    # Static weight of the one-magnon bands, Σ_{n ≤ L} w[n, μ] conj(w[n, ν]) with
    # w = T'ũ, contracted through the measure's combiner exactly as
    # `corrected_intensities` contracts the spectral function. Passing `drop` subtracts
    # |δA|², leaving the cross term 2Re(A conj(δA)) that is the correction of relative
    # order 1/s; the sum rule below wants that, whereas the spectral weight wants the
    # full corrected amplitude. This is used again for the square lattice further down.
    function static_weights(swt, qs, δc; drop=false)
        sys = swt.sys
        cryst = Sunny.orig_crystal(sys)
        L = Sunny.nbands(swt)
        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Sunny.num_observables(swt.measure))
        δu = zero(u)
        Ncells = Sunny.nsites(sys) / Sunny.natoms(cryst)
        return map(qs) do q
            q_reshaped = Sunny.to_reshaped_rlu(sys, q)
            q_global = cryst.recipvecs * q
            Sunny.excitations!(T, H, swt, q)
            Sunny.set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
            fill!(δu, 0)
            isnothing(δc) || Sunny.accum_observable_corrections!(δu, swt, q_reshaped, q_global, δc)
            w = T' * (u + δu)
            δw = T' * δu
            corr = map(swt.measure.corr_pairs) do (μ, ν)
                c = dot(view(w, 1:L, μ), view(w, 1:L, ν))
                drop && (c -= dot(view(δw, 1:L, μ), view(δw, 1:L, ν)))
                c / Ncells
            end
            real(swt.measure.combiner(q_global, corr))
        end
    end
    δc = Sunny.observable_corrections(swt2; v=Sunny.tadpole_correction(swt2; opts...).v, opts...)
    refs = static_weights(swt2, qs2, δc)

    # A Lorentzian tail needs range rather than resolution, so the window is wide and
    # the step is a fraction of η. The residual 0.17% is the truncated tail.
    η = 0.06
    energies = range(-20, 24, 1501)
    chans = Sunny.corrected_channels(swt2, qs2; energies, η, tol=opts.tol,
                                     loop_grid=(12, 12, 1), mean_field_maxevals=opts.maxevals)
    # The two weight identities are properties of the transverse spectral function, so
    # that channel is taken on its own; the pair channel created directly by the
    # observable carries weight of its own, and their interference sums to zero only
    # over the whole zone.
    (; transverse) = chans
    @test all(≥(0), transverse + chans.cross + chans.direct)
    @test all(vec(maximum(transverse; dims=1)) .< refs ./ (π * η))
    @test vec(sum(transverse; dims=1)) * step(energies) ≈ refs rtol=5e-3

    # The interference is invisible to a trace measure integrated over the zone: 𝐒·𝐒 is
    # a scalar, so its expansion in bosons has no term linking an odd number of them to
    # an even one, and the cancellation is exact for every ω once Σ_𝐪 restores momentum
    # conservation. It is the only check here that constrains the *relative phase* of
    # the two routes to a pair, a phase that each channel alone is free of. The window
    # must reach below ω = 0, or a pair at x ≈ 0 contributes just half of its
    # Lorentzian.
    qs4 = vec([[i, j, 0] ./ 3 for i in 0:2, j in 0:2])
    chans4 = Sunny.corrected_channels(swt2, qs4; energies=range(-4, 12, 161), η=0.3,
                                      loop_grid=(3, 3, 1), mean_field_maxevals=opts.maxevals)
    @test abs(sum(chans4.cross)) < 1e-3 * sum(chans4.direct)

    # ---- Gauge invariance, where Σ̂ is genuinely off-diagonal ----

    # With three inequivalent sublattices the off-diagonal elements of Σ̂ are a third of
    # the diagonal, whereas for the 120° structure above they vanish identically: its
    # branches sit in momentum sectors that the cubic vertex cannot connect, as do those
    # of the umbrella that a field along z produces. Only the exact-diagonalization
    # cluster of the derivation tier constrains them otherwise. They are also the one
    # part of Σ̂ sensitive to the per-band phase that `bogoliubov!` fixes independently,
    # since a gauge conjugation Σ̂ → DΣ̂D† with D diagonal and unitary leaves the
    # eigenvalues of the Dyson denominator alone, and with them the poles, and leaves
    # τ₃Σ̂ Hermitian; only the contraction against the observable amplitudes T'u sees it.
    # Rotation about z is an exact symmetry of this model, fixing both diagm([1, 1, Δ])
    # and the trace structure factor while moving every local frame, so the corrected
    # intensities must be invariant under it. Breaking that gauge violates this by 5%.
    # The loop grid is coarse because an invariance holds grid by grid and needs no
    # converged integral.
    swts = map(φ -> let sys3 = triangular(; Δ=0.6, field=[0.7, 0, 1.2], φ, g=1)
                        SpinWaveTheory(sys3; measure=ssf_trace(sys3; apply_g=false))
                    end, (0.0, 0.9))
    qs3 = [[0.23, 0.11, 0], [0.37, 0.09, 0]]
    Σ3 = Sunny.cubic_self_energy(swts[1], qs3[1:1], dispersion(swts[1], qs3[1:1])[1:1];
                                 η=0.05, grid=(6, 6, 1))[1:L, 1:L, 1, 1]
    @test maximum(abs, Σ3 - Diagonal(diag(Σ3))) > 0.2 * maximum(abs, diag(Σ3))
    datas = map(swt3 -> Sunny.corrected_intensities(swt3, qs3; energies=range(0.2, 2.0, 25),
                                                    η=0.15, loop_grid=(4, 4, 1)).data, swts)
    @test datas[1] ≈ datas[2] rtol=1e-6

    # ---- Quantum sum rule on the square lattice ----

    # A field-polarized ferromagnet conserves Sᶻ. Its Bogoliubov transformation is
    # trivial, so the two-magnon channel must carry no weight at all, and nothing else
    # in the 1/s expansion is nonzero for this state either: every mean field vanishes
    # in the empty vacuum, and a collinear structure has no cubic vertex. So the
    # corrected intensities must reduce to those of linear spin wave theory, which is
    # itself exact here, the polarized state and its one-magnon excitations being
    # eigenstates.
    let
        sys = System(square_cryst, [1 => Moment(s=1, g=2)], :dipole)
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        set_field!(sys, [0, 0, 0.5])
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 1))
        @assert energy_per_site(sys) ≈ -2 - 1.0
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        qs = [[0.3, 0.2, 0]]
        (energies, η) = (range(0, 10, 101), 0.2)
        chans = Sunny.corrected_channels(swt, qs; energies, η, loop_grid=(4, 4, 1), mean_field_maxevals=1000)
        @test maximum(abs, chans.direct) < 1e-25
        @test chans.transverse + chans.cross + chans.direct ≈
              intensities(swt, qs; energies, kernel=lorentzian(fwhm=2η)).data atol=1e-12
    end

    # The fast path of `corrected_channels`. A collinear structure has no cubic vertex,
    # but its monomials cancel on merging rather than being absent, so the gate must be on
    # magnitude: the coefficients of the 26 monomials here sum to 1e-12 of the quadratic
    # Hamiltonian, against 5e-2 once a field cants the same model. Its observable
    # consequence is that `cross` comes out identically zero rather than at the 1e-24 of
    # the round-off it replaces, the rest of the result being untouched — checked to be
    # bit-identical, which the two tests below cannot see.
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

        chans = Sunny.corrected_channels(swt, [[0.3, 0.2, 0]]; energies=range(0, 12, 121),
                                         η=0.2, loop_grid=(8, 8, 1), mean_field_maxevals=1000)
        @test iszero(chans.cross) && !iszero(chans.direct)
    end

    # Weights of the three channels into which the quantum sum rule decomposes, all per
    # site and in units where a trace measure is used, so that no local frame projection
    # survives.
    #
    # Because Sᶻ = s - b†b is exact in the local frame, the two-magnon spectrum
    # saturates the longitudinal sum rule ⟨(δSᶻ)²⟩ = n(1+n) + |Δ|², where n = ⟨b†b⟩ and
    # Δ = ⟨bb⟩ follow from Wick's theorem. The transverse channel obeys an identity
    # sharper still: because the truncated S⁺ = σ(b - b†bb/4s) makes S⁻S⁺ = n̂(2s+1-n̂)
    # exact, and because completeness turns a sum of one-magnon weights over bands and
    # wavevectors into the static expectation value ⟨Â†Â⟩, the one-magnon bands must
    # carry ⟨(Sˣ)² + (Sʸ)²⟩ = s + 2s⟨n̂⟩ - ⟨n̂²⟩. Linear spin wave theory produces only
    # the first two terms; the -⟨n̂²⟩ is supplied entirely by `observable_corrections`,
    # so this fixes both the sign and the magnitude of that correction.
    # `tol` is set by the two identities below, which come out 3e-8 from exact here and
    # 4e-7 at `tol=1e-5`. The sum rule itself is limited instead by the energy integral,
    # whose 3.7e-4 does not move with `tol` at all.
    function channel_weights(sys; nq=16, tol=1e-6)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        L = Sunny.nbands(swt)

        # Onsite ⟨b†b⟩ and ⟨bb⟩, from which ⟨n̂²⟩ = ⟨n̂⟩² + ⟨n̂⟩(1+⟨n̂⟩) + |Δ|²
        ckeys = [[(L+i, i, (0, 0, 0)) for i in 1:L]; [(i, i, (0, 0, 0)) for i in 1:L]]
        gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; tol)
        ss = [swt.data.sqrtS[i]^2 for i in 1:L]
        n = real.(gs[1:L])
        n2 = @. n^2 + n * (1 + n) + abs2(gs[L+1:2L])

        # A uniform grid offset by half a step cancels the phase factors of all
        # correlations at distances below `nq`, leaving only the onsite ones above.
        qs = vec([[(a - 0.5)/nq, (b - 0.5)/nq, 0] for a in 1:nq, b in 1:nq])
        δc = Sunny.observable_corrections(swt; tol)
        harm = sum(static_weights(swt, qs, nothing)) / length(qs)
        transverse = sum(static_weights(swt, qs, δc; drop=true)) / length(qs)

        # Elastic weight of the ordered moment, shortened by zero-point fluctuations. In
        # dipole mode the shortening is the boson density `n` already gathered above.
        elastic = sum(i -> (ss[i] - n[i])^2, 1:L) / L

        # Two-magnon continuum, integrated over energy. Its wavevector average converges
        # quickly enough to use a coarser grid, which matters because each point requires
        # its own momentum-space integral. For a one-atom chemical cell this average is
        # the longitudinal weight per site.
        qs2 = vec([[(a - 0.5)/3, (b - 0.5)/3, 0] for a in 1:3, b in 1:3])
        energies = range(-2, 16, 181)
        direct = Sunny.corrected_channels(swt, qs2; energies, η=0.2, loop_grid=(12, 12, 1),
                                          mean_field_maxevals=1000).direct
        longitudinal = sum(direct) * step(energies) / length(qs2)

        return (; harm, transverse, elastic, longitudinal,
                harm_ref = sum(@. ss + 2ss*n) / L,
                transverse_ref = sum(@. ss + 2ss*n - n2) / L,
                casimir = sum(@. ss * (ss + 1)) / L)
    end

    # Both fields are `Vector{Float64}`, so `square_afm` and everything downstream of
    # it compile once rather than once per element type
    for field in ([0.0, 0, 0], [1.5, 0, 0])
        # Canting makes the onsite ⟨bb⟩ nonzero, exercising the anomalous contraction
        w = channel_weights(square_afm(; field))
        @test abs(w.harm / w.harm_ref - 1) < 1e-7
        @test abs(w.transverse / w.transverse_ref - 1) < 1e-7

        # The quantum sum rule, and the point of the whole exercise. Because 𝐒⋅𝐒 is a
        # Casimir, the elastic weight of the ordered moment, the one-magnon bands and the
        # two-magnon continuum must together carry exactly s(s+1), and each is produced by
        # a different part of this module. Linear spin wave theory saturates the rule only
        # through O(s): using its uncorrected one-magnon weights instead overshoots by
        # ⟨n̂²⟩, which is the entire O(s⁰) content of the rule, and some 3% of s(s+1) here.
        # The residual error is that of the energy integral above.
        @test abs((w.elastic + w.transverse + w.longitudinal) / w.casimir - 1) < 1e-3
        @test (w.elastic + w.harm + w.longitudinal) / w.casimir - 1 > 0.02
    end
end


@testitem "1/s corrections, end to end" setup=[CorrectionModels] begin
    using .CorrectionModels: cluster, anisotropic_square

    # Two models pinned through the public entry point, so that any change to
    # any part of the 1/s machinery moves a number here. Between them they carry
    # inequivalent sublattices, an off-diagonal Σ̂, a canted structure with a
    # nonzero tadpole, anisotropy words at two Stevens orders, unequal spins,
    # nonzero bond offsets, and a triclinic cell with no site symmetry — there
    # is no invariant being checked, only reproducibility. The two Stevens orders
    # come from `anisotropic_square`, since the k > 2 words of `cluster_aniso`
    # vanish identically at the small spins that keep the cluster cheap.
    #
    # Neither number may depend on which branches the adaptive cubature happens
    # to take, so `tol` is tightened until the result is a property of the
    # model: both agree with `tol=1e-8` to 4e-9, well inside the `rtol` asserted
    # below, for about 0.1 s each. `mean_field_maxevals` is raised so that `tol`
    # is what stops the integration. The ordered states are deterministic,
    # either from an explicit `polarize_spins!` start or hard coded, so a
    # rebuild reproduces both sums bit-for-bit.
    function pinned(swt, qs, energies, η, loop_grid)
        res = Sunny.corrected_intensities(swt, qs; energies, η, tol=1e-6, loop_grid,
                                          mean_field_maxevals=10_000_000)
        return (sum(res.data), maximum(res.data))
    end

    let sys = cluster((1.0, 3/2, 1.0); aniso=true)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        # Zero bond offsets make every vertex 𝐪-independent, so a single loop
        # wavevector integrates the self-energy exactly
        (s, m) = pinned(swt, [[0.21, -0.33, 0.12], [0.5, 0, 0]], range(0.1, 6.0, 120), 0.1, (1, 1, 1))
        @test s ≈ 136.63238024098126 rtol=1e-6
        @test m ≈ 8.328808756077763 rtol=1e-6
    end

    let sys = anisotropic_square()
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        (s, m) = pinned(swt, [[0.3, 0.1, 0], [0.37, 0.22, 0]], range(0.2, 14.0, 200), 0.2, (8, 8, 1))
        @test s ≈ 74.4824841476884 rtol=1e-6
        @test m ≈ 3.194787806986998 rtol=1e-6
    end
end


@testitem "1/s corrections: boson expansion (derivation)" setup=[CorrectionModels] skip=true begin
    using LinearAlgebra
    using .CorrectionModels: cluster, fock_ops, fock_quadratic, expand, cluster_self_energy

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
    ckeys = [(a, a′, (0, 0, 0)) for a in 1:2L for a′ in 1:2L]
    gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; tol=1e-8)
    @test gs ≈ [expect(bop(a)*bop(a′)) for (a, a′, _) in ckeys] atol=1e-7

    terms4 = Sunny.quartic_monomials(swt)
    ckeys = Sunny.correlation_keys(L, terms4)
    gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; tol=1e-8)
    (terms2, δE) = Sunny.hartree_fock_decoupling(terms4, Sunny.correlation_lookup(ckeys, gs, L))
    H4 = expand(bop, terms4, dim)
    H4mf = expand(bop, terms2, dim)

    # Wick's theorem is exact in a Gaussian state, so the constant subtracted by
    # the decoupling is precisely -⟨H₄⟩. This pins the three pairings, their sign,
    # and the fact that the constant is removed once rather than twice.
    @test -δE ≈ expect(H4) atol=1e-7

    # Defining property of the decoupling: the mean-field operator reproduces the
    # response of H₄ to every quadratic perturbation.
    Qs = [bop(a)*bop(a′) for a in 1:2L, a′ in 1:2L]
    @test [expect(H4*Q - Q*H4) for Q in Qs] ≈ [expect(H4mf*Q - Q*H4mf) for Q in Qs] atol=1e-6

    # The analogous decoupling of H₃, which Wick-contracts down to the linear
    # operator that tadpole relaxation must cancel. Only linear perturbations test
    # anything here, since a Gaussian state gives ⟨[H₃, Q]⟩ = 0 for quadratic Q; a
    # commutator of two linear operators is a c-number, so this pins the three
    # cubic pairings exactly rather than to integration accuracy.
    terms3 = Sunny.cubic_monomials(swt)
    ckeys = Sunny.correlation_keys(L, terms3)
    gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; tol=1e-8)
    ℓ = Sunny.tadpole_vector(terms3, Sunny.correlation_lookup(ckeys, gs, L), L, 1e-6)
    H3raw = expand(bop, terms3, dim)
    H3mf = sum(a -> ℓ[a] * bop(a), 1:2L)
    @test [expect(H3raw*bop(a) - bop(a)*H3raw) for a in 1:2L] ≈
          [expect(H3mf*bop(a) - bop(a)*H3mf) for a in 1:2L] atol=1e-8
    # The linear part must be removed from the perturbation below: normal ordering
    # H₃ in the quasi-particle basis leaves that piece behind, and it is the
    # tadpole correction rather than a loop.
    H3 = H3raw - H3mf

    τ₃ = Diagonal([ones(L); -ones(L)])
    T0 = zeros(ComplexF64, 2L, 2L)
    ε = copy(Sunny.bogoliubov!(T0, H))

    # The whole Nambu matrix of the cubic self-energy, against the exact retarded
    # Green function. This fixes the off-diagonal elements, which mix bands and
    # correct spectral weights, and the anomalous blocks, which admix three magnons
    # into the ground state. Nothing cheaper constrains the off-diagonal ones: a
    # collinear structure has Σ̂ = 0 outright, and the branches of a spiral or an
    # umbrella live in momentum sectors that the cubic vertex cannot mix, so Σ̂ is
    # diagonal there too. They are also the one part of Σ̂ sensitive to the per-band
    # phase that `bogoliubov!` fixes independently; leaving it free corrupts them by
    # 70%, two orders of magnitude above this tolerance.
    ωs = [0.5 + 0.3im, -1.1 + 0.25im]
    Σed = cluster_self_energy(H2, H3, bop, T0, ε, 0.06, ωs)
    Σm = Sunny.cubic_self_energy(swt, [[0, 0, 0]], ωs; η=1e-10, grid=(1, 1, 1))
    @test [Σm[:, :, iω, 1] for iω in eachindex(ωs)] ≈ Σed atol=3e-4
    @test maximum(abs, Σm[:, :, 1, 1] - Diagonal(diag(Σm[:, :, 1, 1]))) >
          0.2 * maximum(abs, diag(Σm[:, :, 1, 1]))

    # Below the three-magnon threshold, unitarity requires τ₃Σ̂ to be Hermitian,
    # which is what makes the Dyson equation preserve spectral weight. The
    # broadening η is what breaks it, by an amount η ∂Σ/∂ω.
    Σh = Sunny.cubic_self_energy(swt, [[0, 0, 0]], [0.5]; η=1e-10, grid=(1, 1, 1))[:, :, 1, 1]
    @test τ₃ * Σh ≈ (τ₃ * Σh)' atol=1e-9

    # Second-order perturbation theory in H₃, evaluated exactly in the truncated
    # Fock space, against the on-shell form of `cubic_self_energy`, which takes no
    # frequencies. The vacuum shift constrains the source channel alone, and the
    # level shifts the diagonal of Σ̂ at the on-shell frequency. Three magnons are
    # created and destroyed in the former, and the unrestricted sum over their
    # bands supplies 3! orderings, cancelling one of the two factors of 3! that
    # relate the symmetrized vertex to Γ₂. Levels are identified by proximity to
    # the harmonic band energy, the Fock spectrum interleaving two-magnon states
    # among the one-magnon ones, so a band permutation would be caught here. The
    # residual is the O(h²) error of the differencing.
    function levels(λ)
        E = eigen(Hermitian(Matrix(H2 + λ*H3))).values
        ΔE = E .- E[1]
        return [E[1]; [ΔE[argmin(abs.(ΔE .- ε[n]))] for n in 1:L]]
    end
    shifts = ((levels(0.01) + levels(-0.01))/2 - levels(0)) / 0.01^2
    Σ = Sunny.cubic_self_energy(swt, [[0, 0, 0]]; η=1e-10, grid=(1, 1, 1))
    @test real(vec(Σ)) ≈ shifts[2:L+1] atol=1e-3
    @test maximum(abs, imag(Σ)) < 1e-8
    U3 = Sunny.vertex(swt, terms3, ntuple(_ -> zero(Sunny.Vec3), 3))
    @test shifts[1] ≈ -6 * sum(abs2(U3[L+n1, L+n2, L+n3]) / (ε[n1] + ε[n2] + ε[n3])
                           for n1 in 1:L, n2 in 1:L, n3 in 1:L) atol=1e-6

    # Both amplitudes for the observable to create a pair of magnons, against exact
    # matrix elements in the same Fock space. This is what pins their relative phase,
    # on which the interference of `corrected_channels` depends and nothing else does:
    # each squared amplitude is invariant under a phase on either one separately, so
    # the two diagonal blocks of the measure — the self-energy and the `direct` channel
    # — are blind to an error here. So are the sum rule and the rotation-invariance
    # checks elsewhere, both being homogeneous in the interference. A sign error was
    # found this way, `pair_amplitude` having returned the amplitude of +b†b where
    # Sᶻ = s - b†b carries a minus.
    #
    # The comparison needs a two-magnon state and the operator that reaches it. For a
    # normalized pair |ab⟩ = y_a†y_b†|0⟩/√(1+δ_ab), each of Sunny's amplitudes is
    # conj(⟨ab|·|·⟩) up to the combinatorial factor √(2/(1+δ_ab)) by which the
    # symmetrized state differs from Sunny's sum over ordered pairs. The readout is
    # the full S^{μν} tensor, so that every component is reachable and nothing is
    # protected by the scalar sum rule that a trace measure enjoys.
    swt2m = SpinWaveTheory(sys; measure=ssf_custom((q, ssf) -> ssf, sys; apply_g=false))
    Nobs = Sunny.num_observables(swt2m.measure)
    u2 = zeros(ComplexF64, 2L, Nobs)
    q2m = [0.23, -0.41, 0.17]
    qr2 = Sunny.to_reshaped_rlu(sys, q2m)
    qg2 = Sunny.orig_crystal(sys).recipvecs * q2m
    Sunny.set_swt_observable_vectors!(u2, swt2m, qr2, qg2)
    words2m = Sunny.observable_pair_words(swt2m, qr2, qg2)

    # Observable A_ν(q) as a linear form in the Nambu vector. The index swap is forced,
    # not chosen: `intensities_bands` forms Avec[μ] = dot(u[:,μ], T[:,n]) for the left
    # amplitude ⟨0|A†|n⟩, so A = Σ_a u[ā] x_a. The `nambu_correlations` block above
    # fixes the state, so nothing here is free, and the transverse amplitude must come
    # out equal to Sunny's own w = T†u.
    Aodd = [sum(a -> u2[Sunny.nambu_conj(a, L), ν] * bop(a), 1:2L) for ν in 1:Nobs]
    # The even part, whose words `observable_pair_words` supplies. They carry the
    # conjugated Fourier phase, being the amplitude to create a pair, so they are
    # conjugated back to describe the same operator as `Aodd`.
    Aeven = [sum(w -> conj(w.c) * bop(w.as[1]) * bop(w.as[2]), words2m[ν]) for ν in 1:Nobs]

    Y = [sum(a -> (τ₃ * T0' * τ₃)[m, a] * bop(a), 1:2L) for m in 1:2L]
    ψ1 = [Y[m]' * ψ for m in 1:L]
    @test [dot(ψ1[n], Aodd[ν] * ψ) for n in 1:L, ν in 1:Nobs] ≈ (T0' * u2)[1:L, :] atol=1e-7

    U3m = zeros(ComplexF64, 2L, 2L, 2L)
    T3m = [conj(T0[Sunny.nambu_conj(a, L), Sunny.nambu_conj(b, L)]) for a in 1:2L, b in 1:2L]
    Sunny.vertex!(U3m, terms3, ntuple(_ -> zero(Sunny.Vec3), 3), (T0, T0, T3m), similar(U3m))
    w2m = T0' * u2

    # Each amplitude separately, then both at once, which is the quantity the
    # interference actually depends on. A phase convention shared by both would cancel
    # from each separately and survive in the coherent sum, so that sum is compared
    # against first-order perturbation theory for ⟨ab|A_ν|0⟩: the even part reaches the
    # pair directly, the odd part through a one-magnon intermediate.
    #
    # The same matrix elements, accumulated into the measure that `corrected_channels`
    # actually consumes, certify it in its own storage convention rather than as
    # rederived here. Linear bin splitting preserves the zeroth and first moments of
    # the measure exactly, so those two moments pin it without replicating the binning.
    # The 12 block is the interference, which no published calculation constrains.
    ρ2m = Matrix{ComplexF64}[]
    onshell2m = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
    Sunny.accum_pair_measure!(ρ2m, swt2m, terms3, qr2, Sunny.LoopGrid([zero(Sunny.Vec3)], [1.0], 1);
                              source_freqs=onshell2m, bin_width=1e-3, words2=words2m)
    (M0, M1) = (sum(ρ2m), sum(((i, r),) -> (i - 1) * 1e-3 * r, enumerate(ρ2m)))
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


@testitem "1/s corrections: vertex contraction (derivation)" setup=[CorrectionModels] skip=true begin
    using LinearAlgebra
    using .CorrectionModels: anisotropic_square

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
        Ts = Sunny.bogoliubov_matrices(swt, qs)
        U = Sunny.vertex!(zeros(ComplexF64, ntuple(_ -> 2L, K)), terms, qs, Ts)
        @test U ≈ reference_vertex(terms, qs, Ts)
    end
end


@testitem "1/s corrections: mean fields and tadpole (derivation)" setup=[CorrectionModels] skip=true begin
    using LinearAlgebra
    using .CorrectionModels: canted_square

    tol = 1e-6

    # ---- Two independent routes to the tadpole, canted ----

    # Thermodynamic consistency, and the observable correction, both on the canted
    # structure and both of the form "two routes agree to leading order, so their
    # difference falls off like 1/s".
    #
    # The correction to the uniform magnetization can be read from
    # `corrected_magnetic_moments`, which combines the tadpole tilt with the reduction of
    # the moment magnitude, or obtained as ∂/∂B of the zero-point energy. Those routes
    # share no machinery, and the second is s-independent. The tilt is norm-preserving,
    # so composing it with the radial shortening is unambiguous; the alternative of
    # subtracting δs along the original axis is a non-geodesic split that differs at
    # 1.6% for s = 1, falling like 1/s.
    #
    # The tilt also corrects the amplitude for creating one magnon, at relative order
    # 1/s. Here `observable_corrections` displaces the boson as b → b + v within the
    # untilted frame, whereas re-expressing the transverse spin components about the
    # tilted axis rotates the observable vectors. The rotation is exact in the tilt
    # angle while the displacement is linear in it, so their difference must be smaller
    # than the correction itself by one more power of the tilt. Since the local frames
    # make `v` complex, this is sensitive to the conjugations and to the Nambu labeling.
    function canted_routes(s)
        # Differencing the zero-point energy loses accuracy, so 1e-5 is as tight as the
        # magnetization can be read anyway; the tilt comparison is a ratio of two
        # corrections from the same quadrature, so it is looser still
        B = 3s
        sys = canted_square(s, B)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        tad = Sunny.tadpole_correction(swt; tol=1e-5)
        # `corrected_magnetic_moments` applies the tilt and the shortening together, which
        # is what makes it comparable to ∂/∂B below; either alone would be an incomplete
        # 1/s. Both sides are the moment μ = -g𝐒, not the spin dipole, so that the
        # comparison is the thermodynamic identity μ = -∂E/∂B.
        corrected = Sunny.corrected_magnetic_moments(swt; tol=1e-5)
        δmz = sum(i -> (corrected[1, 1, 1, i] - magnetic_moments(sys)[1, 1, 1, i])[3], 1:2) / 2
        # Only the 1/s part of the energy, since `δmz` is likewise only the correction to
        # the magnetization; the classical -∂E/∂B would add the classical moment itself
        function zp(B′)
            sys′ = canted_square(s, B′)
            swt′ = SpinWaveTheory(sys′; measure=nothing)
            return Sunny.corrected_energy_per_site(swt′; tol=1e-5) - energy_per_site(sys′)
        end
        δzp = -(zp(B + 1e-4) - zp(B - 1e-4)) / 2e-4

        δc = Sunny.observable_corrections(swt; v=tad.v, tol=1e-4) -
             Sunny.observable_corrections(swt; tol=1e-4)
        (err, mag) = (0.0, 0.0)
        for i in 1:2
            R = swt.data.local_rotations[i]
            σ = √2 * swt.data.sqrtS[i]
            Rot = Sunny.rotation_between_vectors([0, 0, 1], R' * tad.dipoles[i])
            for μ in 1:3
                (O, O′) = (swt.data.observables[μ, i], Rot' * swt.data.observables[μ, i])
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
    @test rs[1][1] ≈ -0.0485824 atol=1e-6
    @test rs[1][2] ≈ rs[2][2] atol=1e-8
    @test rs[1][1] - rs[1][2] ≈ 2 * (rs[2][1] - rs[2][2]) rtol=0.03
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
    Rloc = swt.data.local_rotations[1]
    δE₀(n) = Sunny.anisotropy_correction(anisotropic_site(n)).δE
    # Tilting n by t along a transverse axis of the local frame moves that component of
    # the dipole by s t, to first order
    grad = [(δE₀(normalize(Rloc[:, 3] + 1e-5*Rloc[:, k])) - δE₀(normalize(Rloc[:, 3] - 1e-5*Rloc[:, k]))) / (2e-5 * 3)
            for k in 1:2]
    terms1 = Sunny.anisotropy_monomials(swt, Val{1}())
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
    # not depend on how strongly the iteration is damped. Both paths share the same
    # quadrature, so its accuracy is irrelevant here and a loose `tol` keeps the
    # twenty-odd iterations cheap.
    ress = map((0.0, 0.5)) do damping
        r = Sunny.hartree_fock_correction(swt; maxiters=100, scf_tol=1e-9, damping, tol=1e-3)
        return (r.δE, Sunny.corrected_dispersion(swt, [[0.3, 0.1, 0]], r.terms2))
    end
    @test ress[1][1] ≈ ress[2][1] atol=1e-8
    @test ress[1][2] ≈ ress[2][2] atol=1e-7

    # `static_self_energy` is the term linear in the correction of the shift that
    # `corrected_dispersion` obtains by rediagonalizing.
    scaled(λ) = [Sunny.BosonMonomial(λ*t.c, t.as, t.ns) for t in terms2]
    q = [[0.3, 0.1, 0]]
    δdisp = (Sunny.corrected_dispersion(swt, q, scaled(1e-4)) - Sunny.corrected_dispersion(swt, q, scaled(-1e-4))) / 2e-4
    @test δdisp ≈ Sunny.static_self_energy(swt, q, terms2) atol=1e-6
end


@testitem "1/s corrections: quadrature (derivation)" setup=[CorrectionModels] skip=true begin
    using LinearAlgebra
    using .CorrectionModels: triangular, square_afm

    sys = triangular()
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)
    q = [[1/2, 0, 0]]

    # The loop grid must keep both internal lines, 𝐩 and 𝐪-𝐩, off the zone centre,
    # where the cubic vertex diverges. Offsetting by half a step does that only for 𝐩:
    # 𝐪 = [0, 1/4, 0] is [1/4, 1/4, 0] in the reshaped cell, so at nk = 26 the
    # reflected grid hits the zone centre exactly, and one point out of 26² then
    # dominates the integral. Nothing about this 𝐪 is singular, so grids on either side
    # of it must agree; a fixed half-step offset instead gave +1.07 - 0.23im for the
    # first band, wrong even in sign, and -2.82 for the second.
    Σgrid = [Sunny.cubic_self_energy(swt, [[0, 1/4, 0]]; η=0.02, grid=(nk, nk, 1))[:] for nk in (24, 26)]
    @test all(Σ -> isapprox(Σ, [-0.655 - 0.086im, -0.996 - 0.040im, -0.996 - 0.040im]; atol=0.012), Σgrid)

    # Binning the decay measure in the pair energy is a choice of quadrature, not a
    # change of interface, so it must reproduce the frequency loop it replaces. The
    # error is second order in the bin width relative to the regulator Γ, which here is
    # carried by the imaginary part of the frequencies. The reference is written out
    # term by term, straight from the formula at the head of SelfEnergy.jl, so that
    # nothing but `foreach_cubic_line` is shared with the implementation under test.
    #
    # The grid must be built at 𝐪, as every caller builds it: the multiplicity trick
    # of `loop_wavevectors` visits one point of each pair {𝐩, 𝐪-𝐩} and doubles it,
    # which is exact only on a grid closed under that involution. Omitting 𝐪 gives a
    # grid closed under 𝐩 ↦ -𝐩 instead, and the doubling then lands on the wrong
    # partner: the missing (b, a) contribution is what symmetrizes the band block, so
    # the two sides came out transposed and disagreed by 7% at every bin width.
    terms3 = Sunny.cubic_monomials(swt)
    ε = dispersion(swt, q)[:]
    onshell = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
    ωs = range(0, 2, 6) .+ im*0.06
    k = Sunny.to_reshaped_rlu(sys, q[1])
    grid = Sunny.loop_wavevectors((12, 12, 1), k)
    Σloop = let Σ = zeros(ComplexF64, L, L, length(ωs))
        Sunny.foreach_cubic_line(swt, terms3, k, grid, L) do a, _b, w, u, x, _T1, _T2
            for m′ in 1:L, m in 1:L
                R = 18w * conj(u[m]) * u[m′]
                # The source channel a > L is frozen at its on-shell frequency, hence
                # is ω-independent and enters with the opposite sign
                for iω in eachindex(ωs)
                    Σ[m, m′, iω] += a > L ? -R / (onshell[m, m′] - x) : R / (ωs[iω] - x)
                end
            end
        end
        Σ ./ grid.npts
    end
    Σbin = let ρ = Matrix{ComplexF64}[]
        Σsrc = Sunny.accum_pair_measure!(ρ, swt, terms3, k, grid; source_freqs=onshell, bin_width=0.06/16)
        Sunny.pair_self_energy!(zeros(ComplexF64, L, L, length(ωs)), ρ, Σsrc, ωs, 0.06/16)
    end
    @test maximum(abs, Σbin - Σloop) / maximum(abs, Σloop) < 1e-3

    # Reference for the `direct` channel of `corrected_channels`: the same sum over pairs
    # of magnons, but broadening each pair individually instead of binning its energy
    # first, and contracting the three observable amplitudes β by `contract` rather than
    # through a `measure`. The binning is therefore the only approximation separating the
    # two, so this doubles as the gate on `bin_width`, which `corrected_channels` derives
    # from `η` rather than exposing. Too slow for a converged grid, but an identity holds
    # grid by grid.
    function direct_unbinned(contract, swt, qs, energies, η, grid)
        (; sys) = swt
        L = Sunny.nbands(swt)
        Ncells = Sunny.nsites(sys) / Sunny.natoms(Sunny.orig_crystal(sys))
        ref = zeros(length(energies), length(qs))
        for (iq, q) in enumerate(qs)
            q_reshaped = Sunny.to_reshaped_rlu(sys, Sunny.Vec3(q))
            q_global = Sunny.orig_crystal(sys).recipvecs * Sunny.Vec3(q)
            words2 = Sunny.observable_pair_words(swt, q_reshaped, q_global)
            lg = Sunny.loop_wavevectors(grid, q_reshaped)
            Sunny.foreach_magnon_pair(swt, q_reshaped, lg) do _p, w, T1, T2, ε1, ε2
                for b in 1:L, a in 1:L
                    β = ntuple(μ -> Sunny.pair_amplitude(words2[μ], T1, T2, a, b), 3)
                    x = ε1[a] + ε2[b]
                    for (iω, ω) in enumerate(energies)
                        ref[iω, iq] += w * contract(β) * (η/π) / ((ω - x)^2 + η^2)
                    end
                end
            end
            view(ref, :, iq) ./= Ncells * lg.npts
        end
        return ref
    end

    # Orientation of the μν pair in the channels of `corrected_channels`. Every combiner
    # used elsewhere (`ssf_trace`, `ssf_perp`) puts zero weight on off-diagonal
    # `corr_pairs` and on imaginary parts, so a μν transpose — which for a Hermitian S is
    # a complex conjugation — is invisible to all of them, as are the rotation and
    # sum-rule checks, the latter being homogeneous in the interference. This measure
    # sees it, Im Sˣʸ being genuinely nonzero on the 120° structure, so `direct` is fixed
    # in sign and not just in magnitude. A transpose would deviate by 2 rather than by
    # the 3e-4 of the binning.
    let
        swt5 = SpinWaveTheory(sys; measure=ssf_custom((q, ssf) -> imag(ssf[1, 2]), sys; apply_g=false))
        (η, grid) = (0.15, (8, 8, 1))
        qs = [[0.3, 0.2, 0], [1/6, 1/6, 0]]
        energies = range(0, 4, 61)
        direct = Sunny.corrected_channels(swt5, qs; energies, η, loop_grid=grid).direct
        ref = direct_unbinned(β -> imag(β[1] * conj(β[2])), swt5, qs, energies, η, grid)
        scale = maximum(abs, ref)
        @test scale > 1e-2   # the measure is not trivially zero
        @test maximum(abs, direct - ref) < 1e-3 * scale
    end

    # The sum rules elsewhere constrain the total weight of the continuum, which
    # converges exponentially here, but not its distribution in energy. This gates the
    # shape, and with it the two discretizations that produce it. Measured against a
    # 64×64 grid, the error falls as 3.5e-2, 4.5e-3, 2.8e-4 at 8, 16 and 32 points per
    # direction, so the first figure below is a convergence rate; the second isolates
    # the binning of the pair energy, which is an order of magnitude smaller than the
    # grid error it is paired with and so never the limiting approximation.
    let
        sys = square_afm(; field=[1.5, 0, 0])
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        qs = [[0.3, 0.2, 0]]
        (energies, η) = (range(-2, 16, 181), 0.2)
        direct(grid) = Sunny.corrected_channels(swt, qs; energies, η, loop_grid=grid,
                                                mean_field_maxevals=1000).direct
        ref = direct((32, 32, 1))
        scale = maximum(abs, ref)
        @test 1e-3 < maximum(abs, direct((16, 16, 1)) - ref) / scale < 1e-2
        unbinned = direct_unbinned(β -> sum(abs2, β), swt, qs, energies, η, (32, 32, 1))
        @test maximum(abs, unbinned - ref) < 3e-4 * scale
    end
end
