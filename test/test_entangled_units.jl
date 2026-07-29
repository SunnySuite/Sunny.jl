@testitem "Reshaping entangled units" begin
    # Dimer oriented along x within a square cell, with weak inter-cell exchange.
    latvecs = [1.0 0 0; 0 1 0; 0 0 1]
    positions = [[0.0, 0, 0], [0.5, 0, 0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    # Bare system
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, 0.1, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, 0.1, Bond(2, 2, [0, 1, 0]))

    # Entangle the original cell
    esys = entangle_system(sys, [[1, 2]])
    for u in eachsite(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], u)
    end
    E0 = energy_per_site(esys)

    # Check that energy per site of a q=0 state is invariant under various
    # reshapings reshape. The [-2 0 0; …] shape flips the dimer's orientation so
    # that a unit straddles the cell boundary.
    shapes = ([-2 0 0; 0 1 0; 0 0 1],       # straddling
              [1 1 0; -1 1 0; 0 0 1],       # shear
              [3 0 0; 0 1 0; 0 0 1])        # resize ×3
    for shape in shapes
        r = reshape_supercell(esys, shape)
        @test energy_per_site(r) ≈ E0
    end

    # A system may also be reshaped prior to entangling
    for shape in shapes
        esys2 = entangle_system(reshape_supercell(sys, shape), [[1, 2]])
        for u in eachsite(esys2)
            set_coherent!(esys2, [0, 1/√2, -1/√2, 0], u)
        end
        @test energy_per_site(esys2) ≈ E0
    end

    # State set on the (reshaped) bare system prior to entangling — external
    # field and spin dipoles, including per-site overrides — is transferred to
    # the entangled system rather than reverting to the original chemical cell.
    bare = resize_supercell(sys, (3, 1, 1))
    set_field!(bare, [0, 0, 5])
    set_field_at!(bare, [0, 0, 7], (2, 1, 1, 1))
    set_dipole!(bare, [1, 0, 0], (1, 1, 1, 1))
    esys3 = entangle_system(bare, [[1, 2]])
    unc = esys3.entanglement.uncontracted
    @test unc.extfield[2, 1, 1, 1] ≈ [0, 0, 7]
    @test unc.extfield[1, 1, 1, 2] ≈ [0, 0, 5]
    @test unc.dipoles[1, 1, 1, 1] ≈ [1/2, 0, 0]
end

# TODO: Add test with magnetic unit cell larger than a single unit (i.e. not q=0
# ordering).

@testitem "Dimer Tests" begin
    import LinearAlgebra: norm, I
    import Random

    J = 1.0
    J′ = 0.1
    latvecs = [1 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0], [0.0, 0.5, 0.0]] 

    crystal = Crystal(latvecs, positions, 1; types = ["A", "B"])
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, J′, Bond(2, 2, [1, 0, 0]))

    esys = entangle_system(sys, [[1, 2]])
    bare = esys.entanglement.uncontracted
    interactions = esys.interactions_union[1]

    ### Test on-bond exchange

    onsite_operator = interactions.onsite
    S = spin_matrices(1/2)
    Sl, Su = to_product_space(S, S)
    @test onsite_operator ≈ J * (Sl' * Su)

    ### Test applied polarization

    set_dipole!(esys, [0, 1, 0], (1,1,1,1))
    @test bare.dipoles[1,1,1,1][2] ≈ 1/2
    @test bare.dipoles[1,1,1,2][2] ≈ 1/2

    ### Test external field
    set_field!(esys, [0, 0, -10])
    randomize_spins!(esys)
    minimize_energy!(esys)
    @test bare.dipoles[1][3] ≈ 1/2
    @test bare.dipoles[2][3] ≈ 1/2
    set_field!(esys, [0, 0, 0])
    minimize_energy!(esys)
    @test norm(bare.dipoles[1]) < 1e-10
    @test norm(bare.dipoles[2]) < 1e-10

    ### Test inter-bond exchange

    # Inter-unit bilinear exchange is delegated to the uncontracted system (like
    # Zeeman/Ewald), so it lives as `bilin` on the bare bonds rather than in the
    # unit's product-space tensor decomposition, which stays empty.
    @test all(pc -> isempty(pc.general.data), interactions.pair)
    for a in 1:2
        pc = only(pc for pc in bare.interactions_union[a].pair if !pc.isculled)
        @test pc.bilin ≈ J′
    end

    ### Test dispersion against analytical formula for antisymmetric channel.

    qs = [[0.2, 0.3, 0]]
    ω_ref(q, J, J′) = J*sqrt(1 + 2(J′/J) * cos(2π*q))
    ωs_analytical = ω_ref.([q[1] for q in qs], J, J′)

    set_field!(esys, [0, 0, 0])
    for unit in eachsite(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], unit)
    end
    swt = SpinWaveTheory(esys; measure=Sunny.empty_measurespec(esys), regularization=0.0)
    disp = dispersion(swt, qs)
    ωs_numerical = disp[1,:]

    @test all(both -> isapprox(both[1], both[2]; atol=1e-12), zip(ωs_analytical, ωs_numerical))

    ### Reshaped entangled system produces the correct intensities

    swt = SpinWaveTheory(esys; measure=ssf_perp(esys))
    res = intensities_bands(swt, qs)
    shape = [1 2 0; 0 1 0; 0 0 1]
    esys_sheared = reshape_supercell(esys, shape)
    swt_sheared = SpinWaveTheory(esys_sheared; measure=ssf_perp(esys_sheared))
    res_sheared = intensities_bands(swt_sheared, qs)
    @test res_sheared.disp ≈ res.disp atol=1e-11
    @test sum(res_sheared.data; dims=1) ≈ sum(res.data; dims=1) atol=1e-11

    ### Static structure factor must be zero in dipolar sector

    ssf = SampledCorrelationsStatic(esys; measure=ssf_trace(esys))
    add_sample!(ssf, esys)
    @test all(x -> isapprox(x, 0.0; atol=1e-12), ssf.parent.data)

    ### Golden test for classical dynamics ###

    esys = repeat_periodically(esys, (8, 1, 1))
    for site in eachsite(esys)
        set_coherent!(esys, [1, 0, 0, 0], site)
        i = site[1]
        set_field_at!(esys, [cos(i), sin(2i), cos(3i)], site)
    end

    energies = range(0, 2, 5)
    dt = 0.1
    measure = ssf_trace(esys)
    sc = SampledCorrelations(esys; dt, energies, measure)
    integrator = ImplicitMidpoint(dt)

    for _ in 1:100
        step!(esys, integrator)
    end
    add_sample!(sc, esys)
    res = intensities(sc, qs; energies, kT=0.05)

    @test res.data ≈ [0.6963636938867421; 4.833911994098928; 28.719089468355055; 52.61793172686838; 42.62520819431053;;]
end

@testitem "Entangled dipole-dipole SWT" begin
    import LinearAlgebra: norm, normalize, I

    # Intra-unit dipole-dipole is treated exactly by the entangled unit (folded
    # into the unit onsite), whereas an unentangled SU(N) system would factorize
    # it as ⟨S₁⟩D⟨S₂⟩. So the faithful reference is a *second* entangled system
    # in which the intra-unit dipole is instead supplied as an explicit exchange
    # matrix D. A large cell makes the inter-unit dipole tail negligible, so the
    # two entangled calculations must agree.
    units = Units(:meV, :angstrom)
    latvecs = [100.0 0 0; 0 100 0; 0 0 100]
    positions = [[0, 0, 0.0], [0, 0.005, 0.0]]
    cryst = Crystal(latvecs, positions, 1; types=["A", "B"])
    r = Sunny.global_displacement(cryst, Bond(1, 2, [0, 0, 0]))
    D = (units.vacuum_permeability/4π) * 2^2 * (I - 3*(normalize(r)*normalize(r)')) / norm(r)^3

    function build(; folded)
        sys = System(cryst, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
        set_exchange!(sys, 0.3, Bond(1, 1, [1, 0, 0]))
        set_exchange!(sys, 0.3, Bond(2, 2, [1, 0, 0]))
        set_field!(sys, [0, 0, -5])
        # The intra-unit dipole enters either through Ewald or as explicit exchange.
        folded ? set_exchange!(sys, Sunny.Mat3(D), Bond(1, 2, [0, 0, 0])) :
                 enable_dipole_dipole!(sys, units.vacuum_permeability)
        sys = entangle_system(sys, [[1, 2]])
        randomize_spins!(sys); minimize_energy!(sys)
        return sys
    end

    qs = [[0, 0, 0], [0.2, 0.1, 0], [0.4, 0, 0]]
    sys_e = build(folded=false); swt_e = SpinWaveTheory(sys_e; measure=nothing)
    sys_f = build(folded=true);  swt_f = SpinWaveTheory(sys_f; measure=nothing)

    @test energy_per_site(sys_e) ≈ energy_per_site(sys_f) atol=1e-7
    disp_e = sort(dispersion(swt_e, qs); dims=1)
    disp_f = sort(dispersion(swt_f, qs); dims=1)
    @test disp_e ≈ disp_f atol=1e-6
end

@testitem "General inter-unit coupling" begin
    import LinearAlgebra: I

    # Spin-1 dimers with a *biquadratic* inter-unit coupling. This exercises the
    # general (non-dipole) contraction path: each bare coupling is compressed via
    # an SVD in the small two-atom operator space and then embedded into the
    # units, which must reproduce the full-space operator exactly.
    latvecs = [1 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0.0], [0, 0.5, 0.0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    sys = System(crystal, [1 => Moment(s=1, g=2), 2 => Moment(s=1, g=2)], :SUN)
    set_pair_coupling!(sys, (Si, Sj) -> (Si' * Sj), Bond(1, 2, [0, 0, 0]))
    f(Si, Sj) = 0.3(Si' * Sj) + 0.1(Si' * Sj)^2
    set_pair_coupling!(sys, f, Bond(1, 1, [1, 0, 0]))
    set_pair_coupling!(sys, f, Bond(2, 2, [1, 0, 0]))

    esys = entangle_system(sys, [[1, 2]])
    pc = esys.interactions_union[1].pair[1]
    op = Sunny.bond_operator_in_tensor_space(pc, 9, 9)

    # Analytic reference: the two bare couplings (part1↔part1 and part2↔part2),
    # each embedded into the 81-dimensional two-unit product space.
    S = spin_matrices(1)
    Sl, Su = to_product_space(S, S)         # the two atoms within a 9-dim unit
    Sl1, Sl2 = to_product_space(Sl, Sl)     # part1-of-A ↔ part1-of-B
    Su1, Su2 = to_product_space(Su, Su)     # part2-of-A ↔ part2-of-B
    opref = f(Sl1, Sl2) + f(Su1, Su2)

    # The tensor decomposition retains only the non-bilinear (biquadratic)
    # residual; the bilinear part is peeled off and delegated to the uncontracted
    # system, where it is stored as `bilin` on the bare bonds. Reassembling the
    # residual with the delegated bilinear must reproduce the full operator.
    bare = esys.entanglement.uncontracted
    Jbilin = only(pc for pc in bare.interactions_union[1].pair if !pc.isculled).bilin
    bilin_op = Jbilin * ((Sl1' * Sl2) + (Su1' * Su2))
    @test op + bilin_op ≈ opref
end

@testitem "Inter-unit coupling recompression" begin
    # Two exchange couplings that land on the *same* inter-unit bond have their
    # tensor expansions concatenated during `repopulate_couplings_from_params!`.
    # Recompression must collapse them to minimal Schmidt rank while preserving
    # the bond operator (hence energy and gradient) exactly.
    latvecs = [1 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0.0], [0, 0.5, 0.0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    sys = System(crystal, [1 => Moment(s=1, g=2), 2 => Moment(s=1, g=2)], :SUN)
    set_pair_coupling!(sys, (Si, Sj) -> (Si' * Sj), Bond(1, 2, [0, 0, 0]))
    # Biquadratic couplings (bilinear is delegated to the uncontracted system and
    # never reaches the tensor decomposition). Both share part 1 of unit A, so
    # they overlap on the A-side; each contributes 5 terms and the concatenation
    # of 10 recompresses down to Schmidt rank 5.
    set_pair_coupling!(sys, (Si, Sj) -> 0.3(Si' * Sj)^2, Bond(1, 1, [1, 0, 0]))
    set_pair_coupling!(sys, (Si, Sj) -> 0.5(Si' * Sj)^2, Bond(1, 2, [1, 0, 0]))

    esys = entangle_system(sys, [[1, 2]])

    # The stored expansion on the shared bond is minimal rank (5), not the 10
    # terms that naive concatenation would produce.
    pc = only(pc for pc in esys.interactions_union[1].pair if !pc.isculled && !isempty(pc.general.data))
    @test length(pc.general.data) == 5

    # Recompression is exact: energy matches the unentangled reference in a
    # generic product state.
    randomize_spins!(sys)
    esys = entangle_system(sys, [[1, 2]])
    @test energy_per_site(esys) ≈ energy_per_site(sys) atol=1e-10
end

@testitem "Entangled units with dipole-dipole" begin
    import LinearAlgebra: norm, normalize, I
    import Random

    units = Units(:meV, :angstrom)

    # Two spin-1/2 atoms per cell. Intra-cell exchange makes them one entangled
    # unit; a weak inter-cell exchange gives a nontrivial ground state. The
    # intra-unit dipole-dipole is treated exactly (folded into the unit onsite),
    # so the faithful reference `fsys` is another entangled system that instead
    # carries the intra-unit dipole as an explicit exchange matrix D. A large cell
    # makes the inter-unit dipole tail negligible, so the two must agree.
    latvecs = [30.0 0 0; 0 30 0; 0 0 60]
    positions = [[0, 0, 0], [0, 0.005, 0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])
    r = Sunny.global_displacement(crystal, Bond(1, 2, [0, 0, 0]))
    D = (units.vacuum_permeability/4π) * 2^2 * (I - 3*(normalize(r)*normalize(r)')) / norm(r)^3

    function build(; folded)
        sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
        set_field!(sys, [0, 40, 0])
        set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
        set_exchange!(sys, 0.1, Bond(1, 1, [1, 0, 0]))
        set_exchange!(sys, 0.1, Bond(2, 2, [1, 0, 0]))
        folded ? set_exchange!(sys, Sunny.Mat3(D), Bond(1, 2, [0, 0, 0]), :dip => 1.0) :
                 enable_dipole_dipole!(sys, units.vacuum_permeability)
        esys = entangle_system(sys, [[1, 2]])
        minimize_energy!(esys)
        return esys
    end

    # The residual mismatches are the inter-unit dipole tail (Ewald includes it;
    # the folded exchange does not), which vanishes as the cell grows.
    esys = build(folded=false)
    fsys = build(folded=true)
    @test energy_per_site(esys) ≈ energy_per_site(fsys) atol=1e-5

    esys = resize_supercell(esys, (2, 1, 1))
    fsys = resize_supercell(fsys, (2, 1, 1))
    minimize_energy!(esys); minimize_energy!(fsys)
    @test energy_per_site(esys) ≈ energy_per_site(fsys) atol=1e-5

    # SpinWaveTheory dispersions and total intensities match between the two
    # entangled systems.
    qs = [[0.1, 0, 0], [0.3, 0.2, 0]]
    res1 = intensities_bands(SpinWaveTheory(fsys; measure=ssf_perp(fsys)), qs)
    res2 = intensities_bands(SpinWaveTheory(esys; measure=ssf_perp(esys)), qs)
    @test res1.disp ≈ res2.disp atol=1e-3
    @test vec(sum(res1.data; dims=1)) ≈ vec(sum(res2.data; dims=1)) atol=1e-4

    let sampler = LocalSampler(kT=0.1, propose=Sunny.propose_delta(0.1))
        @test_throws "LocalSampler does not yet support entangled units" step!(esys, sampler)
    end

    let sampler = LocalSampler(kT=0.1, propose=Sunny.propose_flip)
        @test_throws "propose_flip is not supported for general coherent states" step!(esys, sampler)
    end

    # On a state entangled *within* a unit, the intra-unit dipole energy must be
    # the exact ⟨Sᵢ D Sⱼ⟩, not the factorized ⟨Sᵢ⟩ D ⟨Sⱼ⟩ that the Ewald path
    # currently delivers. Use a large box (inter-unit/image dipole negligible) so
    # only the intra-unit dipole matters; the reference folds the direct dipole
    # tensor D into a matrix exchange, which the contraction handles exactly.
    bigcry = Crystal([100.0 0 0; 0 100 0; 0 0 100], [[0,0,0], [0,0.01,0]], 1; types=["A","B"])
    moments = [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)]
    r = Sunny.global_displacement(bigcry, Bond(1, 2, [0, 0, 0]))
    r̂ = normalize(r)
    D = (units.vacuum_permeability/4π) * 2^2 * (I - 3*(r̂ * r̂')) / norm(r)^3
    ewald_sys = System(bigcry, moments, :SUN)
    set_exchange!(ewald_sys, 1.0, Bond(1, 2, [0, 0, 0]))
    enable_dipole_dipole!(ewald_sys, units.vacuum_permeability)
    eewald = entangle_system(ewald_sys, [[1, 2]])
    ref = System(bigcry, moments, :SUN)
    set_exchange!(ref, Sunny.Mat3(1.0*I + D), Bond(1, 2, [0, 0, 0]))
    eref = entangle_system(ref, [[1, 2]])
    Z = normalize(Sunny.CVec{4}(0.1, 0.7, -0.68, 0.2))  # entangled, not a product state
    for s in eachsite(eewald)
        set_coherent!(eewald, Z, s)
        set_coherent!(eref, Z, s)
    end
    @test energy(eewald) ≈ energy(eref)

    # Only the *local* (home-image) intra-unit dipole is folded exactly into the
    # onsite; that direct tensor is box-independent, so energy per site is invariant
    # under reshaping even on a genuinely entangled (non-product) state, where the
    # exact ⟨SᵢDSⱼ⟩ contributes. A small box makes the periodic-image tail (which
    # stays in the delegated factorized Ewald sum) significant, so this would fail
    # if the whole Ewald block were folded instead.
    small = System(Crystal([5.0 0 0; 0 5 0; 0 0 5], [[0,0,0.0],[0,0.2,0]], 1; types=["A","B"]),
                   moments, :SUN)
    set_exchange!(small, 1.0, Bond(1, 2, [0, 0, 0]))
    set_exchange!(small, 0.2, Bond(1, 1, [1, 0, 0]))
    set_exchange!(small, 0.2, Bond(2, 2, [1, 0, 0]))
    enable_dipole_dipole!(small, units.vacuum_permeability)
    esmall = entangle_system(small, [[1, 2]])
    Zent = normalize(Sunny.CVec{4}(0.1, 0.7, -0.68im, 0.2))
    for s in eachsite(esmall); set_coherent!(esmall, Zent, s); end
    E0 = energy_per_site(esmall)
    for shape in ([2 0 0; 0 1 0; 0 0 1], [1 0 0; 0 2 0; 0 0 1], [3 0 0; 0 1 0; 0 0 1])
        rs = reshape_supercell(esmall, shape)
        for s in eachsite(rs); set_coherent!(rs, Zent, s); end
        @test isapprox(energy_per_site(rs), E0; atol=1e-12)
    end

    # SWT dispersion exercises the entangled-Ewald `Aq` path (the q-dependent
    # dipole matrix), which the energy checks above do not. This small box keeps
    # the periodic-image dipole tail significant, so the q≠0 dispersion is only
    # correct if the home-image direct term stripped from the q=0 matrix is also
    # removed from `Aq(q)` consistently.
    box = Crystal([5.0 0 0; 0 5 0; 0 0 6], [[0, 0, 0], [0, 0.2, 0]], 1; types=["A","B"])
    dsys = System(box, moments, :SUN)
    set_exchange!(dsys, 1.0, Bond(1, 2, [0, 0, 0]))
    set_exchange!(dsys, 0.2, Bond(1, 1, [1, 0, 0]))
    set_exchange!(dsys, 0.2, Bond(2, 2, [1, 0, 0]))
    set_field!(dsys, [0, 0, -3])
    enable_dipole_dipole!(dsys, units.vacuum_permeability; demag=Sunny.Mat3(I/3))
    # modify_exchange_with_truncated_dipole_dipole!(dsys, cutoff, units.vacuum_permeability)
    edsys = entangle_system(dsys, [[1, 2]])
    randomize_spins!(edsys); minimize_energy!(edsys)
    qs = [[0.1, 0, 0], [0.3, 0.2, 0], [0.25, 0.25, 0]]
    disp = sort(dispersion(SpinWaveTheory(edsys; measure=nothing), qs); dims=1)

    # To validate these golden values, replace enable_dipole_dipole! with
    # modify_exchange_with_truncated_dipole_dipole! per comment above. At
    # cutoff=50 the max the max deviation is 2.5e-4. At cutoff=70 it is 7.7e-5.
    golden = [4.906945628750579  4.683252377366419  4.744997203967233
              5.800316340774632  5.577887749213549  5.639863059424799
              11.597318534836047 11.597144492278936 11.597192079860365]
    @test isapprox(disp, golden; atol=1e-6)
end

@testitem "Model parameters for entangled units" begin
    latvecs = [1 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0.0], [0, 0.5, 0.0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    # A dimer (intra-unit :J) with an inter-unit exchange (:Jp) along x, so both
    # labeled parameters enter the entangled energy.
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]), :J => 1.0)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]), :Jp => 0.2)
    esys = entangle_system(sys, [[1, 2]])
    randomize_spins!(esys)

    # `get_params` reads the labels off the physical (bare) system.
    @test get_params(esys, [:J, :Jp]) == [1.0, 0.2]

    # `set_params!` regenerates the contracted couplings (changing the energy) and
    # propagates the values into the uncontracted subsystem. A round trip is exact.
    E0 = energy_per_site(esys)
    set_params!(esys, [:J, :Jp], [1.5, 0.3])
    @test get_params(esys, [:J, :Jp]) == [1.5, 0.3]
    @test get_params(esys.entanglement.uncontracted, [:J, :Jp]) == [1.5, 0.3]
    set_params!(esys, [:J, :Jp], [1.0, 0.2])
    @test energy_per_site(esys) ≈ E0

    # `make_loss_fn` evaluates on a clone, so it leaves the original `esys`
    # parameters untouched.
    loss = make_loss_fn(energy_per_site, esys, [:J, :Jp])
    loss([1.5, 0.3])
    @test get_params(esys, [:J, :Jp]) == [1.0, 0.2]

    # `set_params!` must recurse to origin's uncontracted subsystem, from which
    # a subsequent reshaping rebuilds.
    resys = reshape_supercell(esys, [1 1 0; 0 1 0; 0 0 1])
    set_params!(resys, [:J, :Jp], [1.0, 0.9])
    reresys = reshape_supercell(resys, [1 0 0; 1 1 0; 0 0 1])
    @test energy_per_site(reresys) ≈ energy_per_site(resys)
end


@testitem "Entangled Unit Intensity Scaling" begin
    latvecs = lattice_vectors(1, 1, 2, 90, 90, 90)
    positions = [[0, 0, 0], [0, 1/3, 0]]
    crystal = Crystal(latvecs, positions)

    J  = 1
    J′ = 0.12
    J2 = 0.1

    sys = System(crystal, [1 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 2, [0, -1, 0]));
    set_exchange!(sys, J2, Bond(1, 1, [1, 0, 0]));

    set_field!(sys, [0., 0., 10])

    esys = entangle_system(sys, [[1, 2]])

    for s in (sys, esys)
        minimize_energy!(s)
        swt = SpinWaveTheory(s; measure=ssf_trace(s; apply_g=false))
        qs = q_space_path(crystal, [[0, 1, 0], [1/2, 1, 0], [1, 1, 0], [0, 0, 0]], 20)
        res = intensities_bands(swt, qs)
        @test sum(res.data[:,1]) ≈ 1
    end
end


@testitem "Cell offset for straddling dimers" begin
    # A unit specified with cell offsets straddles the cell boundary from the
    # start: atom 2 at offset [0,0,0] pairs with atom 1 in the next cell along x.
    # The strong bond then lives *within* the unit (an onsite operator) and the
    # weak bond couples neighboring units.
    latvecs = [1.0 0 0; 0 1 0; 0 0 1]
    positions = [[0.0, 0, 0], [0.5, 0, 0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    J = 1.0   # Strong coupling within the straddling dimer
    J′ = 0.1  # Weak coupling between dimers

    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, J, Bond(2, 1, [1, 0, 0]))  # strong bond, atom 2 → atom 1 next cell
    set_exchange!(sys, J′, Bond(1, 2, [0, 0, 0])) # weak bond between dimers

    esys = entangle_system(sys, [[(2, [0, 0, 0]), (1, [1, 0, 0])]])

    # One unit per cell, centered at (0.5 + 1.0)/2 = 0.75
    @test length(esys.crystal.positions) == 1
    @test esys.crystal.positions[1][1] ≈ 0.75

    # The strong bond is captured in the onsite operator. Member 1 is atom 2 and
    # member 2 is atom 1, following the unit spec order.
    interactions = esys.interactions_union[1]
    S = spin_matrices(1/2)
    S2, S1 = to_product_space(S, S)
    @test interactions.onsite ≈ J * (S2' * S1)

    # The weak J′ bond couples distinct units; after symmetrization it appears as
    # a pair coupling in both [±1, 0, 0] directions.
    @test length(interactions.pair) == 2
    for pc in interactions.pair
        b = pc.bond
        @test b.i == b.j == 1
        @test abs(b.n[1]) == 1
    end
end

@testitem "Symmetry validation of entangled units" begin
    # A Cmmm cell (spacegroup 47) with a single orbit of 4 atoms forming two
    # dimers related by the mirror x → -x.
    latvecs = lattice_vectors(6, 4, 8, 90, 90, 90)
    positions = [[0.2, 0.3, 0], [0.2, 0.5, 0], [0.8, 0.3, 0], [0.8, 0.5, 0]]
    crystal = Crystal(latvecs, positions)
    @test crystal.classes == [1, 1, 1, 1]

    sys = System(crystal, [1 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))

    # Units are used verbatim (no symmetry propagation), so every unit must be
    # listed. This symmetry-consistent partition entangles.
    esys = entangle_system(sys, [[1, 2], [3, 4]])
    @test length(esys.crystal.positions) == 2

    # Listing only one of the two units leaves atoms uncovered.
    @test_throws "Atoms [3, 4] have not been assigned" entangle_system(sys, [[1, 2]])

    # A P4/n cell where naive (zero-offset) units straddle the cell so their
    # symmetry images collide. This tiles but breaks symmetry: rejected by
    # default (with a compact, symmetry-consistent suggestion), yet accepted with
    # `enforce_symmetry=false`.
    p4 = Crystal(lattice_vectors(8.333, 8.333, 5.008, 90, 90, 90),
                 [[0.1584, 0.5366, 0.6256]], "P 4/n"; choice="2")
    p4_sys = System(p4, [1 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(p4_sys, 1.0, Bond(1, 2, [0, 0, 0]))
    zero_offset = [[1, 2, 3, 4], [5, 6, 7, 8]]
    @test length(entangle_system(p4_sys, zero_offset; enforce_symmetry=false).crystal.positions) == 2

    # The default rejection reports a compacted grouping (from
    # `compactify_cell_offsets`) as a copy-pasteable suggestion.
    compact_offset = [[(1, [0, 1, 0]), (2, [0, 0, 0]), (3, [0, 0, 0]), (4, [1, 0, 0])],
                      [(5, [-1, 0, 0]), (6, [0, 0, 0]), (7, [0, 0, 0]), (8, [0, -1, 0])]]
    @test_throws "consider groupings $compact_offset" entangle_system(p4_sys, zero_offset)

    # That compacted grouping is symmetry-consistent and entangles. The two units
    # are symmetry-equivalent, so they share one class in the contracted crystal.
    ecryst = entangle_system(p4_sys, compact_offset).crystal
    @test length(ecryst.positions) == 2
    @test ecryst.classes == [1, 1]
    @test ecryst.sg.number == 85  # inherits the parent P4/n spacegroup

    # A second P4/n cell where each unit both straddles the cell and sits on the
    # fourfold axis. A symmetry-consistent grouping exists, but is not found by
    # the naive `compactify_cell_offsets` algorithm.
    p4b = Crystal(lattice_vectors(8, 8, 6, 90, 90, 90), [[0.20, 0.55, 0.30]], "P 4/n"; choice="2")
    p4b_sys = System(p4b, [1 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(p4b_sys, 1.0, Bond(1, 2, [0, 0, 0]))
    p4b_offset = [[(1, [-1, 0, 0]), (2, [0, 0, 0]), (3, [0, 0, 0]), (4, [0, -1, 0])],
                  [(5, [0, 1, 0]), (6, [0, 0, 0]), (7, [0, 0, 0]), (8, [1, 0, 0])]]
    ecryst = entangle_system(p4b_sys, p4b_offset).crystal
    @test ecryst.classes == [1, 1] && ecryst.sg.number == 85
    err = try entangle_system(p4b_sys, zero_offset); "" catch e; sprint(showerror, e) end
    @test_broken occursin("consider groupings", err)

    # A partition whose atom sets themselves break symmetry: the fourfold axis maps
    # {1,2,3,5} onto atoms that are not a unit. No offsets can fix this, so it is
    # rejected up front with a distinct message (still overridable).
    @test_throws "split by spacegroup symmetries" entangle_system(p4_sys, [[1, 2, 3, 5], [4, 6, 7, 8]])
    # With the override the units carry no symmetry, so the contracted crystal is P1
    # with every unit in its own class.
    ebroken = entangle_system(p4_sys, [[1, 2, 3, 5], [4, 6, 7, 8]]; enforce_symmetry=false).crystal
    @test length(ebroken.positions) == 2
    @test ebroken.sg.number == 1
    @test ebroken.classes == [1, 2]
end

@testitem "Offsets for different unit orientations" begin
    # A measure on an entangled system records each member's position offset from
    # its unit center. With two dimers of different orientation (one along x, one
    # along y), the offsets must point along the respective bond directions; these
    # offsets carry the intra-unit phase factors into structure-factor sums.

    latvecs = lattice_vectors(2, 2, 1, 90, 90, 90)
    positions = [
        [0.0, 0.0, 0.0],   # Atom 1
        [0.25, 0.0, 0.0],  # Atom 2 (x-dimer with atom 1)
        [0.5, 0.0, 0.0],   # Atom 3
        [0.5, 0.25, 0.0],  # Atom 4 (y-dimer with atom 3)
    ]
    cryst = Crystal(latvecs, positions, 1)
    sys = System(cryst, [i => Moment(s=1/2, g=1) for i in 1:4], :SUN; dims=(1,1,1))

    units = [
        [(1, [0,0,0]), (2, [0,0,0])],  # x-oriented dimer
        [(3, [0,0,0]), (4, [0,0,0])],  # y-oriented dimer
    ]
    esys = entangle_system(sys, units)

    # Check measure offsets reflect the different orientations
    measure = ssf_trace(esys)
    @test size(measure.offsets) == (2, 2)  # 2 unit types × 2 members

    # x-oriented dimer (unit 1): offsets should be along x
    offset_x1 = esys.crystal.latvecs * measure.offsets[1, 1]
    offset_x2 = esys.crystal.latvecs * measure.offsets[1, 2]
    @test abs(offset_x1[1]) > 0.1 && abs(offset_x1[2]) < 0.01
    @test abs(offset_x2[1]) > 0.1 && abs(offset_x2[2]) < 0.01

    # y-oriented dimer (unit 2): offsets should be along y
    offset_y1 = esys.crystal.latvecs * measure.offsets[2, 1]
    offset_y2 = esys.crystal.latvecs * measure.offsets[2, 2]
    @test abs(offset_y1[1]) < 0.01 && abs(offset_y1[2]) > 0.1
    @test abs(offset_y2[1]) < 0.01 && abs(offset_y2[2]) > 0.1
end
