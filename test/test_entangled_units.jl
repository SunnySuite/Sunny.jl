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

    pc = interactions.pair[1]
    op = Sunny.bond_operator_in_tensor_space(pc, 4, 4)
    Sl1, Sl2 = to_product_space(Sl, Sl)
    Su1, Su2 = to_product_space(Su, Su)
    op_ref = J′*((Sl2' * Sl1) .+ (Su2' * Su1))
    @test op ≈ op_ref

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
    # Compare systems with and without entanglement. The result should match if
    # there's no intra-unit exchange.
    latvecs = [1.0 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0.0], [0, 0.5, 0.0]]
    cryst = Crystal(latvecs, positions, 1; types=["A", "B"])

    function build(; entangled)
        sys = System(cryst, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
        set_exchange!(sys, 0.3, Bond(1, 1, [1, 0, 0]))
        set_exchange!(sys, 0.3, Bond(2, 2, [1, 0, 0]))
        enable_dipole_dipole!(sys, Units(:meV, :angstrom).vacuum_permeability)
        set_field!(sys, [0, 0, -5])
        entangled && (sys = entangle_system(sys, [[1, 2]]))
        randomize_spins!(sys); minimize_energy!(sys)
        return sys
    end

    qs = [[0, 0, 0], [0.2, 0.1, 0], [0.5, 0, 0]]
    sys_b = build(entangled=false); swt_b = SpinWaveTheory(sys_b; measure=nothing)
    sys_e = build(entangled=true);  swt_e = SpinWaveTheory(sys_e; measure=nothing)

    @test energy_per_site(sys_b) ≈ energy_per_site(sys_e)
    # Entangled dispersion carries an extra high-energy product-space band; its
    # two lowest bands must match the bare two-band dispersion.
    disp_b = sort(dispersion(swt_b, qs); dims=1)
    disp_e = sort(dispersion(swt_e, qs); dims=1)
    @test disp_e[1:2, :] ≈ disp_b atol=1e-10
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

    @test op ≈ opref
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
    # Both couplings share part 1 of unit A, so they overlap on the A-side and
    # recompress from 6 concatenated terms down to Schmidt rank 3.
    set_pair_coupling!(sys, (Si, Sj) -> 0.3(Si' * Sj), Bond(1, 1, [1, 0, 0]))
    set_pair_coupling!(sys, (Si, Sj) -> 0.5(Si' * Sj), Bond(1, 2, [1, 0, 0]))

    esys = entangle_system(sys, [[1, 2]])

    # The stored expansion on the shared bond is minimal rank (3), not the 6
    # terms that naive concatenation would produce.
    pc = only(pc for pc in esys.interactions_union[1].pair if !pc.isculled && !isempty(pc.general.data))
    @test length(pc.general.data) == 3

    # Recompression is exact: energy matches the unentangled reference.
    randomize_spins!(sys)
    minimize_energy!(sys)
    esys = entangle_system(sys, [[1, 2]])
    Sunny.transfer_uncontracted_state!(esys, sys)
    @test energy_per_site(esys) ≈ energy_per_site(sys) atol=1e-10
end

@testitem "Entangled units with dipole-dipole" begin
    import LinearAlgebra: norm, I
    import Random

    units = Units(:meV, :angstrom)

    # Two spin-1/2 atoms per cell, well separated so the moments are physical
    # point dipoles. Intra-cell exchange makes them one entangled unit; a weak
    # inter-cell exchange gives a nontrivial ground state.
    latvecs = [3.0 0 0; 0 3 0; 0 0 6]
    positions = [[0, 0, 0], [0, 0.15, 0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
    set_field!(sys, [0, 40, 0])
    set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, 0.1, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, 0.1, Bond(2, 2, [1, 0, 0]))
    enable_dipole_dipole!(sys, units.vacuum_permeability)
    minimize_energy!(sys)

    esys = entangle_system(sys, [[1, 2]])
    @test energy_per_site(esys) ≈ energy_per_site(sys)

    sys = resize_supercell(sys, (2, 1, 1))
    esys = resize_supercell(esys, (2, 1, 1))
    @test energy_per_site(esys) ≈ energy_per_site(sys)

    # SpinWaveTheory calculations should be consistent too
    qs = [[0.1, 0, 0], [0.3, 0.2, 0]]
    res1 = intensities_bands(SpinWaveTheory(sys; measure=ssf_perp(sys)), qs)
    res2 = intensities_bands(SpinWaveTheory(esys; measure=ssf_perp(esys)), qs)

    # The two highest energy bands describe intra-unit singlet-triplet
    # excitations and do not carry intensity.
    @test all(<(1e-12), abs.(res2.data[1:2, :]))

    # The remaining bands should match between entangled and non-entangled
    # calculations.
    @test res1.disp ≈ res2.disp[3:end, :]
    @test vec(sum(res1.data; dims=1)) ≈ vec(sum(res2.data; dims=1))

    # LocalSampler still does not support long-range dipole-dipole for entangled units.
    let sampler = LocalSampler(kT=0.1, propose=Sunny.propose_flip)
        @test_throws "does not yet support long-range dipole-dipole" step!(esys, sampler)
    end
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
    for u in eachsite(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], u)
    end

    # `get_params` reads the labels off the physical (bare) system.
    @test get_params(esys, [:J, :Jp]) == [1.0, 0.2]

    # `set_params!` regenerates the contracted couplings (changing the energy) and
    # propagates the values into the uncontracted subsystem. A round trip is exact.
    E0 = energy_per_site(esys)
    set_params!(esys, [:J, :Jp], [1.5, 0.3])
    @test get_params(esys, [:J, :Jp]) == [1.5, 0.3]
    @test get_params(esys.entanglement.uncontracted, [:J, :Jp]) == [1.5, 0.3]
    @test !(energy_per_site(esys) ≈ E0)
    set_params!(esys, [:J, :Jp], [1.0, 0.2])
    @test energy_per_site(esys) ≈ E0

    # `make_loss_fn` evaluates on a clone, so it leaves the original `esys`
    # parameters untouched.
    loss = make_loss_fn(energy_per_site, esys, [:J, :Jp])
    loss([1.5, 0.3])
    @test get_params(esys, [:J, :Jp]) == [1.0, 0.2]
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

# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end
