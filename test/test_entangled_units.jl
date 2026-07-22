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

    # Cannot reshape then entangle
    sys_sheared = resize_supercell(sys, (2, 1, 1))
    @test_throws "Entangle a single-cell system first, then reshape" entangle_system(sys_sheared, [(1, 2)])

    # Entangle the original cell
    esys = entangle_system(sys, [(1, 2)])
    for u in eachsite(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], u)
    end
    E0 = energy_per_site(esys)

    # Check that energy per site of a q=0 state is invariant under various
    # reshapings reshape. The [-2 0 0; …] shape flips the dimer's orientation so
    # that a unit straddles the cell boundary.
    for shape in ([-2 0 0; 0 1 0; 0 0 1],       # straddling
                  [1 1 0; -1 1 0; 0 0 1],       # shear
                  [3 0 0; 0 1 0; 0 0 1])        # resize ×3
        r = reshape_supercell(esys, shape)
        @test energy_per_site(r) ≈ E0
    end
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

    esys = entangle_system(sys, [(1, 2)])
    bare = esys.entanglement.bare_system
    interactions = esys.interactions_union[1]

    # Test on-bond exchange
    onsite_operator = interactions.onsite
    S = spin_matrices(1/2)
    Sl, Su = to_product_space(S, S)
    onsite_ref = J * (Sl' * Su)
    @test onsite_operator ≈ onsite_ref

    # Test apparatus for setting coherent states from dipoles specification
    dipoles = [[0, 1/2, 0], [0, -1/2, 0]] # Dipoles specifying a dimer state
    cs = Sunny.coherent_state_from_dipoles(esys, dipoles, 1)
    set_coherent!(esys, cs, CartesianIndex(1,1,1,1))
    @test bare.dipoles[1,1,1,1][2] ≈ 1/2
    @test bare.dipoles[1,1,1,2][2] ≈ -1/2

    # Test external field. The Zeeman term is not folded into the onsite
    # operator; it lives in `extfield` (per unit) and couples to the unit's total
    # moment via the per-unit total-moment operator. The onsite operator is
    # unchanged.
    set_field!(esys, [0, 0, -10])
    @test esys.interactions_union[1].onsite ≈ onsite_ref
    @test esys.extfield[1] ≈ [0, 0, -10]
    @test bare.extfield[1, 1, 1, 1] ≈ [0, 0, -10]
    @test bare.extfield[1, 1, 1, 2] ≈ [0, 0, -10]

    # Test external field works in action
    randomize_spins!(esys)
    minimize_energy!(esys)
    @test bare.dipoles[1][3] ≈ 1/2
    @test bare.dipoles[2][3] ≈ 1/2
    set_field!(esys, [0, 0, 0])
    minimize_energy!(esys)
    @test norm(bare.dipoles[1]) < 1e-10
    @test norm(bare.dipoles[2]) < 1e-10

    # Test inter-bond exchange
    pc = Sunny.as_general_pair_coupling(interactions.pair[1], esys)
    Sl1, Sl2 = to_product_space(Sl, Sl)
    Su1, Su2 = to_product_space(Su, Su)
    bond_operator = zeros(ComplexF64, 16, 16)
    for (A, B) in pc.general.data
        bond_operator .+= kron(A, B)
    end
    bond_ref = J′*((Sl2' * Sl1) .+ (Su2' * Su1))
    @test bond_operator ≈ bond_ref

    # Test dispersion against analytical formula for antisymmetric channel.
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

    # Reshaped entangled system produces the correct intensities
    swt = SpinWaveTheory(esys; measure=ssf_perp(esys))
    res = intensities_bands(swt, qs)
    shape = [1 2 0; 0 1 0; 0 0 1]
    esys_sheared = reshape_supercell(esys, shape)
    swt_sheared = SpinWaveTheory(esys_sheared; measure=ssf_perp(esys_sheared))
    res_sheared = intensities_bands(swt_sheared, qs)
    @test res_sheared.disp ≈ res.disp atol=1e-11
    @test sum(res_sheared.data; dims=1) ≈ sum(res.data; dims=1) atol=1e-11

    # Test static structure factor is zero (dipolar sector)
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

    esys = entangle_system(sys, [(1, 2)])
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

    # A dimer (intra-unit :J) with an inter-unit exchange (:Jp) along x, so that
    # the entangled system's dispersion is pinned by both labeled parameters.
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]), :J => 1.0)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]), :Jp => 0.2)
    esys = entangle_system(sys, [(1, 2)])
    for u in eachsite(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], u)
    end

    qs = [[0.1, 0, 0], [0.25, 0, 0], [0.4, 0, 0]]
    dispersion_at(s) = let
        swt = SpinWaveTheory(s; measure=Sunny.empty_measurespec(s), regularization=0.0)
        dispersion(swt, qs)[1, :]
    end

    # `get_param` reads the labels off the physical (bare) system.
    @test get_params(esys, [:J, :Jp]) == [1.0, 0.2]

    # `set_params!` mutably regenerates the contracted couplings, changing the
    # dispersion. A round-trip back to the original values reproduces it exactly.
    disp0 = dispersion_at(esys)
    set_params!(esys, [:J, :Jp], [1.5, 0.3])
    @test get_params(esys, [:J, :Jp]) == [1.5, 0.3]
    disp1 = dispersion_at(esys)
    @test !(disp1 ≈ disp0)
    set_params!(esys, [:J, :Jp], [1.0, 0.2])
    @test dispersion_at(esys) ≈ disp0

    # `make_loss_fn` operates on a clone: it sees the trial parameters, and the
    # original `esys` parameters are left untouched. The loss vanishes exactly at
    # the target parameters that generated `disp1`.
    loss = make_loss_fn(esys, [:J, :Jp]) do s
        return squared_error(disp1, dispersion_at(s))
    end
    @test loss([1.5, 0.3]) < 1e-20
    @test loss([1.0, 0.2]) > 1e-6
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

    esys = entangle_system(sys, [(1, 2)])

    for s in (sys, esys)
        minimize_energy!(s)
        swt = SpinWaveTheory(s; measure=ssf_trace(s; apply_g=false))
        qs = q_space_path(crystal, [[0, 1, 0], [1/2, 1, 0], [1, 1, 0], [0, 0, 0]], 20)
        res = intensities_bands(swt, qs)
        @test sum(res.data[:,1]) ≈ 1
    end
end


@testitem "Cell offset for straddling dimers" begin
    import LinearAlgebra: norm

    # This test verifies that cell offsets work correctly when specifying units
    # that straddle cell boundaries from the start.
    #
    # Geometry: 2 atoms per cell at positions [0, 0, 0] and [0.5, 0, 0]
    # We create dimers where atom 2 from one cell is paired with atom 1 from the
    # next cell, specified as: [(2, [0,0,0]), (1, [1,0,0])]

    latvecs = [1.0 0 0; 0 1 0; 0 0 1]
    positions = [[0.0, 0, 0], [0.5, 0, 0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    J = 1.0   # Strong coupling within the straddling dimer
    J′ = 0.1  # Weak coupling between dimers

    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)

    # The strong bond connects atom 2 to atom 1 in the next cell
    set_exchange!(sys, J, Bond(2, 1, [1, 0, 0]))
    # The weak bond connects the dimers
    set_exchange!(sys, J′, Bond(1, 2, [0, 0, 0]))

    # Key: specify the straddling dimer using cell offsets
    # This unit contains atom 2 at offset [0,0,0] and atom 1 at offset [1,0,0]
    esys = entangle_system(sys, [[(2, [0, 0, 0]), (1, [1, 0, 0])]])

    # Verify the contracted crystal has the expected structure
    @test length(esys.crystal.positions) == 1  # One unit per cell

    # The unit center should be at (0.5 + 1.0)/2 = 0.75
    @test esys.crystal.positions[1][1] ≈ 0.75

    # Verify the strong J coupling was captured in the onsite operator
    interactions = esys.interactions_union[1]
    S = spin_matrices(1/2)
    # Note: member 1 is atom 2, member 2 is atom 1 (order from the unit spec)
    S2, S1 = to_product_space(S, S)
    onsite_ref = J * (S2' * S1)
    @test interactions.onsite ≈ onsite_ref

    # The weak J′ bond connects different units, so it appears as a pair coupling
    # It gets symmetrized, so we have both [±1, 0, 0] directions
    @test length(interactions.pair) == 2
    for pc in interactions.pair
        b = pc.bond
        @test b.i == b.j == 1
        @test abs(b.n[1]) == 1  # Cell offset along x
    end

    # Find the ground state
    randomize_spins!(esys)
    minimize_energy!(esys)
    E0 = energy_per_site(esys)

    # Test intensities_bands calculation with full structure factor
    qs = [[0.1, 0, 0], [0.2, 0, 0]]

    swt = SpinWaveTheory(esys; measure=ssf_perp(esys))
    res = intensities_bands(swt, qs)

    @test all(isfinite, res.disp)
    @test all(isfinite, res.data)

    # Test that reshaping preserves the physics
    esys_reshaped = reshape_supercell(esys, [2 0 0; 0 1 0; 0 0 1])
    @test energy_per_site(esys_reshaped) ≈ E0 atol=1e-12

    # Verify intensities_bands still works after reshaping (different number of bands due to larger system)
    swt_reshaped = SpinWaveTheory(esys_reshaped; measure=ssf_perp(esys_reshaped))
    res_reshaped = intensities_bands(swt_reshaped, qs)

    @test all(isfinite, res_reshaped.disp)
    @test all(isfinite, res_reshaped.data)
end


# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end
