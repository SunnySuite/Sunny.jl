@testitem "Crystal contraction and expansion" begin
    crystal = Sunny.diamond_crystal()

    # Specify various entangled units.
    units_all = [
        [(1, 3), (2, 4), (5, 7), (6, 8)], 
        [(1, 3, 5, 6), (2, 4, 7, 8)], 
        [(1, 3)]
    ]

    # Check that re-expansion of a contracted crystal matches original crystal
    # in terms of site-ordering and positions.
    for units in units_all
        contracted_crystal, contraction_info = Sunny.contract_crystal(crystal, units)
        expanded_crystal = Sunny.expand_crystal(contracted_crystal, contraction_info)
        @test expanded_crystal.positions ≈ crystal.positions
    end
end

# TODO: Add reshapings tests.

# TODO: Add test with magnetic unit cell larger than a single unit (i.e. not q=0
# ordering).

@testitem "Dimer Tests" begin
    import LinearAlgebra: norm
    J = 1.0
    J′ = 0.1
    latvecs = [1 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0], [0.0, 0.5, 0.0]] 

    crystal = Crystal(latvecs, positions, 1; types = ["A", "B"])
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN; seed=1)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, J′, Bond(2, 2, [1, 0, 0]))  # Needed because we broke the symmetry equivalence of the two sites

    esys = Sunny.EntangledSystem(sys, [(1, 2)])
    interactions = esys.sys.interactions_union[1]

    # Test on-bond exchange
    onsite_operator = interactions.onsite
    S = spin_matrices(1/2)
    Sl, Su = to_product_space(S, S)
    onsite_ref = J * (Sl' * Su)
    @test onsite_operator ≈ onsite_ref

    # Test external field works as expected
    set_field!(esys, [0, 0, 1])
    onsite_operator = esys.sys.interactions_union[1].onsite
    field_offset = 2*(Sl[3] + Su[3]) # 2 for g-factor
    @test onsite_operator ≈ onsite_ref + field_offset 

    # Test apparatus for setting coherent states from dipoles specification
    dipoles = [[0, 1/2, 0], [0, -1/2, 0]] # Dipoles specifying a dimer state
    cs = Sunny.coherent_state_from_dipoles(esys, dipoles, 1)
    set_coherent!(esys, cs, CartesianIndex(1,1,1,1))
    @test esys.sys_origin.dipoles[1,1,1,1][2] ≈ 1/2
    @test esys.sys_origin.dipoles[1,1,1,2][2] ≈ -1/2

    # Test external field works in action
    set_field!(esys, [0, 0, 10])
    randomize_spins!(esys)
    minimize_energy!(esys)
    @test esys.sys_origin.dipoles[1][3] ≈ -1/2
    @test esys.sys_origin.dipoles[2][3] ≈ -1/2

    set_field!(esys, [0, 0, -10])
    randomize_spins!(esys)
    minimize_energy!(esys)
    @test esys.sys_origin.dipoles[1][3] ≈ 1/2
    @test esys.sys_origin.dipoles[2][3] ≈ 1/2

    set_field!(esys, [0, 0, 0])
    randomize_spins!(esys)
    minimize_energy!(esys)
    @test norm(esys.sys_origin.dipoles[1]) < 1e-10
    @test norm(esys.sys_origin.dipoles[2]) < 1e-10

    # Test inter-bond exchange
    pc = Sunny.as_general_pair_coupling(interactions.pair[1], esys.sys)
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
    for unit in Sunny.eachunit(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], unit)
    end
    swt = SpinWaveTheory(esys; measure=Sunny.empty_measurespec(sys), regularization=0.0)
    disp = dispersion(swt, qs)
    ωs_numerical = disp[1,:]

    @test all(both -> isapprox(both[1], both[2]; atol=1e-12), zip(ωs_analytical, ωs_numerical))

    # Test static structure factor is zero (dipolar sector)
    ssf = SampledCorrelationsStatic(esys; measure=ssf_trace(esys))
    add_sample!(ssf, esys)
    @test all(x -> isapprox(x, 0.0; atol=1e-12), ssf.sc.parent.data)

    # Test classical dynamics and perform golden test.

    # Remove optimization noise for reproducibility
    randomize_spins!(esys)
    set_coherent!(esys, [0, 1/√2, -1/√2, 0], (1, 1, 1, 1))

    esys = repeat_periodically(esys, (8, 1, 1))
    energies = range(0, 2, 5)
    dt = 0.1
    measure = ssf_trace(esys)
    sc = SampledCorrelations(esys; dt, energies, measure)
    integrator = Langevin(dt; damping=0.4, kT=0.05)

    for _ in 1:100
        step!(esys, integrator)
    end
    add_sample!(sc, esys)
    res = intensities(sc, qs; energies, kT=0.05)

    # The revised SVD in Julia 1.11 leads to slight differences in the
    # decomposition of "general" pair interactions, which amplify dynamically.
    @static if v"1.10" <= VERSION < v"1.11"
        @test res.data ≈ [0.048921434630783356; 0.210889224234417; -0.008510613962195454; 0.0015281806623530139; 0.0017370496602666513;;]
    elseif v"1.11" <= VERSION
        @test res.data ≈ [0.047614832748034186; 0.24767988465544913; 0.01779072121898864; 0.0037208693822031644; 0.014790643152215031;;]
    end
end

# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end