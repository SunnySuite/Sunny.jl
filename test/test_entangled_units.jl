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
    latvecs = [
        1  0  0
        0  1  0
        0  0  2
    ]
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
    minimize_energy!(esys; g_tol=1e-14)
    @test norm(esys.sys_origin.dipoles[1]) < 1e-14
    @test norm(esys.sys_origin.dipoles[2]) < 1e-14

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


    # Test classical dynamics and perform golden test.
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

    # Julia 1.11 slightly changes the SVD, which leads to a different
    # decomposition of "general" pair interactions. Small numerical differences
    # can be amplified over a long dynamical trajectory.
    @static if v"1.10" <= VERSION < v"1.11"
        intensities_ref = [0.03781658031882558; 0.16151778881111165; -0.0037168273786800762; 0.003932155182202315; 0.004622186821935108;;]
    elseif v"1.11" <= VERSION
        intensities_ref = [0.06308051766790573; 0.3382120639629263; 0.01923042197493646; -0.0008116785050247188; 0.013219669220764382;;]
    end
    @test isapprox(res.data, intensities_ref)
end

# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end