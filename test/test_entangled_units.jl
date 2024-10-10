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
    qpts = q_space_path(crystal, [[0, 1, 0], [1, 1, 0]], 5)
    qs = qpts.qs

    ω_ref(q, J, J′) = J*sqrt(1 + 2(J′/J) * cos(2π*q))
    ωs_analytical = ω_ref.([q[1] for q in qs], J, J′)

    for unit in Sunny.eachunit(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], unit)
    end
    swt = SpinWaveTheory(esys; measure=Sunny.empty_measurespec(sys), regularization=0.0)
    disp = dispersion(swt, qptss)
    ωs_numerical = disp[1,:]

    @test all(both -> isapprox(both[1], both[2]; atol=1e-12), zip(ωs_analytical, ωs_numerical))


    # Test classical dynamics and perform golden test.
    esys = repeat_periodically(esys, (8, 1, 1))
    energies = range(0, 2, 10)
    dt = 0.1
    measure = ssf_trace(esys)
    sc = SampledCorrelations(esys; dt, energies, measure)
    integrator = Langevin(dt ; damping=0.4, kT=0.05)

    for _ in 1:100
        step!(esys, integrator)
    end
    add_sample!(sc, esys)
    res = intensities(sc, qpts; energies, kT=0.05)
    intensities_ref = [-0.0021805338918195185 -0.00323580926131731 -0.003475460861305824 -0.003235809261317303 -0.0021805338918195185; 0.015394114132126614 0.05504458917464571 0.006363628894456872 0.024858210899731406 0.015394114132126614; -0.062439105112575104 -0.24430343196998291 -0.1124375404098702 -0.02581168756492061 -0.062439105112575104; 0.392194235986946 2.8241763179059833 4.85842932490045 0.4259092656302079 0.392194235986946; 9.761448823359627 18.38878057796522 15.41704091322518 2.3935803890294105 9.761448823359627; 19.84319371801515 22.70779362388382 11.689296089047344 2.883234285853617 19.84319371801515; 9.747007990151367 5.296141862559427 0.8103824155070479 0.6877505519539169 9.747007990151367; -0.42870193438878 -0.9643857517663812 -0.0903901605398698 -0.22212424683216814 -0.42870193438878; 0.14524976514187082 0.2766893084879453 0.20093635005032084 0.032539493064106036 0.14524976514187082; -0.11486352541237474 -0.2600088555205098 -0.0036428642723143265 -0.26000885552050884 -0.11486352541237474]
    @test isapprox(res.data, intensities_ref)
end

# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end