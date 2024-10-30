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
    esys.sys_origin.dipoles[1][3] ≈ -1/2
    esys.sys_origin.dipoles[2][3] ≈ -1/2

    set_field!(esys, [0, 0, -10])
    randomize_spins!(esys)
    minimize_energy!(esys)
    esys.sys_origin.dipoles[1][3] ≈ 1/2
    esys.sys_origin.dipoles[2][3] ≈ 1/2

    set_field!(esys, [0, 0, 0])
    randomize_spins!(esys)
    minimize_energy!(esys; g_tol=1e-14)
    norm(esys.sys_origin.dipoles[1]) < 1e-14
    norm(esys.sys_origin.dipoles[2]) < 1e-14

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

    set_field!(esys, [0, 0, 0])
    for unit in Sunny.eachunit(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], unit)
    end
    swt = SpinWaveTheory(esys; measure=Sunny.empty_measurespec(sys), regularization=0.0)
    disp = dispersion(swt, qpts)
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

    # Julia 1.11 slightly changes the SVD, which leads to a different
    # decomposition of "general" pair interactions. Small numerical differences
    # can be amplified over a long dynamical trajectory.
    @static if VERSION == v"1.10"
        intensities_ref = [-0.0009151914880177183 -0.002340042772876626 -0.0037316702420023395 -0.002340042772876612 -0.0009151914880177183; 0.02116154679963919 0.05967703767242055 0.009161242329689078 0.03233210977077469 0.02116154679963919; -0.04206086791537969 -0.2227189243447692 -0.12040711429107293 -0.025389877189384635 -0.04206086791537969; 0.3962131659090294 2.7150236439903197 4.3968225602687 0.5625477581811109 0.3962131659090294; 8.48497330044607 17.689468342820923 14.169270995083266 3.1979660472637534 8.48497330044607; 17.1145045637742 21.877166922282075 10.751399167760733 3.8125499846147113 17.1145045637742; 8.490961712258231 5.100854031990811 0.6862255516496436 0.8503635473494539 8.490961712258231; -0.3334143274663385 -0.966097718147396 -0.1252420669485946 -0.2868227954981837 -0.3334143274663385; 0.07518063142323421 0.2337350060517672 0.13638239262286792 0.012717547170902647 0.07518063142323421; -0.09925761646769077 -0.29077390054022495 -0.03237049466881591 -0.290773900540224 -0.09925761646769077] 
    elseif VERSION >= v"1.11"
        intensities_ref = [-0.005898711265804917 0.00037194381298070534 -0.007139302988678504 0.0003719438129806763 -0.005898711265804917; -0.021174997374043665 0.04366825611810105 0.007169734714023062 0.05743679447787905 -0.021174997374043665; -0.0898921215977289 -0.026851190523929185 -0.21509883592744028 -0.13293335904345077 -0.0898921215977289; 0.07554637005717268 1.0310198342691788 7.6305129282212 2.0138084437390957 0.07554637005717268; 9.029181600659825 4.6645584970233935 26.007406182009305 12.508876537699672 9.029181600659825; 19.33194294971418 4.561683298052969 21.585685225141788 15.176512232105381 19.33194294971418; 9.190329239958585 0.5342077291704622 2.3267523495368336 3.3977085221103662 9.190329239958585; -0.5730586490317443 -0.36243855204206615 -0.2377418272654857 -0.7846957436907301 -0.5730586490317443; 0.4227065393257319 -0.0978593259945275 0.5472087416056827 0.03159040448673515 0.4227065393257319; -0.08949094320156001 -0.3280980521751531 0.015933918006198686 -0.32809805217515364 -0.08949094320156001]
    end
    @test isapprox(res.data, intensities_ref)
end

# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end