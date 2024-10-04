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
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
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
    field_offset = -2*(Sl[3] + Su[3]) # 2 for g-factor
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

    # Test dispersion against analytical formula
    ω_ref(q, J, J′) = J*sqrt(1 + 2(J′/J) * cos(2π*q))

    qs = [[0, 1.0, 0], []]

end

# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end