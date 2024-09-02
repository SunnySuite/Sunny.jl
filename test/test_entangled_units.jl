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
        contracted_crystal, contraction_info = contract_crystal(crystal, units)
        expanded_crystal = Sunny.expand_crystal(contracted_crystal, contraction_info)
        @test expanded_crystal.positions ≈ crystal.positions
    end
end

# TODO: Add reshapings tests.

# TODO: Add more complicated interaction tests (e.g., for Ba3Mn2O8).

# TODO: Add test with magnetic unit cell larger than a single unit (i.e. not q=0
# ordering).

@testitem "Test Dimer Interactions" begin
    J = 1.0
    J′ = 0.1
    dims = (1, 1, 1)
    latvecs = [
        1  0  0
        0  1  0
        0  0  2
    ]
    positions = [[0, 0, 0], [0.0, 0.5, 0.0]] 

    crystal = Crystal(latvecs, positions, 1; types = ["A", "B"])
    sys = System(crystal, dims, [SpinInfo(1; S=1/2, g=2), SpinInfo(2; S=1/2, g=2)], :SUN)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, J′, Bond(2, 2, [1, 0, 0]))  # Needed because we broke the symmetry equivalence of the two sites

    sys_entangled = EntangledSystem(sys, [(1, 2)])
    interactions = sys_entangled.sys.interactions_union[1]

    # Test on-bond exchange
    onsite_operator = interactions.onsite
    S = spin_matrices(1/2)
    Sl, Su = to_product_space(S, S)
    onsite_ref = J * (Sl' * Su)
    @test onsite_operator ≈ onsite_ref

    # Test inter-bond exchange
    pc = Sunny.as_general_pair_coupling(interactions.pair[1], sys_entangled.sys)
    Sl1, Sl2 = to_product_space(Sl, Sl)
    Su1, Su2 = to_product_space(Su, Su)
    bond_operator = zeros(ComplexF64, 16, 16)
    for (A, B) in pc.general.data
        bond_operator .+= kron(A, B)
    end
    bond_ref = J′*((Sl2' * Sl1) .+ (Su2' * Su1))
    @test bond_operator ≈ bond_ref
end

@testitem "Ba3Mn2O8 Dispersion and Intensities" begin
    function J_of_Q(_, J1, J2, J3, J4, q)
        h, k, l = q
        ω1(h, k, l) = cos((2π/3)*(-h + k + l)) + cos((2π/3)*(-h - 2k + l)) + cos((2π/3)*(2h + k + l))
        ω2(h, k, l) = cos(2π*k) + cos(2π*(h + k)) + cos(2π*h)
        ω4(h, k, l) = cos((2π/3)*(2h - 2k + l)) + cos((2π/3)*(2h + 4k + l)) + cos((2π/3)*(-4h - 2k + l))
        -J1*ω1(h, k, l) + 2(J2 - J3)*ω2(h, k, l) - J4*ω4(h, k, l)
    end

    function disp0(q)
        J0, J1, J2, J3, J4, D = 1.642, 0.118, 0.256, 0.142, 0.037, -0.032
        Δ0 = J0 + 2D/3
        sqrt(Δ0^2 + (8/3)*Δ0*J_of_Q(J0, J1, J2, J3, J4, q) )
    end

    function dispm(q)
        J0, J1, J2, J3, J4, D = 1.642, 0.118, 0.256, 0.142, 0.037, -0.032
        Δm = J0 - D/3
        sqrt(Δm^2 + (8/3)*Δm*J_of_Q(J0, J1, J2, J3, J4, q) )
    end

    function Mn_crystal()
        latvecs = [5.710728 -2.855364 0.0; 0.0 4.945635522103099 0.0; 0.0 0.0 21.44383] 
        positions = [[1/3, 2/3, 0.07374666666666664], [1/3, 2/3, 0.25958666666666663], [0.0, 0.0, 0.40708], [0.0, 0.0, 0.59292], [2/3, 1/3, 0.7404133333333334], [2/3, 1/3, 0.9262533333333333]]
        symprec = 0.01
        return Crystal(latvecs, positions; symprec)
    end

    function Ba3Mn2O8(; dims=(1, 1, 1), g=1.0, D = -0.032, J0 = 1.642, J1 = 0.118, J2 = 0.256, J3 = 0.142, J4 = 0.037)
        sys = System(Mn_crystal(), dims, [SpinInfo(1; S=1, g)], :SUN)
        S = spin_matrices(1)
        S1, S2 = to_product_space(S, S)
        set_pair_coupling!(sys, J0 * (S1' * S2), Bond(1, 2, [0, 0, 0]))
        set_pair_coupling!(sys, J1 * (S1' * S2), Bond(3, 2, [0, 0, 0]))
        set_pair_coupling!(sys, J2 * (S1' * S2), Bond(1, 1, [1, 0, 0]))
        set_pair_coupling!(sys, J3 * (S1' * S2), Bond(1, 2, [1, 0, 0]))
        set_pair_coupling!(sys, J4 * (S1' * S2), Bond(3, 2, [1, 0, 0]))
        set_onsite_coupling!(sys, D*S[3]^2, 1)
        return sys
    end

    # Set up system and reshape into primitive cell
    crystal = Mn_crystal()
    sys = Ba3Mn2O8()
    shape = [-1 -1 2/3; 0 -1 1/3; 0 0 1/3]
    sys_reshaped = reshape_supercell(sys, shape)

    # Make entangled system and spin wave theory
    sys_entangled = EntangledSystem(sys_reshaped, [(1, 2)])
    set_coherent!(sys_entangled.sys, Sunny.CVec{9}(0.0, 0.0, √3/3, 0.0, -√3/3, 0.0, √3/3, 0, 0), (1, 1, 1, 1))
    swt = EntangledSpinWaveTheory(sys_entangled)

    # Calculation dispersion and intensities
    FWHM = 0.295

    points = [
        [0.175, 0.175, 1.5],
        [0.85, 0.85, 1.5],
        [0.85, 0.85, 3],
        [0.0, 0.0, 3],
        [0.0, 0.0, 8],
    ]
    qpts = q_space_path(crystal, points, 20)
    energies = range(0.0, 4.0, 10)

    is = intensities(swt, qpts; energies; kernel=gaussian(; fwhm=FWHM))
    is_bands = intensities_bands(swt, qpts)

    d0_ref = disp0.(path)
    d0 = disp[:,8]
    band_err1 = √sum((d0_ref .- d0) .^ 2)
    @test band_err1 < 1e-7

    dm_ref = dispm.(path)
    dm = disp[:,6]
    band_err2 = √sum((dm_ref .- dm) .^ 2)
    @test band_err2 < 1e-7

    is_ref = [1.5318559039289553e-54 5.4757218095775535e-33 7.357403275541178e-17 3.9795517317276186e-6 0.9020338263369125 0.8314157146337926 2.911058313640675e-6 3.665251511398035e-17 1.6108014399323818e-33 1.2770658651600437e-42; 6.344775971921163e-16 1.8297903880834307e-5 2.186860140613253 1.125414819762803 2.372670614798904e-6 1.8884749639385228e-17 5.364594261005178e-34 5.293626361366888e-56 1.5507874038885882e-56 1.5731739245245433e-42; 1.916445283080627e-22 1.5112603523543528e-9 0.046321847672565045 5.889951708200119 0.0031331962838944093 6.5822336276008166e-12 5.1035843246680587e-26 1.398942793394379e-45 1.4068300720405523e-55 1.0612068475874395e-41; 1.0342024392072121e-16 6.295383373660083e-6 1.5314917741949632 1.5884648898847082 6.9292015386003456e-6 1.1783549059185498e-16 7.29390986501124e-33 1.5826736515812515e-54 7.205550912747091e-56 5.5831854222249985e-42; 7.816369005435853e-65 4.958135881093393e-41 1.1004603393444979e-22 8.837745291621732e-10 0.02729382866823639 3.4792844815914976 0.001870907449751048 4.029919472431617e-12 3.236721779435184e-26 1.258781799195952e-42; 4.1286086857912686e-32 4.109912848294968e-16 1.53865152859326e-5 2.312737091532787 1.4474386749428667 3.6605144080508623e-6 3.5053122616269065e-17 1.205982396647456e-33 1.888239550856115e-55 3.384617657031145e-42; 1.0956905921250642e-15 3.09807215254437e-5 3.6242252907637944 1.8260559377579642 3.77624592541727e-6 2.953182889984101e-17 8.250224084439066e-34 8.009410306835117e-56 4.543889223635611e-56 5.863106256174803e-42; 2.4641979663556544e-21 9.259498870544819e-9 0.1405047201368198 8.991804007003132 0.0023532144516871254 2.345225668006157e-12 8.408440394132043e-27 1.0520889028677954e-46 1.508356309761414e-55 1.1755786936605455e-41; 5.9006127773674436e-21 1.6083917235355553e-8 0.18138815381567364 8.612130243663842 0.0016307094799296721 1.1473990136958224e-12 2.86349056374958e-27 2.4776173823909797e-47 1.4600597514018359e-56 1.2635254358709617e-42; 3.534719908717514e-16 1.5534826875777495e-5 2.933172504328625 2.2501936681347714 6.4727258857098395e-6 6.621955526402235e-17 2.3498631523198714e-33 2.862774376331706e-55 9.224249931423226e-57 1.3867670560912468e-42; 2.9056707554206643e-34 1.1389749648037732e-17 1.8359864206642449e-6 1.1476609410649823 2.6195576292904894 2.104174627891465e-5 5.842573811256141e-16 5.564262837339542e-32 1.8126351339556568e-53 9.514020434069003e-43; 4.95669586498486e-51 4.090111065132001e-30 1.1512703563038524e-14 0.00011053986637697482 3.620428209362293 0.40448305917267885 1.541487919718299e-7 2.0039150465080363e-19 8.886236730407736e-37 8.03538987497683e-43; 1.6256804884099134e-104 1.626131554260588e-73 5.548499522448917e-48 6.457952119615799e-28 2.563973194366388e-13 0.0003472414267885618 1.6041647745239174 0.025279344649153018 1.358882914693268e-9 2.491710931378807e-22; 1.4445317508683965e-67 1.425287573157436e-43 4.79708370999404e-25 5.507461722013346e-12 0.00021568754619639152 0.02881363456352638 1.3130171839110966e-5 2.0409960671735094e-14 1.0822144233620588e-28 1.799465203023977e-45; 9.225789423317937e-80 5.21566646643929e-53 1.0058089054657094e-31 6.616381427047744e-16 1.4846530782284771e-5 1.1363937177899095 0.29670991367977784 2.64261978483415e-7 8.028542283600199e-19 8.320283660276395e-36; 1.0913958451883772e-93 2.04937342166082e-64 1.312681289078253e-40 2.868117121788508e-22 2.1376365239065535e-9 0.05434633134660532 4.713092846351799 0.0013942515956084163 1.4069403885772203e-12 4.8429461011708685e-27]
    @test isapprox(is, is_ref; atol=1e-12)
end