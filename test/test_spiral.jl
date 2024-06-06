@testitem "Ewald energy" begin
    using LinearAlgebra

    latvecs = lattice_vectors(5.254,11.65,16.15, 90, 90.8, 90)
    positions = [[1/4, 0.1178, 0], [3/4, 0.8822, 0],[1/4, 0.3822, 1/2],[3/4, 0.6178, 1/2]]
    cryst = Crystal(latvecs, positions)

    ################ test special case ###################
    La = 2
    Lb = 1
    Lc = 2
    k = [1/La, 1/Lb, 1/Lc]

    Na, Nb, Nc = (1, 2, 3)
    axis = normalize(randn(3))

    # compute ewald energy using J(k)
    sys = System(cryst, (1, 1, 1), [SpinInfo(1, S=1, g=1)], :dipole, seed=0)
    randomize_spins!(sys)
    enable_dipole_dipole!(sys)
    E1 = Sunny.spiral_energy_per_site(sys; k, axis)

    # compute ewald energy using supercell
    sys_large = System(cryst, (La*Na, Lb*Nb, Lc*Nc), [SpinInfo(1, S=1, g=1)], :dipole, seed=0)
    for i in 1:Sunny.natoms(sys.crystal)
        set_spiral_order_on_sublattice!(sys_large, i; k, axis, S0=sys.dipoles[i])
    end
    enable_dipole_dipole!(sys_large)
    E2 = energy_per_site(sys_large)

    @test E1 ≈ E2

    ################ test random case ###################
    La, Lb, Lc = (2, 3, 4)
    k = [1/La, 1/Lb, 1/Lc]

    Na, Nb, Nc = (1, 2, 3)
    axis = normalize(randn(3))

    # compute ewald energy using J(Q)
    sys = System(cryst, (1, 1, 1), [SpinInfo(1, S=1, g=1)], :dipole, seed=0)
    randomize_spins!(sys)
    enable_dipole_dipole!(sys)
    E1 = Sunny.spiral_energy_per_site(sys; k, axis)

    # compute ewald energy using supercell
    sys_large = System(cryst, (La*Na, Lb*Nb, Lc*Nc), [SpinInfo(1, S=1, g=1)], :dipole, seed=0)
    for i in 1:Sunny.natoms(sys.crystal)
        set_spiral_order_on_sublattice!(sys_large, i; k, axis, S0=sys.dipoles[i])
    end
    enable_dipole_dipole!(sys_large)
    E2 = energy_per_site(sys_large)

    @test E1 ≈ E2
end


@testitem "Canted AFM" begin
    function test_canted_afm(S)
        J, D, h = 1.0, 0.54, 0.76
        a = 1
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)

        dims = (1, 1, 1)
        sys = System(cryst, dims, [SpinInfo(1; S, g=1)], :dipole; units=Units.theory)
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
        set_external_field!(sys, [0, 0, h])

        k = Sunny.minimize_luttinger_tisza_exchange(sys; k_guess=randn(3))
        @test k[1:2] ≈ [0.5, 0.5]

        axis = [0, 0, 1]
        randomize_spins!(sys)
        k = Sunny.minimize_energy_spiral!(sys, axis; k_guess=randn(3))
        @test k[1:2] ≈ [0.5, 0.5]
        c₂ = 1 - 1/2S
        @test only(sys.dipoles)[3] ≈ h / (8J + 2D*c₂)

        q = [0.12, 0.23, 0.34]
        swt = SpinWaveTheory(sys)
        formula = Sunny.intensity_formula_spiral(swt, :perp; k, axis, kernel=delta_function_kernel)
        disp, _ = intensities_bands(swt, [q], formula)
        ϵq_num = disp[1,1]

        # Analytical
        c₂ = 1 - 1/(2S)
        θ = acos(h / (2S*(4J+c₂*D)))
        Jq = 2J*(cos(2π*q[1])+cos(2π*q[2]))
        ϵq_ana = real(√Complex(4J*S*(4J*S+2D*S*c₂*sin(θ)^2) + cos(2θ)*(Jq*S)^2 + 2S*Jq*(4J*S*cos(θ)^2 + c₂*D*S*sin(θ)^2)))

        @test ϵq_num ≈ ϵq_ana
    end

    test_canted_afm(1)
    test_canted_afm(2)
end


@testitem "Langasite" begin
    a = b = 8.539
    c = 5.2414
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    crystal = Crystal(latvecs, [[0.24964, 0, 0.5]], 150)
    sys = System(crystal, (1,1,1), [SpinInfo(1; S=5/2, g=2)], :dipole; seed=1)
    set_exchange!(sys, 0.85,  Bond(3, 2, [1,1,0]))   # J1
    set_exchange!(sys, 0.24,  Bond(1, 3, [0,0,0]))   # J2
    set_exchange!(sys, 0.053, Bond(2, 3, [-1,-1,1])) # J3
    set_exchange!(sys, 0.017, Bond(1, 1, [0,0,1]))   # J4
    set_exchange!(sys, 0.24,  Bond(3, 2, [1,1,1]))   # J5

    k_ref = [0, 0, 0.14264604656200577]
    k_ref_alt = [0, 0, 1] - k_ref

    k = Sunny.minimize_luttinger_tisza_exchange(sys; k_guess=randn(3))
    @test isapprox(k, k_ref; atol=1e-6) || isapprox(k, k_ref_alt; atol=1e-6)

    axis = [0, 0, 1]
    randomize_spins!(sys)
    k = Sunny.minimize_energy_spiral!(sys, axis; k_guess=randn(3))
    @test Sunny.spiral_energy(sys; k, axis) ≈ -16.356697120589477

    # There are two possible chiralities. Select just one.
    @test isapprox(k, k_ref; atol=1e-6) || isapprox(k, k_ref_alt; atol=1e-6)
    if isapprox(k, k_ref_alt; atol=1e-6)
        sys.dipoles[[1,2]] = sys.dipoles[[2,1]]
        k = k_ref
    end
    @test Sunny.spiral_energy(sys; k, axis) ≈ -16.356697120589477

    swt = SpinWaveTheory(sys)
    formula = Sunny.intensity_formula_spiral(swt, :perp; k, axis, kernel=delta_function_kernel)
    q = [0.41568,0.56382,0.76414]
    disp, intens = intensities_bands(swt, [q], formula)
    # TODO: why reverse?
    disp_spinw = reverse([2.6267,3.8202,3.9422,2.8767,3.9095,4.4605,3.31724,4.0113,4.7564])
    intens_spinw = [0.484724856017038, 0.962074686407910, 0.0932786148844105, 0.137966379056292, 0.0196590925454593, 2.37155695274281, 2.21507666401705, 0.118744173882554, 0.0547564956435423]
    # Sunny dispersion bands appear in descending order. Sort SpinW calculation
    # in the same way.
    P = sortperm(disp_spinw; rev=true)
    @test isapprox(disp[1,:], disp_spinw[P]; atol=1e-3)
    @test isapprox(intens[1,:], intens_spinw[P]/Sunny.natoms(crystal); atol=1e-3)    
end


@testitem "BFMO" begin
    latvecs = lattice_vectors(5.254, 11.65, 16.15, 90, 90.8, 90)
    positions = [[1/4, 0.1178, 0], [3/4, 0.8822, 0], [1/4, 0.3822, 1/2], [3/4, 0.6178, 1/2]]
    cryst = Crystal(latvecs, positions)
    bond1 = Bond(3, 4, [0, 0, 0])
    bond2 = Bond(1, 1, [1, 0, 0])
    bond3 = Bond(3, 4, [1, 0, 0])
    bond4 = Bond(1, 3, [0, 0, 0])
    bond5 = Bond(1, 2, [0, 0, 0])
    bond6 = Bond(1, 3, [1, 0, 0])
    bond7 = Bond(1, 3, [-1, 0, 0])
    
    ka = 0.4413
    kc = 0.1851
    k_ref = [ka, 0, kc]
    k_ref_alt = [1-ka, 0, 1-kc]
    
    J1 = 0.4442
    J3 = 0.07762
    J4 = -0.00775
    J5 = 0.00100
    J7 = -0.00461
    
    # Constrain J6 and J2 to make k_ref exact
    J6 = (sin((-2*ka+kc)*pi)*J7+sin(kc*pi)*J4)/(-sin((2*ka+kc)*pi))
    J2 = (2*sin((2*ka+kc)*pi)*J6-sin(ka*pi)*(J1+J5)-3*sin(3*ka*pi)*J3-(sin(2*pi*(ka-kc/2+1/2))-sin(2*pi*(-ka+kc/2+1/2)))*J7)/(-2*sin(2*ka*pi))
    
    # Exact energy reference
    E_ref = (-((J1+J5)*cos(pi*ka))+J2*cos(2*pi*ka)-J3*cos(3*pi*ka)+J4*cos(pi*kc)+J6*cos(pi*(2*ka+kc))+J7*cos(pi*(2*ka-kc)))*5/2*(5/2)*4
    
    sys = System(cryst, (1,1,1), [SpinInfo(1,S=5/2,g=2)], :dipole, seed=0)
    set_exchange!(sys, J1, bond1)
    set_exchange!(sys, J2, bond2)
    set_exchange!(sys, J3, bond3)
    set_exchange!(sys, J4, bond4)
    set_exchange!(sys, J5, bond5)
    set_exchange!(sys, J6, bond6)
    set_exchange!(sys, J7, bond7)
    
    # Unfortunately, randomizing the initial guess will lead to occasional
    # optimization failures. TODO: Use ForwardDiff to improve accuracy.
    k_guess = [0.2, 0.4, 0.8]
    k = Sunny.minimize_luttinger_tisza_exchange(sys; k_guess)
    E = Sunny.luttinger_tisza_exchange(sys; k)
    @test isapprox(E, E_ref; atol=1e-12)
    @test isapprox(k, k_ref; atol=1e-5) || isapprox(k, k_ref_alt; atol=1e-5)
    
    axis = [0, 0, 1]
    randomize_spins!(sys)
    k = Sunny.minimize_energy_spiral!(sys, axis; k_guess=randn(3))
    E = Sunny.luttinger_tisza_exchange(sys; k)
    @test isapprox(E, E_ref; atol=1e-12)
    @test isapprox(k, k_ref; atol=1e-6) || isapprox(k, k_ref_alt; atol=1e-6)    
end
