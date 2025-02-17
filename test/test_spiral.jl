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
    sys = System(cryst, [1 => Moment(s=1, g=1)], :dipole, seed=0)
    randomize_spins!(sys)
    enable_dipole_dipole!(sys, 1.0)
    E1 = spiral_energy_per_site(sys; k, axis)

    # compute ewald energy using supercell
    sys_large = repeat_periodically_as_spiral(sys, (La*Na, Lb*Nb, Lc*Nc); k, axis)
    enable_dipole_dipole!(sys_large, 1.0)
    E2 = energy_per_site(sys_large)

    @test E1 ≈ E2

    ################ test random case ###################
    La, Lb, Lc = (2, 3, 4)
    k = [1/La, 1/Lb, 1/Lc]

    Na, Nb, Nc = (1, 2, 3)
    axis = normalize(randn(3))

    # compute ewald energy using J(Q)
    sys = System(cryst, [1 => Moment(s=1, g=1)], :dipole, seed=0)
    randomize_spins!(sys)
    enable_dipole_dipole!(sys, 1.0)
    E1 = spiral_energy_per_site(sys; k, axis)

    sys_large = repeat_periodically_as_spiral(sys, (La*Na, Lb*Nb, Lc*Nc); k, axis)
    enable_dipole_dipole!(sys_large, 1.0)
    E2 = energy_per_site(sys_large)

    @test E1 ≈ E2
end


@testitem "Canted AFM" begin
    s = 3/2
    J, D, h = 1.0, 0.54, 0.76
    rcs = Sunny.rcs_factors(s)[2]
    a = 1
    latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    
    sys = System(cryst, [1 => Moment(; s, g=-1)], :dipole)
    set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
    set_onsite_coupling!(sys, S -> (D/rcs)*S[3]^2, 1)
    set_field!(sys, [0, 0, h])
    
    k = Sunny.minimize_luttinger_tisza_exchange(sys; k_guess=randn(3))
    @test k[1:2] ≈ [0.5, 0.5]
    
    axis = [0, 0, 1]
    randomize_spins!(sys)
    k = minimize_spiral_energy!(sys, axis; k_guess=randn(3))
    @test k[1:2] ≈ [0.5, 0.5]
    @test isapprox(only(sys.dipoles)[3], h / (8J + 2D); atol=1e-6)
    
    q = [0.12, 0.23, 0.34]
    swt = SpinWaveTheorySpiral(sys; measure=ssf_trace(sys), k, axis)
    res = intensities_bands(swt, [q])
    
    # Analytical check on dispersion
    
    θ = acos(h / (2s*(4J+D)))
    Jq = 2J*(cos(2π*q[1])+cos(2π*q[2]))
    disp_ref = real(√Complex(4J*s*(4J*s+2D*s*sin(θ)^2) + cos(2θ)*(Jq*s)^2 + 2s*Jq*(4J*s*cos(θ)^2 + D*s*sin(θ)^2)))
    @test res.disp[1] ≈ res.disp[2] ≈ disp_ref
    
    # Same calculation with supercell
    
    sys2 = repeat_periodically_as_spiral(sys, (2, 2, 1); k, axis)
    swt2 = SpinWaveTheory(sys2; measure=ssf_trace(sys2))
    res2 = intensities_bands(swt2, [q])
    @test res.disp[1] ≈ res.disp[2] ≈ res2.disp[2]
    @test res.data[1] + res.data[2] ≈ res2.data[2]
    @test res.disp[end] ≈ res2.disp[end]
    @test res.data[end] ≈ res2.data[end]    
end


@testitem "Langasite" begin
    a = b = 8.539
    c = 5.2414
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    crystal = Crystal(latvecs, [[0.24964, 0, 0.5]], 150)
    sys = System(crystal, [1 => Moment(s=5/2, g=2)], :dipole; seed=1)
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
    k = minimize_spiral_energy!(sys, axis; k_guess=randn(3))
    @test spiral_energy(sys; k, axis) ≈ -16.356697120589477

    # There are two possible chiralities. Select just one.
    @test isapprox(k, k_ref; atol=1e-6) || isapprox(k, k_ref_alt; atol=1e-6)
    if isapprox(k, k_ref_alt; atol=1e-6)
        sys.dipoles[[1,2]] = sys.dipoles[[2,1]]
        k = k_ref
    end
    @test spiral_energy(sys; k, axis) ≈ -16.356697120589477

    measure = ssf_perp(sys; apply_g=false)
    swt = SpinWaveTheorySpiral(sys; measure, k, axis)
    q = [0.41568,0.56382,0.76414]
    res = intensities_bands(swt, [q])

    disp_spinw = [4.7564, 4.0113, 3.31724, 4.4605, 3.9095, 2.8767, 3.9422, 3.8202, 2.6267]
    data_spinw = [0.484724856017038, 0.962074686407910, 0.0932786148844105, 0.137966379056292, 0.0196590925454593, 2.37155695274281, 2.21507666401705, 0.118744173882554, 0.0547564956435423]
    # Sunny dispersion is sorted in descending order
    P = sortperm(disp_spinw; rev=true)
    @test isapprox(res.disp[:, 1], disp_spinw[P]; atol=1e-3)
    @test isapprox(res.data[:, 1], data_spinw[P]; atol=2e-3)

    disp = dispersion(swt, [q])
    @test res.disp ≈ disp
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
    
    sys = System(cryst, [1 => Moment(s=5/2, g=2)], :dipole, seed=0)
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
    k = minimize_spiral_energy!(sys, axis; k_guess=randn(3))
    E = Sunny.luttinger_tisza_exchange(sys; k)
    @test isapprox(E, E_ref; atol=1e-12)
    @test isapprox(k, k_ref; atol=1e-6) || isapprox(k, k_ref_alt; atol=1e-6)    
end
