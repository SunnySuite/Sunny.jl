@testitem "Canted AFM" begin
    function test_canted_afm(S)
        J, D, h = 1.0, 0.54, 0.76
        a = 1
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)
        q = [0.12, 0.23, 0.34]

        dims = (1, 1, 1)
        sys = System(cryst, dims, [SpinInfo(1; S, g=1)], :dipole; units=Units.theory)
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
        set_external_field!(sys, [0, 0, h])

        n = [0.,0.,1]
        Sunny.check_rotational_symmetry(sys; n, θ=0.01)

        c₂ = 1 - 1/(2S)
        θ = acos(h / (2S*(4J+D*c₂)))
        sys.dipoles[1]= [sin(θ), 0, cos(θ)]
        k = [0.5,0.5,0]
        swt = SpinWaveTheory(sys)
        formula = Sunny.intensity_formula_SingleQ(swt,k,n, :perp; kernel=delta_function_kernel)
        disp, _ = Sunny.intensities_bands_SingleQ(swt, [q], formula)
        ϵq_num = disp[1,1]

        #  Analytical
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
    crystal = Crystal(latvecs, [[0.24964,0,0.5]], 150)
    latsize = (1,1,1)
    sys = System(crystal, latsize, [SpinInfo(1; S=5/2, g=2)], :dipole; seed=5)
    set_exchange!(sys, 0.85,  Bond(3, 2, [1,1,0]))   # J1
    set_exchange!(sys, 0.24,  Bond(1, 3, [0,0,0]))   # J2
    set_exchange!(sys, 0.053, Bond(2, 3, [-1,-1,1])) # J3
    set_exchange!(sys, 0.017, Bond(1, 1, [0,0,1]))   # J4
    set_exchange!(sys, 0.24,  Bond(3, 2, [1,1,1]))   # J5
    
    n = [0.,0.,1]
    randomize_spins!(sys)
    xmin = [-1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6] # Minimum value of x
    xmax = [2π 2π 2π 1e-6 1e-6 1]  # Maximum value of x
    x0 = [3.2 0.4 3.2 0. 0. 0.3]  # Initial value of x
    k = Sunny.optimagstr(x -> Sunny.gm_planar!(sys, n, x), xmin, xmax, x0)

    swt = SpinWaveTheory(sys)
    formula = Sunny.intensity_formula_SingleQ(swt,k,n, :perp; kernel=delta_function_kernel)
    q = [0.41568,0.56382,0.76414]
    disp, intensity = Sunny.intensities_bands_SingleQ(swt, [q], formula);
    SpinW_energies = [2.6267,3.8202,3.9422,2.8767,3.9095,4.4605,3.31724,4.0113,4.7564]
    SpinW_intensities = [0.484724856017038, 0.962074686407910, 0.0932786148844105, 0.137966379056292, 0.0196590925454593, 2.37155695274281, 2.21507666401705, 0.118744173882554, 0.0547564956435423]
    @test isapprox(disp[:], reverse(SpinW_energies); atol=1e-3)
    @test isapprox(SpinW_intensities/Sunny.natoms(crystal), intensity[:]; atol=1e-3)    
end
