@testitem "Kitchen Sink" begin
    using LinearAlgebra

    # Pyrochlore with nonstandard, primitive lattice vectors
    latvecs = [[1, 1, 0] [1, 0, 1] [0, 1, 1]] / 2
    positions = [
        [5, 5, 1]/8,
        [5, 1, 5]/8,
        [1, 5, 5]/8,
        [5, 5, 5]/8,
    ]
    cryst = @test_warn "Lattice vectors are not right-handed." Crystal(latvecs, positions)
    natoms = Sunny.natoms(cryst)

    g = 7.2 # For comparison with previous test data
    infos = [SpinInfo(1; S=5/2, g)]
    sys = System(cryst, (1, 1, 1), infos, :SUN; seed=0)

    A,B,C,D = 2.6, -1.3, 0.2, -5.7
    set_exchange!(sys, [A C -D; C A D; D -D B], Bond(1, 2, [0, 0, 0]))

    A,B,C,D,E,F,G,H,K = 2.6, -1.3, 0.2, -5.7, 8.2, 0.3, 2.5, -0.6, 1.3
    set_exchange!(sys, [A F+K E-H; F-K B D+G; E+H D-G C], Bond(1, 4, [1, 0, 0]))

    A,B,C,D = 2.6, -1.3, 0.2, -5.7
    set_exchange!(sys, [A D D; D B C; D C B], Bond(4, 4, [1, 1, 0]))

    O = stevens_matrices(spin_label(sys, 3))
    c1, c2, c3 = 2.6, -1.3, 0.2, -5.7
    Λ = c1 * (O[2,-2] - 2O[2,-1] - 2O[2,1]) +
        c2 * (-7O[4,-3] + 2O[4,-2] + O[4,-1] + O[4,1] + 7O[4,3]) +
        c3 * (O[4,0] + 5O[4,4])
    set_onsite_coupling!(sys, Λ, 3)

    A = [1 3 1; -1 1 0; 0 0 1]
    sys = reshape_supercell(sys, A)

    # minimize_energy!(sys; maxiters=1000)
    # println(round(reinterpret(reshape, ComplexF64, sys.coherents); digits=3))
    ground_state = ComplexF64[0.452 + 0.069im; 0.586 + 0.166im; 0.31 + 0.351im; 0.035 + 0.39im; -0.136 + 0.145im; -0.03 - 0.076im;;;;; -0.026 - 0.064im; 0.247 - 0.209im; -0.444 + 0.307im; 0.501 - 0.122im; -0.478 - 0.012im; 0.32 - 0.046im;;;;; -0.078 - 0.451im; -0.029 - 0.609im; 0.234 - 0.406im; 0.358 - 0.158im; 0.181 + 0.083im; -0.063 + 0.053im;;;;; -0.061 + 0.033im; -0.235 - 0.222im; 0.354 + 0.407im; -0.177 - 0.485im; 0.042 + 0.476im; -0.082 - 0.312im;;;;; -0.782 + 0.332im; 0.184 - 0.453im; 0.058 + 0.138im; 0.121 + 0.049im; -0.007 + 0.003im; -0.006 + 0.015im;;;;; 0.057 - 0.005im; -0.023 - 0.013im; -0.101 + 0.009im; -0.247 - 0.257im; -0.631 - 0.193im; -0.617 + 0.209im;;;;; -0.695 - 0.489im; 0.482 - 0.083im; -0.087 + 0.122im; 0.022 + 0.129im; -0.006 - 0.004im; -0.016 + 0.003im;;;;; 0.04 + 0.041im; -0.004 - 0.026im; -0.07 - 0.072im; 0.044 - 0.354im; -0.248 - 0.611im; -0.551 - 0.347im;;;;; 0.347 + 0.149im; -0.216 + 0.511im; -0.386 - 0.36im; 0.434 - 0.153im; -0.024 + 0.24im; 0.026 - 0.021im;;;;; 0.061 + 0.082im; 0.096 + 0.219im; 0.361 + 0.203im; 0.487 - 0.18im; 0.228 - 0.554im; -0.282 - 0.229im;;;;; -0.205 + 0.317im; -0.468 - 0.298im; 0.419 - 0.32im; 0.078 + 0.454im; -0.233 - 0.064im; 0.016 + 0.029im;;;;; 0.101 + 0.018im; 0.22 + 0.095im; 0.402 - 0.097im; 0.233 - 0.464im; -0.211 - 0.561im; -0.363 + 0.025im;;;;; 0.059 - 0.582im; 0.442 - 0.355im; 0.51 + 0.063im; 0.158 + 0.213im; -0.009 + 0.024im; 0.03 - 0.026im;;;;; 0.582 - 0.058im; 0.436 + 0.363im; 0.04 + 0.512im; -0.177 + 0.197im; -0.025 - 0.004im; 0.031 + 0.024im;;;;; 0.561 + 0.165im; 0.268 + 0.5im; -0.156 + 0.49im; -0.238 + 0.116im; -0.022 - 0.013im; 0.02 + 0.034im;;;;; -0.337 + 0.478im; -0.56 + 0.092im; -0.413 - 0.306im; -0.033 - 0.263im; 0.019 - 0.017im; -0.039 + 0.008im]
    
    sys.coherents .= reinterpret(reshape, Sunny.CVec{6}, ground_state)
    minimize_energy!(sys)
    @assert energy_per_site(sys) ≈ -328.38255

    # Verify that this is a local minimum of energy
    @test norm(Sunny.proj.(Sunny.energy_grad_coherents(sys), sys.coherents)) < 1e-8

    # Test energies at an arbitrary wave vector
    ks = [[0.24331089495721447, 0.2818361515716459, 0.21954858411037714]]
    swt = SpinWaveTheory(sys)

    formula = intensity_formula(swt, :perp, kernel=delta_function_kernel)
    disps, is = intensities_bands(swt, ks, formula)

    # println(round.(disps; digits=12))
    # println(round.(is / (natoms * g^2); digits=12))
    disps_golden = [1394.440092589898 1393.728009961772 1393.008551261241 1392.919524984054 1279.23991907861 1279.094568482213 1278.22451852538 1277.691761492894 1194.366336265056 1193.750083635324 1191.583519669771 1189.794451350786 1131.422439597908 1131.202770084574 1065.242927860663 1065.095892456642 1026.649340932941 1024.028348568104 1022.830406309672 1020.767349649408 945.202397540707 944.795817861582 835.545028404179 832.001588705254 827.939501419137 827.307586957877 821.216582176438 820.430993577092 820.294548786087 818.594571008014 810.207001100209 808.553158283705 766.524411081004 766.516102760229 766.513825862473 766.508655568197 758.579854177684 754.683765895715 750.572578901226 750.471006262524 665.954573018229 662.421047663209 651.46556256004 651.417940134413 581.258189162714 568.105209810117 559.053702306466 558.493005833015 552.043762746846 550.131096080954 539.733572957862 530.698033203026 499.661483520139 494.928560833195 435.233706072008 427.70227707436 408.128705863823 399.856401759966 370.069343073402 369.845327696313 365.049514250289 363.639416679443 354.648012601512 346.609483937092 341.98916517756 339.373361078069 318.363717394716 276.219249213429 263.1610538422 257.409506256945 230.539454204132 229.778324183075 203.971681289205 197.504237163938 193.879371544746 189.866421885131 189.815806977935 167.944134441876 154.923566511974 146.21953885867]
    is_golden = [9.6662565e-5 0.0 0.001807884262 0.0 0.002166250472 0.0 0.003835132693 0.0 0.0 0.013550066967 0.018327409336 0.0 0.001319003523 0.0 0.0 0.006677115951 0.0 0.015583193495 0.028006929969 0.0 0.007782414773 0.0 0.028892276057 0.0 0.0 0.009748743725 0.0 0.007788123595 0.007278285719 0.0 0.0 0.00116773942 0.000190731866 0.000235139044 0.0 0.0 0.002193317995 0.0 0.008315770402 0.0 0.016520647062 0.0 0.0 0.014989172648 0.08363192972 0.001902190356 0.0 0.0 0.0 0.017111545439 0.007341641476 0.0 0.0 0.151419922476 0.094666992635 0.0 0.214639628105 0.0 0.0 0.072514184691 0.087279329707 0.0 0.0 0.23891458934 0.0 0.401546184942 0.0 0.050427766439 0.04544662783 0.0 0.0 0.078427781361 0.192350495231 0.0 0.002594420277 0.0 0.037800182163 0.081921558936 0.0 0.0]

    @test isapprox(disps, disps_golden; atol=1e-9)
    @test isapprox(is / (natoms * g^2), is_golden; atol=1e-9)

    # Test first 50 components of :full
    formula = intensity_formula(swt, :full, kernel=delta_function_kernel, formfactors=[FormFactor("Fe2")])
    _, is_full = intensities_bands(swt, ks, formula)
    is_full_flattened = reinterpret(reshape, ComplexF64, is_full)[1:50]

    # println(round.(is_full_flattened / (natoms * g^2); digits=12))
    is_full_golden = ComplexF64[0.000768755803 + 0.0im, 0.000453313199 - 4.8935388e-5im, 0.000468535469 + 8.5812792e-5im, 0.000453313199 + 4.8935388e-5im, 0.000270420761 + 0.0im, 0.000270819458 + 8.0426107e-5im, 0.000468535469 - 8.5812792e-5im, 0.000270819458 - 8.0426107e-5im, 0.000295138352 + 0.0im, 0.0 + 0.0im, 0.0 - 0.0im, 0.0 - 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 - 0.0im, 0.0 + 0.0im, 0.00021794048 + 0.0im, -0.000114399504 - 0.000211293781im, -0.000126018935 + 0.000199895685im, -0.000114399504 + 0.000211293781im, 0.000264899428 + 0.0im, -0.000127650501 - 0.000227103219im, -0.000126018935 - 0.000199895685im, -0.000127650501 + 0.000227103219im, 0.000256212416 + 0.0im, 0.0 + 0.0im, -0.0 - 0.0im, -0.0 + 0.0im, -0.0 + 0.0im, 0.0 + 0.0im, -0.0 - 0.0im, -0.0 - 0.0im, -0.0 + 0.0im, 0.0 + 0.0im, 7.5017149e-5 + 0.0im, 0.000206625047 + 0.000158601355im, 0.000240055384 - 0.000110167141im, 0.000206625047 - 0.000158601355im, 0.000904437193 + 0.0im, 0.000428286033 - 0.000810966567im, 0.000240055384 + 0.000110167141im, 0.000428286033 + 0.000810966567im, 0.000929965845 + 0.0im, 0.0 + 0.0im, 0.0 - 0.0im, -0.0 - 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]
    
    @test isapprox(is_full_flattened / (natoms * g^2), is_full_golden; atol=1e-9)
end


@testitem "Lanczos Bounds" begin
    using LinearAlgebra, Random

    Random.seed!(100)
    A = randn(ComplexF64, 500, 500)
    A = 0.5(A' + A)
    lo, hi = Sunny.eigbounds(A, 20)
    vals = eigvals(A)

    @test (abs(lo/vals[1] - 1) < 0.025) && (abs(hi/vals[end] - 1) < 0.025)
end


@testitem "Single Ion" begin
    # Tetragonal crystal
    a = 1.0
    c = 1.5
    latvecs = lattice_vectors(a, a, c, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)

    # System
    J, J′, D = 1.0, 0.1, 5.0
    infos = [SpinInfo(1, S=1, g=2)]
    sys = System(cryst, (1, 1, 1), infos, :SUN; seed=0)
    set_exchange!(sys, J,  Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 1, [0, 0, 1]))
    set_onsite_coupling!(sys, S -> D * S[3]^2, 1)

    # Reshape to sheared supercell and minimize energy
    A = [1 1 1; -1 1 0; 0 0 1]
    sys = reshape_supercell(sys, A)
    randomize_spins!(sys)
    @test minimize_energy!(sys) > 0

    k = rand(Float64, 3)
    swt = SpinWaveTheory(sys)
    ωk_num = dispersion(swt, [k])[1, :]

    function single_ion_analytical_disp(k)
        γkxy = cos(2π*k[1]) + cos(2π*k[2])
        γkz  = cos(2π*k[3])
        x = 1/2 - D/(8*(2*J+J′))
        Ak₊ = -8 * (x-1) * x * (2*J+J′) - (x-1) * D + 2 * (2*x-1) * (J *γkxy + J′*γkz)
        Bk₊ = -2 * (J * γkxy + J′ * γkz)
        Ak₋ = -16 * (x-1) * x * (2*J+J′) - (2*x-1) * D - 2 * (1-2*x)^2*(J*γkxy + J′*γkz)
        Bk₋ = 2 * (1-2*x)^2 * (J*γkxy + J′*γkz)
        ωk₊ = √(Ak₊^2-Bk₊^2)
        ωk₋ = √(Ak₋^2-Bk₋^2)
        return ωk₊, ωk₋
    end
    ωk1, ωk2 = single_ion_analytical_disp(k)
    ωk3, ωk4 = single_ion_analytical_disp(k + [0.5, 0.5, 0.5])
    ωk_ref = sort([ωk1, ωk2, ωk3, ωk4]; rev=true)

    @test ωk_num ≈ ωk_ref
end


@testitem "Intensities" begin
    using LinearAlgebra

    # Crystal
    a = 8.289
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [[0, 0, 0]]
    fcc = Crystal(latvecs, positions, 225)

    units = Units(:meV)
    S = 5/2
    g = 2
    J = 22.06 * units.K
    K = 0.15  * units.K
    C = J + K
    J₁ = diagm([J, J, C])
    D = 25/24

    dims = (1, 1, 1)
    infos = [SpinInfo(1; S, g)]

    function compute(mode)
        sys = System(fcc, dims, infos, mode)
        set_exchange!(sys, J₁, Bond(1, 2, [0, 0, 0]))
        set_onsite_coupling!(sys, S -> D * (S[1]^4 + S[2]^4 + S[3]^4), 1)
        set_dipole!(sys, (1, 1, 1), position_to_site(sys, (0, 0, 0)))
        set_dipole!(sys, (1, -1, -1), position_to_site(sys, (1/2, 1/2, 0)))
        set_dipole!(sys, (-1, -1, 1), position_to_site(sys, (1/2, 0, 1/2)))
        set_dipole!(sys, (-1, 1, -1), position_to_site(sys, (0, 1/2, 1/2)))
        swt = SpinWaveTheory(sys)
        k = [0.8, 0.6, 0.1]
        _, Sαβs =  Sunny.dssf(swt, [k])

        sunny_trace = [real(tr(Sαβs[1,a])) for a in axes(Sαβs)[2]]
        sunny_trace = filter(x -> abs(x) > 1e-12, sunny_trace)

        return sunny_trace
    end

    reference = [18.78918915723918, 19.679676835786527, 16.76890645406461]
    @test compute(:SUN) ≈ reference
    @test compute(:dipole) ≈ reference
end


@testitem "Scalar biquadratic" begin
    # Cubic crystal
    a = 2.0
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    
    function test_biquad(mode, k, S)
        # System
        dims = (2, 2, 2)
        infos = [SpinInfo(1; S, g=2)]
        sys = System(cryst, dims, infos, mode)        
        α = -0.4π
        J = 1.0
        JL, JQ = J * cos(α), J * sin(α) / S^2
        set_pair_coupling!(sys, (Si, Sj) -> Si'*JL*Sj + JQ*(Si'*Sj)^2, Bond(1, 1, [1, 0, 0]))

        # Initialize Néel order
        sys = reshape_supercell(sys, [1 1 1; -1 1 0; 0 0 1])
        set_dipole!(sys, ( 1, 0, 0), position_to_site(sys, (0, 0, 0)))
        set_dipole!(sys, (-1, 0, 0), position_to_site(sys, (0, 1, 0)))

        # Numerical result
        swt = SpinWaveTheory(sys)
        disp = dispersion(swt, [k])

        # Analytical result
        γk = 2 * (cos(2π*k[1]) + cos(2π*k[2]) + cos(2π*k[3]))
        disp_ref = J * (S*cos(α) - (2*S-2+1/S) * sin(α)) * √(36 - γk^2)
        
        @test disp[end-1] ≈ disp[end] ≈ disp_ref
    end

    k = [0.12, 0.23, 0.34]
    for mode in (:SUN, :dipole), S in (1, 3/2)
        test_biquad(mode, k, S)
    end
end


@testitem "General biquadratic" begin
    using LinearAlgebra

    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0], [0.4,0,0]]; types=["A", "B"])
    
    sys = System(cryst, (1,1,1), [SpinInfo(1; S=1, g=2), SpinInfo(2; S=2, g=2)], :dipole)
    set_pair_coupling!(sys, (S1, S2) -> +(S1'*diagm([2,-1,-1])*S1)*(S2'*diagm([2,-1,-1])*S2), Bond(1, 2, [0,0,0]))
    
    θ = randn()
    set_dipole!(sys, [1, 0, 0], (1, 1, 1, 1))
    set_dipole!(sys, [0, cos(θ), sin(θ)], (1, 1, 1, 2))
    energy(sys)
    @test energy(sys) ≈ -3
    
    swt = SpinWaveTheory(sys; apply_g=false)
    formula = intensity_formula(swt, :trace; kernel=delta_function_kernel)
    disp, intens = intensities_bands(swt, [[0,0,0]], formula)
    @test disp[1] ≈ 9
    @test intens[1] ≈ 1    
end


@testitem "Canted AFM" begin

    function test_canted_afm(S)
        J, D, h = 1.0, 0.54, 0.76
        a = 1
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)
        q = [0.12, 0.23, 0.34]
        
        sys = System(cryst, (1, 1, 1), [SpinInfo(1; S, g=-1)], :dipole)
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
        set_field!(sys, [0, 0, h])

        # Numerical
        sys_swt_dip = reshape_supercell(sys, [1 -1 0; 1 1 0; 0 0 1])
        c₂ = 1 - 1/(2S)
        θ = acos(h / (2S*(4J+D*c₂)))
        set_dipole!(sys_swt_dip, ( sin(θ), 0, cos(θ)), position_to_site(sys_swt_dip, (0,0,0)))
        set_dipole!(sys_swt_dip, (-sin(θ), 0, cos(θ)), position_to_site(sys_swt_dip, (1,0,0)))
        swt_dip = SpinWaveTheory(sys_swt_dip)
        ϵq_num = dispersion(swt_dip, [q])[1,:]

        # Analytical
        c₂ = 1 - 1/(2S)
        θ = acos(h / (2S*(4J+c₂*D)))
        Jq = 2J*(cos(2π*q[1])+cos(2π*q[2]))
        ωq₊ = real(√Complex(4J*S*(4J*S+2D*S*c₂*sin(θ)^2) + cos(2θ)*(Jq*S)^2 + 2S*Jq*(4J*S*cos(θ)^2 + c₂*D*S*sin(θ)^2)))
        ωq₋ = real(√Complex(4J*S*(4J*S+2D*S*c₂*sin(θ)^2) + cos(2θ)*(Jq*S)^2 - 2S*Jq*(4J*S*cos(θ)^2 + c₂*D*S*sin(θ)^2)))
        ϵq_ana = [ωq₊, ωq₋]

        ϵq_num ≈ ϵq_ana
    end

    @test test_canted_afm(1)
    @test test_canted_afm(2)
end


@testitem "Local Stevens expansion" begin
    using LinearAlgebra
    a = 1
    latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
    positions = [[0, 0, 0]]
    # P1 point group for most general single-ion anisotropy
    cryst = Crystal(latvecs, positions, 1)

    dims = (1, 1, 1)
    S = 2
    sys_dip = System(cryst, dims, [SpinInfo(1; S, g=-1)], :dipole)
    sys_SUN = System(cryst, dims, [SpinInfo(1; S, g=-1)], :SUN)

    # The strengths of single-ion anisotropy (must be negative to favor the dipolar ordering under consideration)
    Ds = -rand(3)
    h  = 0.1*rand()
    M = normalize(rand(3))
    SM = M' * spin_matrices(S)
    aniso = Ds[1]*SM^2 + Ds[2]*SM^4 + Ds[3]*SM^6

    set_onsite_coupling!(sys_dip, aniso, 1)
    set_onsite_coupling!(sys_SUN, aniso, 1)

    set_field!(sys_dip, h*M)
    set_field!(sys_SUN, h*M)

    set_dipole!(sys_dip, M, (1,1,1,1))
    set_dipole!(sys_SUN, M, (1,1,1,1))

    energy(sys_dip)
    energy(sys_SUN)

    swt_dip = SpinWaveTheory(sys_dip)
    swt_SUN = SpinWaveTheory(sys_SUN)

    q = rand(3)
    disp_dip = dispersion(swt_dip, [q])
    disp_SUN = dispersion(swt_SUN, [q])

    @test only(disp_dip) ≈ disp_SUN[end-1]
end


@testitem "Intensities interface" begin
    sys = System(Sunny.diamond_crystal(),(1,1,1),[SpinInfo(1,S=1/2,g=2)],:SUN;seed = 0)
    randomize_spins!(sys)

    swt = SpinWaveTheory(sys)
    
    # Just testing that nothing throws errors
    # TODO: accuracy check
    path, _ = reciprocal_space_path(Sunny.diamond_crystal(),[[0.,0.,0.],[0.5,0.5,0.]],50)
    energies = collect(0:0.1:5)

    # Bands
    formula = intensity_formula(swt, :perp, kernel=delta_function_kernel)
    intensities_bands(swt, path, formula)
    @test_throws "without a finite-width kernel" intensities_broadened(swt,path,energies,formula)

    # Broadened
    formula = intensity_formula(swt, :perp, kernel=lorentzian(fwhm=0.1))
    intensities_broadened(swt, path, energies, formula)
    @test_throws "broadening kernel" intensities_bands(swt,path,formula)

    formula = intensity_formula(swt, :perp, kernel=lorentzian(fwhm=0.1))
    intensities_broadened(swt, path, energies, formula)
    @test_throws "broadening kernel" intensities_bands(swt,path,formula)

    # Full
    formula = intensity_formula(swt, :full, kernel=lorentzian(fwhm=0.1))
    intensities_broadened(swt, path, energies, formula)
end


@testitem "Dipole-dipole" begin
    latvecs = lattice_vectors(10, 10, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])

    for mode in (:dipole, :SUN)
        sys = System(cryst, (1,1,1), [SpinInfo(1; S=1, g=1)], mode)
        enable_dipole_dipole!(sys, 1.0)

        polarize_spins!(sys, (0,0,1))
        @test energy_per_site(sys) ≈ -0.1913132980155851
        
        swt = SpinWaveTheory(sys)
        formula = intensity_formula(swt, :perp; kernel=delta_function_kernel)
        
        qpoints = [[0, 0, 0], [0, 0, 1/2], [0, 1/2, 1/2], [0, 0, 0]]
        disps, is = intensities_bands(swt, qpoints, formula)
        disps_ref = [0.5689399140467553, 0.23914164251944922, 0.23914164251948083, 0.5689399140467553]
        @test isapprox(disps[:,end], disps_ref; atol=1e-7)
        @test is[:,end] ≈ [1, 1, 201/202, 1]
    end

    begin
        cryst = Sunny.bcc_crystal()
        sys = System(cryst, (1, 1, 1), [SpinInfo(1, S=1, g=2)], :dipole, seed=2)
        enable_dipole_dipole!(sys, Units(:meV).vacuum_permeability)
        polarize_spins!(sys, (1,2,3)) # arbitrary direction
        
        R = hcat([1,1,-1], [-1,1,1], [1,-1,1]) / 2
        sys_reshape = reshape_supercell(sys, R)
        @test energy_per_site(sys_reshape) ≈ energy_per_site(sys) ≈ -0.89944235377
        
        swt1 = SpinWaveTheory(sys, energy_ϵ=1e-8)
        swt2 = SpinWaveTheory(sys_reshape, energy_ϵ=1e-8)
        path = [[0.5, -0.1, 0.3]]
        disp1 = dispersion(swt1, path)
        disp2 = dispersion(swt2, path)
        
        @test disp1 ≈ [1.3236778213351734 0.9206655419611791]
        @test disp2 ≈ [0.9206655419611772]        
    end
end


@testitem "SW15-Langasite" begin
    # Ba3NbFe3Si2O14
    a = b = 8.539
    c = 5.2414
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    crystal = Crystal(latvecs, [[0.24964,0,0.5]], 150)
    latsize = (1,1,7)
    sys = System(crystal, latsize, [SpinInfo(1; S=5/2, g=2)], :dipole; seed=5)
    set_exchange!(sys, 0.85,  Bond(3, 2, [1,1,0]))   # J1
    set_exchange!(sys, 0.24,  Bond(1, 3, [0,0,0]))   # J2
    set_exchange!(sys, 0.053, Bond(2, 3, [-1,-1,1])) # J3
    set_exchange!(sys, 0.017, Bond(1, 1, [0,0,1]))   # J4
    set_exchange!(sys, 0.24,  Bond(3, 2, [1,1,1]))   # J5
    
    for i in 1:3
        θ = -2π*(i-1)/3
        set_spiral_order_on_sublattice!(sys, i; k=[0,0,1/7], axis=[0,0,1], S0=[cos(θ),sin(θ),0])
    end

    swt = SpinWaveTheory(sys, apply_g=false)
    formula = intensity_formula(swt, :full; kernel=delta_function_kernel)
    q = [0.41568,0.56382,0.76414]
    disp, intensities = intensities_bands(swt, [q], formula)
    SpinW_energies = [2.6267,2.6541,2.8177,2.8767,3.2458,3.3172,3.4727,3.7767,3.8202,3.8284,3.8749,3.9095,3.9422,3.9730,4.0113,4.0794,4.2785,4.4605,4.6736,4.7564,4.7865]
    SpinW_intensities = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2999830079, -0.2999830079im, 0,0.2999830079im, 0.2999830079, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3591387785, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5954018134, -0.5954018134im, 0,0.5954018134im, 0.5954018134, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.3708506016,1.3708506016im, 0, -1.3708506016im, 1.3708506016, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0511743697, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0734875342, 0.0 + 0.0734875342im, 0, 0.0 - 0.0734875342im, 0.0734875342, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0577275935, -0.0577275935im, 0,0.0577275935im, 0.0577275935, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.1733740706, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0338873034,0.0338873034im, 0, -0.0338873034im, 0.0338873034, 0, 0, 0, 0]
    
    @test isapprox(disp[:], reverse(SpinW_energies); atol=1e-3)
    
    intensities_reshaped = reinterpret(reshape, ComplexF64, intensities)[:]
    @test isapprox(SpinW_intensities, intensities_reshaped; atol=1e-7)
end


@testitem "3Q Pyrochlore" begin
    tol = 1e-7
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, (7,7,1), [SpinInfo(1; S=1, g=-1)], :dipole; seed=0)

    J₁ = -1
    J₃ = 1.6234898018587323
    set_exchange!(sys, J₁, Bond(1, 1, [1,0,0]))
    set_exchange!(sys, J₃, Bond(1, 1, [2,0,0]))

    D = 0.066*J₁
    set_onsite_coupling!(sys, S -> D*S[3]^2, 1)

    h = 2.0*J₃
    set_field!(sys, [0, 0, h])
    
    dipole_array = zeros(Float64,3,7,7)
    dipole_array[:,:,1]= [0.856769   0.811926  -0.485645   -0.950288   -0.429787  -0.0800935  0.217557;
                          0.35333   -0.2906     0.0685251  -0.297433   -0.854291  -0.004611   0.97169;
                          -0.375638   0.506288   0.871466    0.0921234   0.292353   0.996777   0.0921244]
    dipole_array[:,:,2]= [-0.2422     0.428786   0.365354  -0.0893441   0.0440378   4.84897e-7  -0.524944;
                          0.368806  -0.652928  -0.556337   0.136048   -0.0670579  -7.39508e-7   0.799351;
                         -0.897397  -0.624361   0.746328   0.986665    0.996777    1.0          0.292356]
    dipole_array[:,:,3]= [-0.978133   -0.664683  -0.0559526  0.130092  0.650549   0.954731   0.0360572;
                          0.186465   -0.645827  -0.860547   0.472888  0.753856   0.0549402  0.0716667;
                          0.0921244  -0.375638   0.506288   0.871466  0.0921234  0.292353   0.996777]
    dipole_array[:,:,4]= [-0.0731477  -0.474581  -0.179152  -0.245875  -0.226963   0.440493    0.73273;
                          -0.145399   -0.123781  -0.387346   0.348821   0.898546   0.0253494  -0.674256;
                          0.986665    0.871465   0.90436    0.904361  -0.375636  -0.897398    0.0921208]
    dipole_array[:,:,5]= [0.183479  0.299124  0.773231  0.425028   -0.717278  -0.779847   -0.122394;
                          -0.454847  0.594572  0.381818  0.0385228   0.478729  -0.0448744  -0.918648;
                          0.871464  0.74633   0.506289  0.90436     0.506291  -0.624359   -0.375639]
    dipole_array[:,:,6]= [-0.20424   -0.154295   0.351059   0.891644   0.34449   -0.664477   -0.657631;
                          -0.374726   0.848447   0.697806  -0.252717  -0.349109  -0.0382366  -0.557847;
                          0.90436    0.506291  -0.624359  -0.375639   0.871464   0.74633     0.506289]
    dipole_array[:,:,7]= [0.426641   -0.2224    -0.734377  -0.198296   0.327583   0.162493    0.302168;
                          0.0104859   0.364238   0.565321  -0.394152  -0.940321   0.00934933  0.38632;
                          0.90436     0.904361  -0.375636  -0.897398   0.0921208  0.986665    0.871465]
    for x in 1:7, y in 1:7
        set_dipole!(sys, dipole_array[:,x,y], (x,y,1,1))
    end
    
    q1 = [0.0391,0.0415,0.5530]
    q2 = [0.2360,0.7492,0.9596]
    q3 = [0.1131,0.7654,0.2810]
    q = [q1,q2,q3]
    swt = SpinWaveTheory(sys);
    formfactors = [FormFactor("Cr4")]
    formula = intensity_formula(swt, :full; kernel=delta_function_kernel, formfactors)
    disp, intensities = intensities_bands(swt, q, formula)
    disp_inds = [27, 119, 60, 126, 42, 46, 15, 77, 132, 52]
    int_inds = [25, 147, 99, 121, 43, 140, 142, 21, 93, 49]
    calc_disp = disp[disp_inds]
    calc_int = intensities[int_inds]
    disp_ref = [8.464621970889235,2.965829202488746,6.539681848582543,2.524276472373584,7.536305768861917,7.21157510322424,9.267100207705882,5.603801899767303,2.2012141464553636,6.933800585478572]
    int_ref = [[1.2252940671236579e-5 + 0.0im 2.274979028773366e-5 - 2.9311941591805888e-6im 1.653212548779866e-5 + 1.1781816967947205e-5im; 2.274979028773366e-5 + 2.9311941591805888e-6im 4.294029257562512e-5 + 0.0im 2.7876377103396445e-5 + 2.582994102312879e-5im; 1.653212548779866e-5 - 1.1781816967947205e-5im 2.7876377103396445e-5 - 2.582994102312879e-5im 3.363456946935e-5 + 0.0im], [0.0023849068155162296 + 0.0im 0.00012123326788175223 + 0.0005597432044950685im 2.9024438353305692e-5 + 0.0019488022511109527im; 0.00012123326788175223 - 0.0005597432044950685im 0.00013753575531155366 + 0.0im 0.0004588638588955403 + 9.225242336242961e-5im; 2.9024438353305692e-5 - 0.0019488022511109527im 0.0004588638588955403 - 9.225242336242961e-5im 0.0015927970884407876 + 0.0im], [3.343948808906739e-5 + 0.0im -0.00020469198105989745 - 4.273264350562476e-5im 0.00014834426130997533 + 1.9072039037346926e-5im; -0.00020469198105989745 + 4.273264350562476e-5im 0.0013075823952430512 + 0.0im -0.0009324269345411454 + 7.282554609603465e-5im; 0.00014834426130997533 - 1.9072039037346926e-5im -0.0009324269345411454 - 7.282554609603465e-5im 0.0006689624696724299 + 0.0im], [0.000371215929258512 + 0.0im 0.00015333075423969357 + 0.0003656727495761781im -0.0005009238108982899 + 0.0003075677944476983im; 0.00015333075423969357 - 0.0003656727495761781im 0.0004235456174856762 + 0.0im 9.606844029120137e-5 + 0.0006204846586781491im; -0.0005009238108982899 - 0.0003075677944476983im 9.606844029120137e-5 - 0.0006204846586781491im 0.0009307860608149384 + 0.0im], [2.9077299015972345e-6 + 0.0im 7.621667422789561e-6 + 2.799505642406002e-6im -2.266990384026427e-6 - 2.148518183351213e-6im; 7.621667422789561e-6 - 2.799505642406002e-6im 2.2673029606106315e-5 + 0.0im -8.010728755108931e-6 - 3.4490269087400155e-6im; -2.266990384026427e-6 + 2.148518183351213e-6im -8.010728755108931e-6 + 3.4490269087400155e-6im 3.3549800413375375e-6 + 0.0im], [4.79789862126743e-5 + 0.0im -2.0507170782740867e-5 - 1.8427040499771927e-5im 1.1851577631841944e-7 + 3.659464580060244e-5im; -2.0507170782740867e-5 + 1.8427040499771927e-5im 1.5842349642893138e-5 + 0.0im -1.4105371891520816e-5 - 1.5595760044615497e-5im; 1.1851577631841944e-7 - 3.659464580060244e-5im -1.4105371891520816e-5 + 1.5595760044615497e-5im 2.7911847518508492e-5 + 0.0im], [0.004511188103971374 + 0.0im -0.0013855170428121722 - 0.003178210621811203im 0.0012245622270242976 + 0.0004875348349735645im; -0.0013855170428121722 + 0.003178210621811203im 0.0026646373317783553 + 0.0im -0.0007195754536780329 + 0.0007129893012763896im; 0.0012245622270242976 - 0.0004875348349735645im -0.0007195754536780329 - 0.0007129893012763896im 0.0003850965251566534 + 0.0im], [1.1287527510791432e-5 + 0.0im 9.784635886704085e-6 + 4.074063598919808e-5im 1.177917483079537e-5 + 1.7619012333705406e-5im; 9.784635886704085e-6 - 4.074063598919808e-5im 0.00015552905794128413 + 0.0im 7.380400215595438e-5 - 2.7242055744441945e-5im; 1.177917483079537e-5 - 1.7619012333705406e-5im 7.380400215595438e-5 + 2.7242055744441945e-5im 3.979423792148188e-5 + 0.0im], [9.457308666169535e-5 + 0.0im -1.7915569710789528e-6 + 2.743855511588096e-5im 2.3814745247182995e-5 - 0.00010529497917934995im; -1.7915569710789528e-6 - 2.743855511588096e-5im 7.994705575514363e-6 + 0.0im -3.100044278968906e-5 - 4.914720059098736e-6im; 2.3814745247182995e-5 + 0.00010529497917934995im -3.100044278968906e-5 + 4.914720059098736e-6im 0.00012322929432616498 + 0.0im], [2.4114075174233513e-6 + 0.0im -2.4607194255187506e-6 - 1.9339332806224738e-6im 2.7253619853740792e-6 + 6.205353429165704e-7im; -2.4607194255187506e-6 + 1.9339332806224738e-6im 4.062041755385664e-6 + 0.0im -3.278759427148752e-6 + 1.5524978029108779e-6im; 2.7253619853740792e-6 - 6.205353429165704e-7im -3.278759427148752e-6 - 1.5524978029108779e-6im 3.2398762990830975e-6 + 0.0im]]
    @test isapprox(calc_disp, disp_ref; atol=1e-6)
    @test isapprox(calc_int, int_ref; atol=1e-7)
end


@testitem "Invariance to reshaping" begin    
    # Diamond-cubic with antiferromagnetic exchange
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]], 227, setting="1")
    S = 3/2
    sys = System(cryst, (1, 1, 1), [SpinInfo(1; S, g=2)], :dipole; seed=0)
    set_exchange!(sys, 1.0, Bond(1, 3, [0,0,0]))
    randomize_spins!(sys)
    minimize_energy!(sys)
    @test energy_per_site(sys) ≈ -2S^2
    
    # Reshaped system
    shape = [0 1 1; 1 0 1; 1 1 0] / 2
    sys_prim = reshape_supercell(sys, shape)
    @test energy_per_site(sys_prim) ≈ -2S^2
    
    # Both systems should produce the same intensities
    swt1 = SpinWaveTheory(sys_prim)
    swt2 = SpinWaveTheory(sys)
    kernel = lorentzian(fwhm=0.8)
    formfactors = [FormFactor("Co2")]
    formula1 = intensity_formula(swt1, :perp; kernel, formfactors)
    formula2 = intensity_formula(swt2, :perp; kernel, formfactors)
    q = randn(3)
    energies = collect(0:0.01:6)
    is1 = intensities_broadened(swt1, [q], energies, formula1)
    is2 = intensities_broadened(swt2, [q], energies, formula2)
    @test is1 ≈ is2        
end


@testitem "Invariance to spin rotation" begin
    using LinearAlgebra

    function build_system(R, D1, D2, J, K1, K2, h, g)
        latvecs = lattice_vectors(1, 1, 1, 92, 93, 94)
        cryst = Crystal(latvecs, [[0,0,0], [0.4,0,0]]; types=["A", "B"])
        infos = [SpinInfo(1; S=1, g=2), SpinInfo(2; S=2, g=R*g*R')]
        sys = System(cryst, (1,1,1), infos, :dipole)

        set_onsite_coupling!(sys, S -> S'*R*(D1+D1')*R'*S, 1)
        set_onsite_coupling!(sys, S -> S'*R*(D2+D2')*R'*S, 2)

        K1 = Sunny.tracelesspart(K1 + K1')
        K2 = Sunny.tracelesspart(K2 + K2')

        set_pair_coupling!(sys, (S1, S2) -> S1'*R*J*R'*S2 + (S1'*R*K1*R'*S1)*(S2'*R*K2*R'*S2), Bond(1, 2, [0,0,0]))
        set_field!(sys, R*h)

        return sys
    end

    g = randn(3,3)
    D1 = randn(3,3)
    D2 = randn(3,3)
    J = randn(3,3)
    K1 = randn(3,3)
    K2 = randn(3,3)
    h = randn(3)
    R1 = Sunny.Mat3(I)
    R2 = Sunny.axis_angle_to_matrix(randn(3), 0.46)
    sys1 = build_system(R1, D1, D2, J, K1, K2, h, g)
    sys2 = build_system(R2, D1, D2, J, K1, K2, h, g)

    randomize_spins!(sys1)
    minimize_energy!(sys1)

    for site in eachsite(sys1)
        sys2.dipoles[site] = R2 * R1' * sys1.dipoles[site]
    end
    @assert energy_per_site(sys1) ≈ energy_per_site(sys2)

    swt = SpinWaveTheory(sys1)
    formula = intensity_formula(swt, :trace; kernel=delta_function_kernel)
    res1 = intensities_bands(swt, [[0,0,0]], formula)

    swt = SpinWaveTheory(sys2)
    formula = intensity_formula(swt, :trace; kernel=delta_function_kernel)
    res2 = intensities_bands(swt, [[0,0,0]], formula)

    @assert res1[1] ≈ res2[1]
    @assert res1[2] ≈ res2[2]
end


@testitem "Generalized interaction consistency" begin
    using LinearAlgebra

    # Construct LSWT Hamiltonian for Heisenberg ferromagnet at a particular wave
    # vector k
    function make_fm_lswt_hamiltonian(k; extract_parts=true)
        k = Sunny.Vec3(k)
        dims = (1, 1, 1)
        cryst = Crystal(I(3), [[0,0,0]]) 
        sys = System(cryst, dims, [SpinInfo(1; S=1, g=1)], :SUN)
        J = -1.0
        S = spin_matrices(1)
        S1, S2 = Sunny.to_product_space(S, S)
        exchange = J*(S1'*S2)
        set_pair_coupling!(sys, exchange, Bond(1,1,[1,0,0]); extract_parts)

        swt = SpinWaveTheory(sys)
        nmodes = (swt.sys.Ns[1] - 1) * Sunny.natoms(swt.sys.crystal)
        H = zeros(ComplexF64, 2nmodes, 2nmodes)
        Sunny.swt_hamiltonian_SUN!(H, swt, k)

        return H
    end

    # Test whether results are identical when using both conventional approach and
    # generalized interactions (tensor decomposition)
    k = [0.23, 0, 0]
    H_conventional = make_fm_lswt_hamiltonian(k; extract_parts=true)
    H_generalized = make_fm_lswt_hamiltonian(k; extract_parts=false)

    @test H_conventional ≈ H_generalized
end


@testitem "Equivalence of dense and sparse Hamiltonian constructions" begin
    import LinearAlgebra: diagm

    function onehot(i, n)
        out = zeros(ComplexF64, 1, n)
        out[1, i] = 1.0
        return out
    end

    # Build System and SpinWaveTheory with exchange, field and single-site anisotropy
    function simple_swt(mode)
        dims = (8, 1, 1)
        cryst = Crystal(diagm([1, 1, 2.0]), [[0,0,0]], "P1")
        sys = System(cryst, dims, [SpinInfo(1; S=1, g=1)], mode)

        K1 = diagm([2, -1, -1])
        K2 = diagm([-1, -1, 2])
        set_pair_coupling!(sys, (Si, Sj) -> -Si'*Sj + (Si'*K1*Si)*(Sj'*K2*Sj), Bond(1,1,[1,0,0]); extract_parts=true)
        set_onsite_coupling!(sys, S -> S[3]^2, 1)
        set_field!(sys, [0,0,0.1])

        randomize_spins!(sys)
        minimize_energy!(sys; maxiters=1_000)

        return SpinWaveTheory(sys)
    end

    # Construct dipole Hamiltonian standard way
    swt = simple_swt(:dipole)
    H1 = zeros(ComplexF64, 16, 16)
    q = Sunny.Vec3([0.5,0,0])
    Sunny.swt_hamiltonian_dipole!(H1, swt, q)

    # Construct dipole Hamiltonian from sparse matrix-vector multiplies
    H2 = zero(H1)
    L = Sunny.natoms(swt.sys.crystal)
    for i in 1:2L
        H2[:, i] = Sunny.multiply_by_hamiltonian_dipole(onehot(i, 2L), swt, [q])
    end

    @test isapprox(H1, H2; atol=1e-12)

    # Construct SU(N) Hamiltonian standard way
    swt = simple_swt(:SUN)
    N = swt.sys.Ns[1]
    L = (N-1) * Sunny.natoms(swt.sys.crystal)
    H1 = zeros(ComplexF64, 2L, 2L)
    Sunny.swt_hamiltonian_SUN!(H1, swt, q)

    # Construct SU(N) Hamiltonian from sparse matrix-vector multiplies
    H2 = zero(H1)
    for i in 1:2L
        H2[:, i] = Sunny.multiply_by_hamiltonian_SUN(onehot(i, 2L), swt, [q])
    end

    @test isapprox(H1, H2; atol=1e-12)
end


@testitem "Spin ladder intensities reference test" begin
    using LinearAlgebra

    function dimer_model(; J=1.0, J′=0.0, h=0.0, dims=(2,1,1), fast=false)
        cryst = Crystal(I(3), [[0,0,0]], 1) 
        sys = System(cryst, dims, [SpinInfo(1; S=3/2, g=1)], :SUN)

        S = spin_matrices(1/2)
        S1, S2 = Sunny.to_product_space(S, S)
        onsite = J*(S1' * S2 - 0.25I(4)) - h*(S1[3] + S2[3])
        set_onsite_coupling!(sys, onsite, 1)

        S1i, S1j = Sunny.to_product_space(S1, S1)
        S2i, S2j = Sunny.to_product_space(S2, S2)
        scalar, bilin, biquad, tensordec = Sunny.decompose_general_coupling(J′*(S1i' * S1j + S2i' * S2j), 4, 4; extract_parts=fast)
        Sunny.set_pair_coupling_aux!(sys, scalar, Sunny.Mat3(I)*bilin, biquad, tensordec, Bond(1, 1, [-1,0,0]))  

        return sys, cryst
    end

    # Set up dimer model and observables (0 and π channels manually specified)
    sys, cryst = dimer_model(; J=1.0, J′=0.2, h=0.0, dims=(2,1,1), fast=false) 
    S = spin_matrices(1/2)
    S1, S2 = Sunny.to_product_space(S, S)
    observables = [
        1/√2*(S1[1] - S2[1]),
        1/√2*(S1[2] - S2[2]),
        1/√2*(S1[3] - S2[3]),
    ]

    # Set up SpinWaveTheory
    randomize_spins!(sys)
    minimize_energy!(sys)
    swt = SpinWaveTheory(sys; observables)

    qs = [[0,0,0], [0.5,0,0], [1,0,0]]
    broadened_formula = intensity_formula(swt, :trace; kernel=lorentzian(fwhm=0.3))
    energies = collect(0:0.1:5)
    is = intensities_broadened(swt, qs, energies, broadened_formula);
    is_ref = [0.042551643530614226 0.05061618774203612 0.06118972834822155 0.07541981401191966 0.09518339306332815 0.12371078181463013 0.1669136576105103 0.23644634156679015 0.3574143016120729 0.589319004913063 1.0795750388213625 2.0570916947398903 2.6569442073245604 1.6749367428248116 0.8709898932656854 0.4927038808971025 0.30849624153704924 0.20903618563778026 0.1502266820186879 0.11286975075833217 0.08777050885259896 0.07013929292833648 0.057300827672198455 0.047672208009393216 0.04027090049655778 0.0344619791000266 0.02982086993597283 0.02605519995256082 0.022958441606058772 0.020381422739929135 0.01821426297223774 0.016374601990583864 0.014799738769001156 0.013441266708211906 0.01226133976759786 0.01123002734145714 0.010323410070056573 0.009522188815697288 0.008810654796662341 0.008175917662895354 0.007607320304517221 0.007095990541228231 0.006634494316422953 0.006216564975067639 0.005836890143659928 0.00549094262866035 0.005174845247781923 0.0048852620341421574 0.004619310095626147 0.004374487768740016 0.00414861571474443; 0.14853118271457438 0.19360218524421224 0.2621797495216911 0.3732135013522173 0.5678610986270936 0.944407657666259 1.7450675448755357 3.2945578443571577 3.9947874291087895 2.418790606759786 1.2612861617795654 0.7201697967624341 0.45442367736983386 0.30970051154679945 0.22353509755861065 0.1685055354270234 0.1313752454162916 0.10520387618639207 0.0860938574888819 0.0717287252840356 0.060665219355115055 0.05196772774655005 0.045008912709399315 0.03935576182418519 0.03470177881993486 0.030825189142763776 0.027562374631854493 0.024790523018178495 0.02241601889071837 0.020366506674885078 0.018585357752204528 0.017027745220739847 0.015657814448345492 0.014446613655329999 0.013370560099471452 0.01241028925551664 0.011549781566933818 0.010775692876605122 0.010076836041215703 0.009443775967574439 0.008868510590520351 0.008344217576743833 0.007865051732026552 0.007425981842385972 0.007022658419610481 0.006651305841367626 0.006308633878300496 0.005991764727412878 0.005698172523177506 0.0054256329470820505 0.005172181054614204; 0.042551643530614205 0.0506161877420361 0.06118972834822153 0.07541981401191963 0.0951833930633281 0.12371078181463005 0.16691365761051016 0.23644634156678992 0.3574143016120724 0.5893190049130621 1.0795750388213603 2.0570916947398867 2.656944207324563 1.6749367428248163 0.8709898932656877 0.4927038808971036 0.30849624153704985 0.20903618563778062 0.15022668201868813 0.1128697507583323 0.08777050885259907 0.07013929292833657 0.05730082767219852 0.04767220800939327 0.04027090049655782 0.03446197910002663 0.02982086993597286 0.026055199952560844 0.02295844160605879 0.020381422739929156 0.018214262972237757 0.016374601990583878 0.014799738769001168 0.013441266708211913 0.01226133976759787 0.011230027341457147 0.010323410070056582 0.009522188815697295 0.008810654796662348 0.00817591766289536 0.007607320304517226 0.007095990541228234 0.006634494316422957 0.006216564975067644 0.005836890143659932 0.005490942628660353 0.005174845247781927 0.00488526203414216 0.00461931009562615 0.0043744877687400185 0.0041486157147444325]

    @test is ≈ is_ref
end


@testitem "LSWT correction to classical energy" begin
    J = 1
    S = 1
    δE_afm1_ref = 0.488056/(2S) * (-2*J*S^2)

    # The results are taken from Phys. Rev. B 102, 220405(R) (2020) for the AFM1
    # phase on the FCC lattice
    function correction(mode)
        a = 1
        latvecs = lattice_vectors(a, a, a, 90, 90, 90)
        positions = [[0, 0, 0]]
        fcc = Crystal(latvecs, positions, 225)
        sys_afm1 = System(fcc, (1, 1, 1), [SpinInfo(1; S, g=1)], mode)
        set_exchange!(sys_afm1, J, Bond(1, 2, [0, 0, 0]))
        set_dipole!(sys_afm1, (0, 0,  1), position_to_site(sys_afm1, (0, 0, 0)))
        set_dipole!(sys_afm1, (0, 0, -1), position_to_site(sys_afm1, (1/2, 1/2, 0)))
        set_dipole!(sys_afm1, (0, 0, -1), position_to_site(sys_afm1, (1/2, 0, 1/2)))
        set_dipole!(sys_afm1, (0, 0,  1), position_to_site(sys_afm1, (0, 1/2, 1/2)))
        swt_afm1 = SpinWaveTheory(sys_afm1)
        # Calculate at low accuracy for faster testing
        δE_afm1 = Sunny.energy_per_site_lswt_correction(swt_afm1; atol=5e-4)
        return isapprox(δE_afm1_ref, δE_afm1; atol=1e-3)
    end

    for mode in (:dipole, :SUN)
        @test correction(mode)
    end
end


@testitem "LSWT correction to the ordered moments (S maximized)" begin
    # Test example 1: The magnetization is maximized to `S`. Reference result
    # comes from Phys. Rev. B 79, 144416 (2009) Eq. (45) for the 120° order on
    # the triangular lattice.
    J = 1
    S = 1/2
    a = 1
    δS_ref = -0.261302

    function δS_triangular(mode)
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 120)
        cryst = Crystal(latvecs, [[0, 0, 0]])
        sys = System(cryst, (3, 3, 1), [SpinInfo(1, S=S, g=2)], mode)
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        set_spiral_order!(sys; k=[2/3, -1/3, 0], axis=[0, 0, 1], S0=[0, 1, 0])
        swt = SpinWaveTheory(sys)
        # Calculate first 3 digits for faster testing
        δS = Sunny.magnetization_lswt_correction(swt; atol=1e-3)[1]

        return isapprox(δS_ref, δS, atol=1e-3)
    end

    for mode in (:dipole, :SUN)
        @test δS_triangular(mode)
    end
end


@testitem "LSWT correction to the ordered moments (S not maximized)" begin
    using LinearAlgebra
    # Test example 2: The magnetization is smaller than `S` due to easy-plane
    # single-ion anisotropy The results are derived in the Supplemental
    # Information (Note 12) of Nature Comm. 12.1 (2021): 5331.
    a = b = 8.3193
    c = 5.3348
    lat_vecs = lattice_vectors(a, b, c, 90, 90, 90)
    types = ["Fe"]
    positions = [[0, 0, 0]]
    cryst = Crystal(lat_vecs, positions, 113; types)

    S = 1
    J₁  = 0.266
    J₁′ = 0.1J₁
    Δ = Δ′ = 1/3
    D = 1.42
    gab, gcc = 2.18, 1.93
    g = diagm([gab, gab, gcc])
    x = 1/2 - D/(8*(2J₁+J₁′))

    sys = System(cryst, (1, 1, 2), [SpinInfo(1; S, g)], :SUN, seed=0)
    set_exchange!(sys, diagm([J₁, J₁, J₁*Δ]),  Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, diagm([J₁′, J₁′, J₁′*Δ′]), Bond(1, 1, [0, 0, 1]))
    set_onsite_coupling!(sys, S -> D*S[3]^2, 1)

    randomize_spins!(sys)
    minimize_energy!(sys; maxiters=1000)
    sys_swt = reshape_supercell(sys, diagm([1,1,2]))
    minimize_energy!(sys_swt)
    swt = SpinWaveTheory(sys_swt)

    δS = Sunny.magnetization_lswt_correction(swt; atol=1e-2)[1]

    M_cl  = 2*√((1-x)*x)
    # Paper reported M_ref = 2.79, but actual result is closer to 2.78
    M_ref = 2.78
    @test isapprox(M_ref, (M_cl+δS)*√3*gab, atol=1e-2)
end
