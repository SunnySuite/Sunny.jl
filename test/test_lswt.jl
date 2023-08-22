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
    S = spin_operators(sys, 1)
    set_onsite_coupling!(sys, D * S[3]^2, 1)

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
    types = ["MzR1"]
    positions = [[0, 0, 0]]
    fcc = Crystal(latvecs, positions, 225; types)

    S = 5/2
    J = 22.06 * meV_per_K
    K = 0.15  * meV_per_K
    C = J + K
    J₁ = diagm([J, J, C])
    D_ST = 0.2
    # Undo Sunny-applied renormalization of quartic anisotropy
    D = D_ST / Sunny.anisotropy_renormalization(2S+1)[4]

    dims = (1, 1, 1)
    infos = [SpinInfo(1; S, g=2)]

    function compute(mode)
        sys = System(fcc, dims, infos, mode)
        set_exchange!(sys, J₁, Bond(1, 2, [0, 0, 0]))
        S = spin_operators(sys, 1)
        Λ = D * (S[1]^4 + S[2]^4 + S[3]^4)
        set_onsite_coupling!(sys, Λ, 1)
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

    reference = [1.1743243223274487, 1.229979802236658, 1.048056653379038]
    @test compute(:SUN) ≈ compute(:dipole) ≈ reference
end

@testitem "Biquadratic interactions" begin
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
        set_exchange!(sys, JL,  Bond(1, 1, [1, 0, 0]); biquad=JQ)

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

@testitem "Canted AFM" begin

    function test_canted_afm(S)
        J, D, h = 1.0, 0.54, 0.76
        a = 1
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)
        q = [0.12, 0.23, 0.34]
        
        dims = (2, 2, 1)
        sys = System(cryst, dims, [SpinInfo(1; S, g=1)], :dipole; units=Units.theory)
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        Sz = spin_operators(sys,1)[3]
        set_onsite_coupling!(sys, D*Sz^2, 1)
        set_external_field!(sys, [0, 0, h])

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

@testitem "Local stevens expansion" begin
    using LinearAlgebra
    a = 1
    latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
    positions = [[0, 0, 0]]
    # P1 point group for most general single-ion anisotropy
    cryst = Crystal(latvecs, positions, 1)

    dims = (1, 1, 1)
    S = 3
    sys_dip = System(cryst, dims, [SpinInfo(1; S, g=1)], :dipole)
    sys_SUN = System(cryst, dims, [SpinInfo(1; S, g=1)], :SUN)

    # The strengths of single-ion anisotropy (must be negative to favor the dipolar ordering under consideration)
    Ds = -rand(3)
    h  = rand()
    # Random magnetic moment
    M = normalize(rand(3))
    θ, ϕ = Sunny.dipole_to_angles(M)

    s_mat = Sunny.spin_matrices(N=2S+1)
    
    s̃ᶻ = M' * spin_operators(sys_dip,1)
    
    U_mat = exp(-1im * ϕ * s_mat[3]) * exp(-1im * θ * s_mat[2])
    hws = zeros(2S+1)
    hws[1] = 1.0
    Z = U_mat * hws

    aniso = Ds[1]*s̃ᶻ^2 + Ds[2]*s̃ᶻ^4 + Ds[3]*s̃ᶻ^6

    set_onsite_coupling!(sys_dip, aniso, 1)
    
    s̃ᶻ = M' * spin_operators(sys_SUN,1)
    aniso = Ds[1]*s̃ᶻ^2 + Ds[2]*s̃ᶻ^4 + Ds[3]*s̃ᶻ^6
    set_onsite_coupling!(sys_SUN, aniso, 1)

    set_external_field!(sys_dip, h*M)
    set_external_field!(sys_SUN, h*M)

    set_dipole!(sys_dip, M, position_to_site(sys_dip, (0, 0, 0)))
    set_coherent!(sys_SUN, Z, position_to_site(sys_SUN, (0, 0, 0)))

    energy(sys_dip)
    energy(sys_SUN)

    q = rand(3)

    swt_dip = SpinWaveTheory(sys_dip)
    swt_SUN = SpinWaveTheory(sys_SUN)

    disp_dip = dispersion(swt_dip, [q])
    disp_SUN = dispersion(swt_SUN, [q])

    @test disp_dip[1] ≈ disp_SUN[end-1]
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
    formula = intensity_formula(swt,:perp,kernel = delta_function_kernel)
    intensities_bands(swt,path,formula)
    @test_throws "without a finite-width kernel" intensities_broadened(swt,path,energies,formula)

    # Broadened
    formula = intensity_formula(swt,:perp,kernel = lorentzian(0.05))
    intensities_broadened(swt,path,energies,formula)
    @test_throws "broadening kernel" intensities_bands(swt,path,formula)

    formula = intensity_formula(swt,:perp,kernel = (w,dw) -> lorentzian(dw,w.^2))
    intensities_broadened(swt,path,energies,formula)
    @test_throws "broadening kernel" intensities_bands(swt,path,formula)

    # Full
    formula = intensity_formula(swt,:full,kernel = lorentzian(0.05))
    intensities_broadened(swt,path,energies,formula)
end

@testitem "Dipole-dipole unimplemented" begin
    sys = System(Sunny.diamond_crystal(),(1,1,1),[SpinInfo(1,S=1/2,g=2)],:SUN;seed = 0)
    enable_dipole_dipole!(sys)
    @test_throws "SpinWaveTheory does not yet support long-range dipole-dipole interactions." SpinWaveTheory(sys)
end

@testitem "langasite" begin
    a = b = 8.539
    c = 5.2414
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    crystal = Crystal(latvecs, [[0.24964,0,0.5]], 150)
    latsize = (1,1,7)
    rn_seed = 5
    S = 5/2
    J₁ = 0.85
    J₂ = 0.24
    J₃ = 0.053
    J₄ = 0.017
    J₅ = 0.24
    sys = System(crystal, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=rn_seed)
    set_exchange!(sys, J₁, Bond(3, 2, [1,1,0]))
    set_exchange!(sys, J₄, Bond(1, 1, [0,0,1]))
    set_exchange!(sys, J₂, Bond(1, 3, [0,0,0]))
    set_exchange!(sys, J₃, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₅, Bond(3, 2, [1,1,1]))
    R1=[0.5 0.5im 0;
    -0.5im 0.5 0;
    0 0 0];
    R2=[0 0 0;
    0 0 0;
    0 0 1];
    function R(site)
        return exp((site-1)*2π*im/7)*R1+exp(-(site-1)*2π*im/7)*conj(R1)+R2
    end
    S1=[1.0, 0, 0]*(5/2)
    S2=[-0.5, -sqrt(3)/2, 0]*(5/2)
    S3=[-0.5, sqrt(3)/2, 0]*(5/2)

    for site ∈ 1:7
        set_dipole!(sys,R(site)*S1,(1,1,site,1))
        set_dipole!(sys,R(site)*S2,(1,1,site,2))
        set_dipole!(sys,R(site)*S3,(1,1,site,3))
    end
    swt = SpinWaveTheory(sys);
    q =  [[0.41568,0.56382,0.76414]];
    formula = intensity_formula(swt, :perp; kernel=delta_function_kernel)
    disp, intensities = intensities_bands(swt, q, formula)
    SpinW_energies=[2.6267,2.6541,2.8177,2.8767,3.2458,3.3172,3.4727,3.7767,3.8202,3.8284,3.8749,3.9095,3.9422,3.9730,4.0113,4.0794,4.2785,4.4605,4.6736,4.7564,4.7865]
    SpinW_intensities = [0.0000,0.0000,0.0000,0.0197,0.0197,0.0548,0.0548,0.0933,0.0933,0.1187,0.1187,0.1380,0.1380,0.4847,0.4847,0.9621,0.9621,2.2151,2.2151,2.3716,2.3716]
    vec(disp)
    @test sum((vec(disp) .- reverse(SpinW_energies)).^2) < 1e-6
  end

  @testitem "langasite" begin
    tol = 1e-6
    a = b = 8.539
    c = 5.2414
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    crystal = Crystal(latvecs, [[0.24964,0,0.5]], 150)
    latsize = (1,1,7)
    rn_seed = 5
    S = 5/2
    J₁ = 0.85
    J₂ = 0.24
    J₃ = 0.053
    J₄ = 0.017
    J₅ = 0.24
    sys = System(crystal, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=rn_seed)
    set_exchange!(sys, J₁, Bond(3, 2, [1,1,0]))
    set_exchange!(sys, J₄, Bond(1, 1, [0,0,1]))
    set_exchange!(sys, J₂, Bond(1, 3, [0,0,0]))
    set_exchange!(sys, J₃, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₅, Bond(3, 2, [1,1,1]))
    R1=[0.5 0.5im 0;
    -0.5im 0.5 0;
    0 0 0];
    R2=[0 0 0;
    0 0 0;
    0 0 1];
    function R(site)
        return exp((site-1)*2π*im/7)*R1+exp(-(site-1)*2π*im/7)*conj(R1)+R2
    end
    S1=[1.0, 0, 0]*(5/2)
    S2=[-0.5, -sqrt(3)/2, 0]*(5/2)
    S3=[-0.5, sqrt(3)/2, 0]*(5/2)

    for site ∈ 1:7
        set_dipole!(sys,R(site)*S1,(1,1,site,1))
        set_dipole!(sys,R(site)*S2,(1,1,site,2))
        set_dipole!(sys,R(site)*S3,(1,1,site,3))
    end
    swt = SpinWaveTheory(sys);
    q =  [[0.41568,0.56382,0.76414]];
    formula = intensity_formula(swt, :full; kernel=delta_function_kernel)
    disp, intensities = intensities_bands(swt, q, formula)
    SpinW_energies=[2.6267,2.6541,2.8177,2.8767,3.2458,3.3172,3.4727,3.7767,3.8202,3.8284,3.8749,3.9095,3.9422,3.9730,4.0113,4.0794,4.2785,4.4605,4.6736,4.7564,4.7865]
    SpinW_intensities = [4.90888351357109e-28+0im,3.02833409310332e-29-4.90183847043242e-28im,-3.47607398564791e-29+9.84332735031517e-30im,3.02833409310332e-29+4.90183847043242e-28im,4.91348560162905e-28+0im,-1.1973621676136e-29-3.4103608905346e-29im,-3.47607398564791e-29-9.84332735031517e-30im,-1.1973621676136e-29+3.4103608905346e-29im,2.65885333210072e-30+0im,0.299983007942544+0im,-6.89924308159919e-16-0.299983007942545im,8.90247894275424e-16+6.34020574445034e-16im,-6.89924308159919e-16+0.299983007942545im,0.299983007942545+0im,-6.34020574445037e-16+8.90247894275424e-16im,8.90247894275424e-16-6.34020574445034e-16im,-6.34020574445037e-16-8.90247894275424e-16im,3.98197021316029e-30+0im,1.6246106384345e-29+0im,1.56739745027586e-30-1.39495896878909e-29im,1.0422183650694e-29-2.28452778255029e-30im,1.56739745027586e-30+1.39495896878909e-29im,1.21289238520266e-29+0im,2.96710658770636e-30+8.72851741843924e-30im,1.0422183650694e-29+2.28452778255029e-30im,2.96710658770636e-30-8.72851741843924e-30im,7.00727771595391e-30+0im,1.69824547461359e-30+0im,-7.42073704424708e-31-2.49624079042798e-30im,-7.22532316736406e-16-2.96399825511172e-16im,-7.42073704424708e-31+2.49624079042798e-30im,3.99346947657154e-30+0im,7.51397596306989e-16-9.3252957161112e-16im,-7.22532316736406e-16+2.96399825511172e-16im,7.51397596306989e-16+9.3252957161112e-16im,0.359138778467998+0im,2.54310795170859e-30+0im,-1.48431420813924e-30-2.81846090718612e-30im,1.87970762572194e-31-1.29075824872833e-30im,-1.48431420813924e-30+2.81846090718612e-30im,3.9899645420098e-30+0im,1.32080275586077e-30+9.61689043605461e-31im,1.87970762572194e-31+1.29075824872833e-30im,1.32080275586077e-30-9.61689043605461e-31im,6.69019914431516e-31+0im,1.21044229386368e-28+0im,-8.71795913355384e-29-1.42476032941553e-28im,1.04145194812635e-29+5.64988100146455e-30im,-8.71795913355384e-29+1.42476032941553e-28im,2.30491790064103e-28+0im,-1.4151076781575e-29+8.18928014087119e-30im,1.04145194812635e-29-5.64988100146455e-30im,-1.4151076781575e-29-8.18928014087119e-30im,1.15976921880538e-30+0im,0.595401813444212+0im,-2.53765262771464e-15-0.595401813444206im,-1.89275462543388e-15-1.1570502575895e-15im,-2.53765262771464e-15+0.595401813444206im,0.5954018134442+0im,1.1570502575895e-15-1.89275462543386e-15im,-1.89275462543388e-15+1.1570502575895e-15im,1.1570502575895e-15+1.89275462543386e-15im,8.26548602904178e-30+0im,1.85037289986544e-27+0im,-1.09744411303626e-27+7.62063238692483e-28im,1.34726290101707e-28+1.70282589058599e-28im,-1.09744411303626e-27-7.62063238692483e-28im,9.64737411110074e-28+0im,-9.77558233462628e-30-1.56479581984755e-28im,1.34726290101707e-28-1.70282589058599e-28im,-9.77558233462628e-30+1.56479581984755e-28im,2.54799091493924e-29+0im,1.37085060162835+0im,-1.90323947078598e-16+1.37085060162835im,1.15600795474678e-14+1.82058944590381e-16im,-1.90323947078598e-16-1.37085060162835im,1.37085060162836+0im,1.8205894459038e-16-1.15600795474679e-14im,1.15600795474678e-14-1.82058944590381e-16im,1.8205894459038e-16+1.15600795474679e-14im,9.75077695879568e-29+0im,1.82708956583492e-27+0im,3.9638753194611e-28+2.48130506786242e-27im,-9.64771586309842e-15+6.49411633741611e-16im,3.9638753194611e-28-2.48130506786242e-27im,3.45576814259609e-27+0im,-1.2111315960298e-15+1.32431082698793e-14im,-9.64771586309842e-15-6.49411633741611e-16im,-1.2111315960298e-15-1.32431082698793e-14im,0.051174369660632+0im,1.86554689156241e-28+0im,4.35843408231783e-29+1.08862288401565e-28im,4.48522989416621e-30+3.93726340360949e-29im,4.35843408231783e-29-1.08862288401565e-28im,7.37081049166241e-29+0im,2.40234156011529e-29+6.58122246190838e-30im,4.48522989416621e-30-3.93726340360949e-29im,2.40234156011529e-29-6.58122246190838e-30im,8.41748661074193e-30+0im,3.35502349365754e-28+0im,1.32802912014496e-29+4.59974670321186e-28im,-2.85886622958724e-29+3.06843423298937e-29im,1.32802912014496e-29-4.59974670321186e-28im,6.31152252351689e-28+0im,4.09366864727777e-29+4.04097245267084e-29im,-2.85886622958724e-29-3.06843423298937e-29im,4.09366864727777e-29-4.04097245267084e-29im,5.24240882190698e-30+0im,0.0734875342470339+0im,2.94209101525666e-15+0.0734875342470341im,-3.74449208872516e-16-2.75527516005471e-16im,2.94209101525666e-15-0.0734875342470341im,0.0734875342470342+0im,-2.75527516005487e-16+3.74449208872506e-16im,-3.74449208872516e-16+2.75527516005471e-16im,-2.75527516005487e-16-3.74449208872506e-16im,2.94101066685499e-30+0im,1.79675507168869e-30+0im,-2.2361697451431e-30+2.77057403478318e-30im,2.16990399841372e-30+1.82024453868899e-30im,-2.2361697451431e-30-2.77057403478318e-30im,7.05523853031007e-30+0im,1.06218475593571e-31-5.61136885100296e-30im,2.16990399841372e-30-1.82024453868899e-30im,1.06218475593571e-31+5.61136885100296e-30im,4.46458934184035e-30+0im,8.6816880879254e-30+0im,-2.2226284801167e-30+1.30416839575628e-30im,6.07940421909521e-31+6.89645420429688e-30im,-2.2226284801167e-30-1.30416839575628e-30im,7.64935632086538e-31+0im,8.80348596110482e-31-1.85690985990034e-30im,6.07940421909521e-31-6.89645420429688e-30im,8.80348596110482e-31+1.85690985990034e-30im,5.52089313312445e-30+0im,0.0577275935458864+0im,-1.1498738469332e-16-0.0577275935458863im,-7.01842617793228e-16+3.18606320181212e-16im,-1.1498738469332e-16+0.0577275935458863im,0.0577275935458862+0im,-3.1860632018121e-16-7.01842617793227e-16im,-7.01842617793228e-16-3.18606320181212e-16im,-3.1860632018121e-16+7.01842617793227e-16im,1.02913184305532e-29+0im,3.76351473089694e-30+0im,1.30564293586995e-30-3.6730097801532e-30im,-1.68432478419235e-29+2.9721226034633e-30im,1.30564293586995e-30+3.6730097801532e-30im,4.03763646687431e-30+0im,-8.74392829715495e-30-1.54071093956788e-29im,-1.68432478419235e-29-2.9721226034633e-30im,-8.74392829715495e-30+1.54071093956788e-29im,7.77274785808427e-29+0im,1.04357889700868e-30+0im,-1.65839764014816e-31+2.7832095108826e-31im,2.16517419221095e-15-1.32454657176033e-15im,-1.65839764014816e-31-2.7832095108826e-31im,1.005821212407e-31+0im,-6.97332075980267e-16-3.66959173530849e-16im,2.16517419221095e-15+1.32454657176033e-15im,-6.97332075980267e-16+3.66959173530849e-16im,6.1733740705613+0im,2.72585119014601e-30+0im,-1.99264966031779e-30+1.38215738123196e-30im,-1.03209897887589e-29+4.61760000501311e-29im,-1.99264966031779e-30-1.38215738123196e-30im,2.15749550691488e-30+0im,3.09586291455806e-29-2.85222314694335e-29im,-1.03209897887589e-29-4.61760000501311e-29im,3.09586291455806e-29+2.85222314694335e-29im,8.21301551215441e-28+0im,5.68531286772134e-29+0im,-3.0817217739414e-30+5.80647999884397e-29im,-1.03325316117828e-28+3.46865485641125e-30im,-3.0817217739414e-30-5.80647999884397e-29im,5.94693394269539e-29+0im,9.14332489769348e-30+1.05339398627098e-28im,-1.03325316117828e-28-3.46865485641125e-30im,9.14332489769348e-30-1.05339398627098e-28im,1.87995854688041e-28+0im,0.0338873034127333+0im,-2.57730345002268e-17+0.0338873034127331im,-1.10230561446743e-15+1.1574183812994e-15im,-2.57730345002268e-17-0.0338873034127331im,0.033887303412733+0im,1.1574183812994e-15+1.10230561446743e-15im,-1.10230561446743e-15-1.1574183812994e-15im,1.1574183812994e-15-1.10230561446743e-15im,7.53879689375408e-29+0im]
    SpinW_intensities_reshape = reshape(SpinW_intensities,(9,21))./3 # 9 x 21 to match static array
    sqdiffs= []
    for i in 1:length(intensities)
        push!(sqdiffs,(SpinW_intensities_reshape[:,i].-intensities[i]).^2)
    end
    @test sum((vec(disp) .- reverse(SpinW_energies)).^2) < tol && real(sum(sum(sqdiffs))) < tol && imag(sum(sum(sqdiffs))) < tol
  end
