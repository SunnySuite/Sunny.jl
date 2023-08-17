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
