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
    J, J′, D = 1.0, 0.1, 5.0

    a = b = 1.0
    c = 1.5
    lat_vecs = lattice_vectors(a, b, c, 90, 90, 90)
    types = ["A"]
    basis_vecs = [[0, 0, 0]]

    cryst = Crystal(lat_vecs, basis_vecs; types)

    # Spin System
    dims = (2, 2, 2)
    infos = [SpinInfo(1, S=1)]
    sys = System(cryst, dims, infos, :SUN)

    set_exchange!(sys, J,  Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 1, [0, 0, 1]))
    S = spin_operators(sys, 1)
    set_onsite_coupling!(sys, D * S[3]^2, 1)

    Δt  = abs(0.05 / D)
    λ = 0.1
    langevin = Langevin(Δt; kT=0, λ)

    randomize_spins!(sys)
    A = [1 1 1; -1 1 0; 0 0 1]
    sys_swt = reshape_geometry(sys, A)

    langevin.kT = 0
    for i in 1:50_000
        step!(sys_swt, langevin)
    end

    swt = SpinWaveTheory(sys_swt)

    function sion_analytical_disp(k :: Vector{Float64})
        # analytical solutions
        γkxy = cos(2*π*k[1]) + cos(2*π*k[2])
        γkz  = cos(2*π*k[3])
        x = 1/2 - D/(8*(2*J+J′))
        Ak₊ = -8 * (x-1) * x * (2*J+J′) - (x-1) * D + 2 * (2*x-1) * (J *γkxy + J′*γkz)
        Bk₊ = -2 * (J * γkxy + J′ * γkz)
        Ak₋ = -16 * (x-1) * x * (2*J+J′) - (2*x-1) * D - 2 * (1-2*x)^2*(J*γkxy + J′*γkz)
        Bk₋ = 2 * (1-2*x)^2 * (J*γkxy + J′*γkz)
        ωk₊ = √(Ak₊^2-Bk₊^2)
        ωk₋ = √(Ak₋^2-Bk₋^2)
        return ωk₊, ωk₋
    end

    k = rand(Float64, 3)
    ωk1, ωk2 = sion_analytical_disp(k)
    ωk3, ωk4 = sion_analytical_disp(k .+= 0.5)
    ωk_ana = [ωk1, ωk2, ωk3, ωk4]
    index  = sortperm(ωk_ana, rev=true)
    ωk_ana = ωk_ana[index]

    ωk_num = dispersion(swt, [k])'

    @test isapprox(ωk_ana, ωk_num)
end

@testitem "Intensities" begin
    using LinearAlgebra

    a = 8.289
    lat_vecs = lattice_vectors(a, a, a, 90, 90, 90)
    types = ["MzR1"]
    basis_vecs = [[0, 0, 0]]
    fcc = Crystal(lat_vecs, basis_vecs, 225; types)
    S = 5/2

    # According to a renormalized classical theory for spins (the details will be presented in a manuscript in preparation), the large-S expansion and the :dipole mode should produce the same results when apply the proper renormalization factor for the single-ion interaction strength.
    cov_factor = (1 - 3/S + 11/(4*S^2)- 3/(4*S^3))

    dims = (1, 1, 1)
    infos = [SpinInfo(1, S=S)]
    sys = System(fcc, dims, infos, :SUN)

    J = 22.06 * Sunny.meV_per_K
    K = 0.15  * Sunny.meV_per_K
    C = J + K
    J₁ = diagm([J, J, C])
    D_ST = 0.2
    D = D_ST / cov_factor

    set_exchange!(sys, J₁, Bond(1, 2, [0, 0, 0]))
    S = spin_operators(sys, 1)
    Λ = D * (S[1]^4 + S[2]^4 + S[3]^4)
    set_onsite_coupling!(sys, Λ, 1)

    polarize_spin!(sys, (1, 1, 1), position_to_site(sys, (0, 0, 0)))
    polarize_spin!(sys, (1, -1, -1), position_to_site(sys, (1/2, 1/2, 0)))
    polarize_spin!(sys, (-1, -1, 1), position_to_site(sys, (1/2, 0, 1/2)))
    polarize_spin!(sys, (-1, 1, -1), position_to_site(sys, (0, 1/2, 1/2)))

    swt = SpinWaveTheory(sys)

    k = [0.8, 0.6, 0.1]
    _, Sαβs =  Sunny.dssf(swt, [k])

    sunny_trace = [real(tr(Sαβs[1,a])) for a in axes(Sαβs)[2]]
    sunny_trace = filter(x -> abs(x) > 1e-12, sunny_trace)
    spintools_trace = [1.1743243223274487, 1.229979802236658, 1.048056653379038]

    @test isapprox(sunny_trace, spintools_trace)
end

@testitem "Biquadratic interactions" begin
    function test_biquad(k :: Vector{Float64}, S)

        a = 1.0
        lat_vecs = lattice_vectors(a, a, a, 90, 90, 90)
        types = ["A"]
        basis_vecs = [[0, 0, 0]]

        cryst = Crystal(lat_vecs, basis_vecs; types)

        # Spin System
        dims = (2, 2, 2)
        infos = [SpinInfo(1, S=S)]
        sys_SUN = System(cryst, dims, infos, :SUN)
        sys_dip = System(cryst, dims, infos, :dipole)

        α = -0.4 * π
        J = 1.0
        JL, JQ = J * cos(α), J * sin(α) / S^2
        set_exchange!(sys_SUN, JL,  Bond(1, 1, [1, 0, 0]))
        set_biquadratic!(sys_SUN, JQ,  Bond(1, 1, [1, 0, 0]))
        set_exchange!(sys_dip, JL,  Bond(1, 1, [1, 0, 0]))
        set_biquadratic!(sys_dip, JQ,  Bond(1, 1, [1, 0, 0]))

        sys_swt_SUN = reshape_geometry(sys_SUN, [1 1 1; -1 1 0; 0 0 1])
        polarize_spin!(sys_swt_SUN, ( 1, 0, 0), position_to_site(sys_swt_SUN, (0, 0, 0)))
        polarize_spin!(sys_swt_SUN, (-1, 0, 0), position_to_site(sys_swt_SUN, (0, 1, 0)))

        sys_swt_dip = reshape_geometry(sys_dip, [1 1 1; -1 1 0; 0 0 1])
        polarize_spin!(sys_swt_dip, ( 1, 0, 0), position_to_site(sys_swt_dip, (0, 0, 0)))
        polarize_spin!(sys_swt_dip, (-1, 0, 0), position_to_site(sys_swt_dip, (0, 1, 0)))

        swt_SUN = SpinWaveTheory(sys_swt_SUN)
        swt_dip = SpinWaveTheory(sys_swt_dip)

        γk(k :: Vector{Float64}) = 2 * (cos(2π*k[1]) + cos(2π*k[2]) + cos(2π*k[3]))
        ϵk₁(k :: Vector{Float64}) = J * (S*cos(α) - (2*S-2+1/S) * sin(α)) * √(36 - γk(k)^2) 

        ϵk_num_SUN = dispersion(swt_SUN, [k])
        ϵk_num_dip = dispersion(swt_dip, [k])
        ϵk_ana = ϵk₁(k)

        ϵk_num_SUN[end-1] ≈ ϵk_num_SUN[end] ≈ ϵk_num_dip[end] ≈ ϵk_ana
    end

    k = [0.12, 0.23, 0.34]
    @test test_biquad(k, 1)
    @test test_biquad(k, 3/2)
end

@testitem "Canted afm" begin

    function disp_analytical_canted_afm(J, D, h, S, qs)
        c₂ = 1 - 1/(2S)
        θ = acos(h / (2S*(4J+c₂*D)))
        disp = zeros(Float64, 2, length(qs))
        for (iq, q) in enumerate(qs)
            Jq = 2J*(cos(2π*q[1])+cos(2π*q[2]))
            ωq₊ = real(√Complex(4J*S*(4J*S+2D*S*c₂*sin(θ)^2) + cos(2θ)*(Jq*S)^2 + 2S*Jq*(4J*S*cos(θ)^2 + c₂*D*S*sin(θ)^2)))
            ωq₋ = real(√Complex(4J*S*(4J*S+2D*S*c₂*sin(θ)^2) + cos(2θ)*(Jq*S)^2 - 2S*Jq*(4J*S*cos(θ)^2 + c₂*D*S*sin(θ)^2)))
            disp[:, iq] = [ωq₊, ωq₋]
        end
        return Sunny.reshape_dispersions(disp)
    end

    function test_canted_afm(q :: Vector{Float64}, S)
        J, D, h = 1.0, 0.54, 0.76
        a = 1
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)

        dims = (2, 2, 1)
        sys_dip = System(cryst, dims, [SpinInfo(1; S=S, g=1)], :dipole; units=Units.theory)

        set_exchange!(sys_dip, J, Bond(1, 1, [1, 0, 0]))
        S = spin_operators(sys_dip,1)
        set_onsite_coupling!(sys_dip, D*S[3]^2, 1)
        set_external_field!(sys_dip, [0, 0, h])
        sys_swt_dip = reshape_geometry(sys_dip, [1 -1 0; 1 1 0; 0 0 1])
        c₂ = 1 - 1/(2S)
        θ = acos(h / (2S*(4J+D*c₂)))
        Sunny.polarize_spin!(sys_swt_dip, ( sin(θ), 0, cos(θ)), position_to_site(sys_swt_dip, (0,0,0)))
        Sunny.polarize_spin!(sys_swt_dip, (-sin(θ), 0, cos(θ)), position_to_site(sys_swt_dip, (1,0,0)))
        swt_dip = SpinWaveTheory(sys_swt_dip)
        ϵq_num = dispersion(swt_dip, [q])
        ϵq_ana = disp_analytical_canted_afm(J, D, h, S, [q])
        ϵq_num ≈ ϵq_ana
    end
    q = [0.12, 0.23, 0.34]
    @test test_canted_afm(q, 1)
    @test test_canted_afm(q, 2)
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
    𝐌 = normalize(rand(3))
    θ, ϕ = Sunny.dipole_to_angles(𝐌)
    s_mat = Sunny.spin_matrices(N=2S+1)
    
    
    S = spin_operators(sys_dip,1)
    s̃ᶻ = 𝐌[1] * S[1] + 𝐌[2] * S[2] + 𝐌[3] * S[3]
    
    U_mat = exp(-1im * ϕ * s_mat[3]) * exp(-1im * θ * s_mat[2])
    hws = zeros(2S+1)
    hws[1] = 1.0
    Z = U_mat * hws

    aniso = Ds[1]*s̃ᶻ^2 + Ds[2]*s̃ᶻ^4 + Ds[3]*s̃ᶻ^6
    set_onsite_coupling!(sys_dip, aniso, 1)
    
    
    S = spin_operators(sys_SUN,1)
    s̃ᶻ = 𝐌[1] * S[1] + 𝐌[2] * S[2] + 𝐌[3] * S[3]
    aniso = Ds[1]*s̃ᶻ^2 + Ds[2]*s̃ᶻ^4 + Ds[3]*s̃ᶻ^6
    set_onsite_coupling!(sys_SUN, aniso, 1)
    set_external_field!(sys_dip, h*𝐌)
    set_external_field!(sys_SUN, h*𝐌)

    Sunny.polarize_spin!(sys_dip, 𝐌, position_to_site(sys_dip, (0, 0, 0)))
    Sunny.set_coherent_state!(sys_SUN, Z, position_to_site(sys_SUN, (0, 0, 0)))

    energy(sys_dip)
    energy(sys_SUN)

    q = rand(3)

    swt_dip = SpinWaveTheory(sys_dip)
    swt_SUN = SpinWaveTheory(sys_SUN)

    disp_dip = dispersion(swt_dip, [q])
    disp_SUN = dispersion(swt_SUN, [q])

    @test disp_dip[1] ≈ disp_SUN[end-1]
end