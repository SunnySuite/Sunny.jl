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
    set_anisotropy!(sys, D * 𝒮[3]^2, 1)

    Δt  = abs(0.05 / D)
    λ = 0.1
    langevin = Langevin(Δt; kT=0, λ)

    randomize_spins!(sys)
    A = [1 1 1; -1 1 0; 0 0 1]
    sys_lswt = construct_magnetic_supercell(sys, [1 1 1; -1 1 0; 0 0 1])

    langevin.kT = 0
    for i in 1:50_000
        step!(sys_lswt, langevin)
    end

    sw_fields = SpinWaveFields(sys_lswt)

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
    ωk_num = lswt_dispersion_relation(sw_fields, k)

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
    sys = System(fcc, dims, infos, :dipole)

    J = 22.06 * Sunny.meV_per_K
    K = 0.15  * Sunny.meV_per_K
    C = J + K
    J₁ = diagm([J, J, C])
    D_ST = 0.2
    D = D_ST / cov_factor

    set_exchange!(sys, J₁, Bond(1, 2, [0, 0, 0]))
    Λ = D * (𝒮[1]^4 + 𝒮[2]^4 + 𝒮[3]^4)
    set_anisotropy!(sys, Λ, 1)

    Δt = abs(0.05 / D)
    λ  = 0.1
    langevin = Langevin(Δt; kT=0, λ)

    polarize_spin!(sys, (1, 1, 1), position_to_site(sys, (0, 0, 0)))
    polarize_spin!(sys, (1, -1, -1), position_to_site(sys, (1/2, 1/2, 0)))
    polarize_spin!(sys, (-1, -1, 1), position_to_site(sys, (1/2, 0, 1/2)))
    polarize_spin!(sys, (-1, 1, -1), position_to_site(sys, (0, 1/2, 1/2)))
    sw_fields = SpinWaveFields(sys)

    disp = zeros(Float64, 4)
    Sαβ_matrix = zeros(Float64, 4, 9)

    k = [0.8, 0.6, 0.1]
    lswt_dynamical_spin_structure_factor!(sw_fields, k, disp, Sαβ_matrix)
    tmp = Sαβ_matrix[:, 1:3]
    sunny_trace = sum(tmp, dims=2)

    spintools_trace = [0.0, 1.1743243223274487, 1.229979802236658, 1.048056653379038]

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
        sys = System(cryst, dims, infos, :dipole)

        α = -0.4 * π
        J = 1.0
        JL, JQ = J * cos(α), J * sin(α) / S^2
        set_exchange!(sys, JL,  Bond(1, 1, [1, 0, 0]))
        set_biquadratic!(sys, JQ,  Bond(1, 1, [1, 0, 0]))

        Δt  = abs(0.05 / JL)
        λ = 0.1
        langevin = Langevin(Δt; kT=0, λ)

        randomize_spins!(sys)

        langevin.kT = 0
        for i in 1:100_000
            step!(sys, langevin)
        end

        sys_lswt = construct_magnetic_supercell(sys, [1 1 1; -1 1 0; 0 0 1])
        langevin.kT = 0

        for i in 1:10_000
            step!(sys_lswt, langevin)
        end

        sw_fields = SpinWaveFields(sys_lswt)

        @inline γk(k :: Vector{Float64}) = 2 * (cos(2π*k[1]) + cos(2π*k[2]) + cos(2π*k[3]))
        @inline ϵk₁(k :: Vector{Float64}) = J * (S*cos(α) - (2*S-2+1/S) * sin(α)) * √(36 - γk(k)^2) 

        ϵk_num = lswt_dispersion_relation(sw_fields, k)
        ϵk_ana = ϵk₁(k)

        isapprox(ϵk_num[1], ϵk_ana)
    end

    k = rand(Float64, 3)
    @test test_biquad(k, 1)
    @test test_biquad(k, 3/2)

end
