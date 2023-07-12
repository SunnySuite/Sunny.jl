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
    set_onsite!(sys, D * S[3]^2, 1)

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
    set_onsite!(sys, Λ, 1)

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
        sys = System(cryst, dims, infos, :SUN)

        α = -0.4 * π
        J = 1.0
        JL, JQ = J * cos(α), J * sin(α) / S^2
        set_exchange!(sys, JL,  Bond(1, 1, [1, 0, 0]))
        set_biquadratic!(sys, JQ,  Bond(1, 1, [1, 0, 0]))

        sys_swt = reshape_geometry(sys, [1 1 1; -1 1 0; 0 0 1])
        polarize_spin!(sys_swt, (1, 0, 0), position_to_site(sys_swt, (0, 0, 0)))
        polarize_spin!(sys_swt, (-1, 0, 0), position_to_site(sys_swt, (0, 1, 0)))

        swt = SpinWaveTheory(sys_swt)

        γk(k :: Vector{Float64}) = 2 * (cos(2π*k[1]) + cos(2π*k[2]) + cos(2π*k[3]))
        ϵk₁(k :: Vector{Float64}) = J * (S*cos(α) - (2*S-2+1/S) * sin(α)) * √(36 - γk(k)^2) 

        ϵk_num = dispersion(swt, [k])
        ϵk_ana = ϵk₁(k)

        ϵk_num[end-1] ≈ ϵk_num[end] ≈ ϵk_ana
    end

    k = [0.12, 0.23, 0.34]
    @test test_biquad(k, 1)
    @test test_biquad(k, 3/2)
end
