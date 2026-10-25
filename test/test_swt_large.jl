@testitem "Block Lanczos moments" begin
    using LinearAlgebra, Random

    # Lanczos with A = Ĩ in the inner product of S builds a block tridiagonal T
    # whose moments R₀†[Tᵏ]₁₁R₀ reproduce V†S(ĨS)ᵏV for k < 2m, where m is the
    # number of steps. A repeated starting row is linearly dependent, and must
    # be handled by a rank-revealing first step.
    Random.seed!(0)
    n = 40
    L = n ÷ 2
    S = let X = randn(ComplexF64, n, n); Hermitian(X'*X/n + I) end
    Ĩ = Diagonal([ones(L); -ones(L)])
    V = randn(ComplexF64, 3, n)
    V = [V; V[1:1, :]]
    mulS!(Y, X) = (Y .= transpose(S * transpose(X)))
    lz = Sunny.BlockLanczos(mulS!, V)
    Sunny.extend!(lz, 6)
    T = Sunny.lanczos_matrix(lz)
    r = size(lz.R₀, 1)
    @test r == 3
    Vc = transpose(V)
    for k in 0:2Sunny.nsteps(lz)-1
        @test lz.R₀' * (Matrix(T)^k)[1:r, 1:r] * lz.R₀ ≈ Vc' * S * (Ĩ * S)^k * Vc
    end
end


@testitem "SpinWaveTheoryLarge vs dense SWT" begin
    J, D, h = 1.0, 0.54, 0.76
    s, g = 1, -1.5
    a = 1
    latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    q = Sunny.Vec3(0.12, 0.23, 0.34)

    for mode in (:dipole, :SUN)
        sys = System(cryst, [1 => Moment(; s, g)], mode)
        sys = reshape_supercell(sys, [1 -1 0; 1 1 0; 0 0 1])
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
        set_field!(sys, [0, 0, -h/g])

        c₂ = 1 - 1/(2s)
        θ = acos(h / (2s*(4J+D*c₂)))
        set_dipole!(sys, ( sin(θ), 0, cos(θ)), position_to_site(sys, (0,0,0)))
        set_dipole!(sys, (-sin(θ), 0, cos(θ)), position_to_site(sys, (1,0,0)))

        formfactors = [1 => FormFactor("Fe2")]
        measure = ssf_perp(sys; formfactors)
        energies = range(0, 6, 100)
        η = 0.05

        # The two calculators agree at the same η, with or without a kernel,
        # including one whose width varies with energy
        swt = SpinWaveTheory(sys; measure)
        swt_large = SpinWaveTheoryLarge(sys; measure, tol=1e-10)
        σ(ϵ) = (0.1 + 0.05ϵ) / 2√(2log(2))
        kernels = [nothing, lorentzian(fwhm=0.3), Sunny.NonstationaryBroadening((ϵ, ω) -> exp(-(ω-ϵ)^2/2σ(ϵ)^2) / √(2π*σ(ϵ)^2))]
        for kernel in kernels, kT in (0.0, 0.3)
            res1 = intensities(swt, [q]; energies, η, kernel, kT)
            res2 = intensities(swt_large, [q]; energies, η, kernel, kT)
            @test isapprox(res1.data, res2.data, rtol=1e-4)
        end

        # Away from ω = 0, where the thermal factor is smooth, a Lorentzian kernel
        # of width Γ combines with η into a Lorentzian of width Γ + 2η
        far = energies .> 1.0
        res1 = intensities(swt, [q]; energies, kernel=lorentzian(fwhm=0.3+2η), kT=0.3)
        res2 = intensities(swt, [q]; energies, η, kernel=lorentzian(fwhm=0.3), kT=0.3)
        @test isapprox(res1.data[far], res2.data[far], rtol=1e-2)
    end
end


@testitem "SpinWaveTheoryLarge at Goldstone modes and negative energies" begin
    # A disordered antiferromagnet with continuous symmetry has exact zero modes
    # at the ordering wavevector of its supercell
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(1, 1, 1), seed=0)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    sys = to_inhomogeneous(repeat_periodically(sys, (6, 6, 1)))
    for (s1, s2, offset) in symmetry_equivalent_bonds(sys, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys, 1.0 + randn(sys.rng)/3, s1, s2; offset)
    end
    for site in eachsite(sys)
        θ = 2π/3 * mod(site[1] - site[2], 3)
        set_dipole!(sys, [cos(θ), sin(θ), 0], site)
    end
    minimize_energy!(sys; maxiters=10_000)

    measure = ssf_perp(sys)
    energies = range(-2, 3, 201)
    η = 0.1
    qs = [[1/3, 1/3, 0], [0.2, 0.05, 0]]

    # The broadened susceptibility is insensitive to the gap that the
    # regularization opens at the exact zero modes
    res = map((1e-8, 1e-6)) do regularization
        swt = SpinWaveTheoryLarge(sys; measure, tol=1e-8, regularization)
        intensities(swt, qs; energies, η, kT=0.5).data
    end
    @test isapprox(res[1], res[2]; rtol=1e-3)

    # Exact reference: the poles of χ'' of either sign, as Lorentzians, times 1 +
    # n_B(ω). The poles at negative energy describe the opposite wavevector, so
    # this tests detailed balance at kT > 0.
    swt = SpinWaveTheory(sys; measure)
    (disp, residues) = Sunny.susceptibility_poles(swt, Sunny.QPoints(qs); threaded=false)
    K(x) = (η/π) / (x^2 + η^2)
    ref = [sum(K.(ω .- disp[:, iq]) .* residues[:, iq]) / (1 - exp(-ω/0.5)) for ω in energies, iq in eachindex(qs)]
    nonzero = .!iszero.(energies)
    @test isapprox(res[1][nonzero, :], ref[nonzero, :]; rtol=1e-8)

    # Dense SWT interpolates a tabulated line shape, here also with a kernel
    swt_large = SpinWaveTheoryLarge(sys; measure, tol=1e-8)
    for kernel in (nothing, gaussian(fwhm=0.3)), kT in (0.0, 0.5)
        res1 = intensities(swt, qs; energies, η, kernel, kT).data
        res2 = intensities(swt_large, qs; energies, η, kernel, kT).data
        @test isapprox(res1, res2; rtol=1e-3)
    end
end


@testitem "SpinWaveTheoryLarge detects instability" begin
    # An antiferromagnetic chain is unstable in the ferromagnetic state
    cryst = Crystal(lattice_vectors(1, 10, 10, 90, 90, 90), [[0, 0, 0]], "P1")
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2, 1, 1))
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    measure = ssf_trace(sys)
    energies = range(0, 5, 11)

    polarize_spins!(sys, [0, 0, 1])
    swt = SpinWaveTheoryLarge(sys; measure, tol=1e-2)
    @test_throws Sunny.InstabilityError intensities(swt, [[0.3, 0, 0]]; energies, η=0.1)

    set_dipole!(sys, [0, 0, -1], (2, 1, 1, 1))
    swt = SpinWaveTheoryLarge(sys; measure, tol=1e-2)
    @test all(isfinite, intensities(swt, [[0.3, 0, 0]]; energies, η=0.1).data)

    # The deprecated interface infers η from the kernel width
    swt_kpm = @test_logs (:warn, r"deprecated") SpinWaveTheoryKPM(sys; measure, tol=1e-2)
    kernel = gaussian(fwhm=0.5)
    @test intensities(swt_kpm, [[0.3, 0, 0]]; energies, kernel).data ≈ intensities(swt, [[0.3, 0, 0]]; energies, η=0.05, kernel).data

    # A soft mode, unresolved by the kernel, is regularized by η
    q = [0.01, 0, 0]
    @test_throws "η` is required" intensities(swt, [q]; energies, kernel)
    @test all(isfinite, intensities(swt, [q]; energies, kernel, η=0.05).data)
end
