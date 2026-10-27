@testitem "Lanczos resolvent" begin
    using LinearAlgebra, Random

    # Lanczos with A = Ĩ in the inner product of S approximates χ = v† S (z -
    # ĨS)⁻¹ v / β₀² by the leading entry x₁ of the resolvent of T, and the probe
    # products vₐ† S (z - ĨS)⁻¹ v / β₀ by the sums sₐ. The error estimates are
    # approximate for x₁, and a Cauchy–Schwarz bound for sₐ, which shares the
    # approximation, and so holds within the factor 2 margin of `tol`. Once the
    # Krylov space is exhausted, the estimates vanish and the results are exact.
    Random.seed!(0)
    n = 12
    L = n ÷ 2
    S = let X = randn(ComplexF64, n, n); Hermitian(X'*X/n + I) end
    Ĩ = Diagonal([ones(L); -ones(L)])
    (v, probes) = (randn(ComplexF64, n), randn(ComplexF64, n, 2))
    mulS!(y, x) = (y .= S * x)
    zs = [0.3 + 0.1im, -1.2 + 0.01im]

    r = Sunny.extend!(Sunny.LanczosResolvent(mulS!, v, probes, zs), 4)
    r_exact = Sunny.extend!(Sunny.LanczosResolvent(mulS!, v, probes, zs), n + 1)
    @test (r.nrows, r_exact.nrows) == (4, n)
    @test Sunny.exhausted(r_exact)
    for (i, z) in enumerate(zs)
        R = S * inv(z*I - Ĩ * S) * v
        @test r.β₀^2 * r_exact.x₁[i] ≈ v' * R
        @test r.β₀ * r_exact.s[:, i] ≈ probes' * R
        @test Sunny.error_estimates(r_exact, i) == (0, 0)

        (δx, δe) = Sunny.error_estimates(r, i)
        @test 0.1 < δx / abs(r_exact.x₁[i] - r.x₁[i]) < 10
        for a in 1:2
            normS = sqrt(real(probes[:, a]' * S * probes[:, a]))
            @test abs(r_exact.s[a, i] - r.s[a, i]) <= 2 * normS * δe
        end
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
        energies = range(0, 6, 61)

        # Without a kernel, both calculators give each excitation the Lorentzian
        # line shape of half width η. Away from ω = 0, where the thermal factor
        # is smooth, this is a Lorentzian kernel of full width 2η.
        swt = SpinWaveTheory(sys; measure)
        swt_large = SpinWaveTheoryLarge(sys; measure, tol=1e-10)
        res1 = intensities(swt, [q]; energies, η=0.05, kT=0.3)
        res2 = intensities(swt_large, [q]; energies, η=0.05, kT=0.3)
        res3 = intensities(swt, [q]; energies, kernel=lorentzian(fwhm=0.1), kT=0.3)
        far = energies .> 1.0
        @test isapprox(res1.data, res2.data, rtol=1e-6)
        @test isapprox(res1.data[far], res3.data[far], rtol=1e-2)

        # With a kernel for instrument resolution, both calculators apply it
        # after the Lorentzian line shape and the thermal factor
        kernel = gaussian(fwhm=0.3)
        res1 = intensities(swt, [q]; energies, η=0.05, kernel, kT=0.3)
        res2 = intensities(swt_large, [q]; energies, η=0.05, kernel, kT=0.3)
        @test isapprox(res1.data, res2.data, rtol=1e-6)

        # The limit η → 0⁺, including ω = 0 at kT > 0, where the pole of the
        # thermal factor is evaluated close to z = 0
        energies0 = range(-5, 5, 41)
        res1 = intensities(swt, [q]; energies=energies0, kernel, kT=0.3)
        res2 = intensities(swt_large, [q]; energies=energies0, kernel, kT=0.3)
        @test maximum(abs, res1.data - res2.data) < 1e-6 * maximum(abs, res1.data)

        # A matrix-valued combiner needs all correlations, including the
        # off-diagonal ones
        measure = ssf_custom((q, ssf) -> ssf, sys; formfactors)
        swt = SpinWaveTheory(sys; measure)
        swt_large = SpinWaveTheoryLarge(sys; measure, tol=1e-10)
        res1 = intensities(swt, [q]; energies, η=0.05, kernel, kT=0.3)
        res2 = intensities(swt_large, [q]; energies, η=0.05, kernel, kT=0.3)
        @test isapprox(res1.data, res2.data, rtol=1e-6)
        @test_throws "Specify a broadening η > 0" intensities(swt_large, [q]; energies, kernel, kT=0.3)

        # A temperature far below the kernel width would enclose too many
        # Matsubara poles
        @test_throws "too small relative to the kernel width" intensities(swt_large, [q]; energies, η=0.05, kernel, kT=1e-9)
        @test_throws "too small relative to the kernel width" intensities(swt, [q]; energies, η=0.05, kernel, kT=1e-9)
    end
end


@testitem "SpinWaveTheoryLarge at Goldstone modes and negative energies" begin
    using LinearAlgebra

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
    energies = range(-2, 3, 51)
    qs = [[1/3, 1/3, 0], [0.2, 0.05, 0]]
    kT = 0.5

    # Without a kernel, both calculators give the poles of χ'' of either sign
    # their Lorentzian line shapes, times 1 + n_B(ω). The poles at negative energy
    # describe the opposite wavevector, so this tests detailed balance. The result
    # is insensitive to the gap that the regularization opens at the zero modes.
    swt = SpinWaveTheory(sys; measure)
    swt_large = SpinWaveTheoryLarge(sys; measure, tol=1e-8)
    res1 = intensities(swt, qs; energies, η=0.1, kT).data
    res2 = intensities(swt_large, qs; energies, η=0.1, kT).data
    res3 = intensities(SpinWaveTheory(sys; measure, regularization=1e-6), qs; energies, η=0.1, kT).data
    @test isapprox(res1, res2; rtol=1e-6)
    @test isapprox(res1, res3; rtol=1e-3)

    # Exceeding the iteration cap is an error
    swt_capped = SpinWaveTheoryLarge(sys; measure, tol=1e-6, maxiters=5)
    @test_throws "Exceeded maxiters" intensities(swt_capped, qs; energies, η=0.1)

    # The line shape of a single pole, explicitly
    (disp, residues) = Sunny.susceptibility_poles(swt, Sunny.QPoints(qs[2:2]); threaded=false)
    (ω, η) = (energies[40], 0.1)
    ref = sum(@. (η/π) / ((ω - disp)^2 + η^2) * residues) / (1 - exp(-ω/kT))
    @test res1[40, 2] ≈ ref

    # Likewise with a kernel for instrument resolution
    kernel = gaussian(fwhm=0.3)
    res1 = intensities(swt, qs; energies, η=0.01, kernel, kT).data
    res2 = intensities(swt_large, qs; energies, η=0.01, kernel, kT).data
    @test isapprox(res1, res2; rtol=1e-6)

    # A matrix-valued combiner obtains cross correlations from probes, whose
    # convergence is slower. Near the soft modes, at kT > 0, the thermal factor
    # amplifies their error, which a check of the diagonal alone would miss.
    sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(1, 1, 1), seed=0)
    J = diagm([1, 1, 0.95])
    set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
    sys = to_inhomogeneous(repeat_periodically(sys, (9, 9, 1)))
    for (s1, s2, offset) in symmetry_equivalent_bonds(sys, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys, (1.0 + randn(sys.rng)/3) * J, s1, s2; offset)
    end
    for site in eachsite(sys)
        θ = 2π/3 * mod(site[1] - site[2], 3)
        set_dipole!(sys, [cos(θ), sin(θ), 0], site)
    end
    minimize_energy!(sys; maxiters=10_000)
    measure = ssf_custom((q, ssf) -> ssf, sys)
    kernel = gaussian(fwhm=0.4)
    qs = [[0.1, 0.2, 0], [0.3, 0.3, 0], [0.2, 0.05, 0]]
    swt = SpinWaveTheory(sys; measure)
    swt_large = SpinWaveTheoryLarge(sys; measure, tol=1e-2)
    res1 = intensities(swt, qs; energies, η=0.05, kernel, kT).data
    res2 = intensities(swt_large, qs; energies, η=0.05, kernel, kT).data
    @test maximum(x -> maximum(abs, x), res1 - res2) < 1e-2 * maximum(x -> maximum(abs, x), res1)
end


@testitem "SpinWaveTheoryLarge detects instability" begin
    # An antiferromagnetic chain is unstable in the ferromagnetic state
    cryst = Crystal(lattice_vectors(1, 10, 10, 90, 90, 90), [[0, 0, 0]], "P1")
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2, 1, 1))
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    measure = ssf_trace(sys)
    energies = range(0, 5, 51)
    qs = [[0.3, 0, 0]]

    polarize_spins!(sys, [0, 0, 1])
    swt = SpinWaveTheoryLarge(sys; measure, tol=1e-2)
    @test_throws Sunny.InstabilityError intensities(swt, qs; energies, η=0.1)

    set_dipole!(sys, [0, 0, -1], (2, 1, 1, 1))
    swt = SpinWaveTheoryLarge(sys; measure, tol=1e-2)
    @test all(isfinite, intensities(swt, qs; energies, η=0.1).data)

    # SpinWaveTheoryKPM forwards with a warning
    swt_kpm = @test_logs (:warn, r"deprecated") match_mode=:any SpinWaveTheoryKPM(sys; measure, tol=1e-2)
    kernel = gaussian(fwhm=0.5)
    @test intensities(swt_kpm, qs; energies, kernel).data ≈ intensities(swt, qs; energies, kernel).data

    # With a kernel, the default η = 0 is the limit η → 0⁺
    qs = [[0.01, 0, 0]]
    @test_throws "Specify a `kernel`" intensities(swt, qs; energies)
    res0 = intensities(swt, qs; energies, kernel).data
    @test isapprox(res0, intensities(swt, qs; energies, η=1e-9, kernel).data; rtol=1e-4)
    @test all(isfinite, intensities(swt, qs; energies, η=0.05, kernel).data)

    # A kernel much narrower than η is an error
    @test_throws "Kernel width must exceed η/4" intensities(swt, qs; energies, η=2.5, kernel)
end
