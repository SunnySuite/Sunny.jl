@testitem "Lanczos consistency" begin
    using LinearAlgebra

    N = 40
    A = hermitianpart(randn(N, N))

    S = randn(N, N)
    S = S' * S
    S = S / eigmax(S)
    S = S + 1e-2I

    # Check that fast and slow Lanczos implementations match

    lhs = randn(N, 2)
    v = randn(N)
    v /= sqrt(dot(v, S, v))
    (; Q, T) = Sunny.lanczos_ref(A, S, v; niters=10)
    ev1 = eigvals(T)

    mulA!(w, v) = mul!(w, A, v)
    mulS!(w, v) = mul!(w, S, v)
    T, lhs_adj_Q = Sunny.lanczos(mulA!, mulS!, copy(v); lhs, min_iters=10)
    ev2 = eigvals(T)

    @test ev1 ≈ ev2
    @test lhs' * Q ≈ lhs_adj_Q

    # Check that extremal eigenvalues match to reasonable accuracy

    T, _ = Sunny.lanczos(mulA!, mulS!, copy(v); min_iters=20)
    @test isapprox(eigmin(T), eigmin(A * S); atol=1e-3)
    @test isapprox(eigmax(T), eigmax(A * S); atol=1e-3)
end


@testitem "Lanczos eigenbounds" begin
    using LinearAlgebra

    fcc = Sunny.fcc_crystal()
    s = 5/2
    g = 2
    J = 1.9
    K = 0.0129
    C = J + K
    J₁ = diagm([J, J, C])
    D = 25/24

    for mode in (:dipole, :SUN)
        sys = System(fcc, [1 => Moment(; s, g)], mode)
        set_exchange!(sys, J₁, Bond(1, 2, [0, 0, 0]))
        set_onsite_coupling!(sys, S -> D * (S[1]^4 + S[2]^4 + S[3]^4), 1)
        set_dipole!(sys, (1, 1, 1), position_to_site(sys, (0, 0, 0)))
        set_dipole!(sys, (1, -1, -1), position_to_site(sys, (1/2, 1/2, 0)))
        set_dipole!(sys, (-1, -1, 1), position_to_site(sys, (1/2, 0, 1/2)))
        set_dipole!(sys, (-1, 1, -1), position_to_site(sys, (0, 1/2, 1/2)))
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        q = [0.8, 0.6, 0.1]
        res = intensities_bands(swt, [q])
        energies, _ = excitations(swt, q)
        @test all(extrema(energies) .≈ Sunny.eigbounds(swt, q, 30))
    end

end


@testitem "KPM vs Lanczos" begin
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

        swt = SpinWaveTheory(sys; measure=nothing)
        energies, T = excitations(swt, q)
        q_reshaped = Sunny.to_reshaped_rlu(sys, q)
        bounds = Sunny.eigbounds(swt, q_reshaped, 50)
        @test all(extrema(energies) .≈ bounds)

        formfactors = [1 => FormFactor("Fe2")]
        measure = ssf_perp(sys; formfactors)
        kernel = lorentzian(fwhm=0.1)
        energies = range(0, 6, 100)
        kT = 0.0

        # Reference calculation
        swt = SpinWaveTheory(sys; measure)
        res1 = intensities(swt, [q]; energies, kernel, kT)

        # Note that KPM accuracy is currently limited by Gibbs ringing
        # introduced in the thermal occupancy (Heaviside step) function.
        swt_kry = SpinWaveTheoryKPM(sys; measure, tol=1e-8, method=:kpm)
        res2 = intensities(swt_kry, [q]; energies, kernel, kT)
        @test isapprox(res1.data, res2.data, rtol=1e-5)

        # Default Lanczos method does this task well
        swt_kry = SpinWaveTheoryKPM(sys; measure, tol=1e-2)
        res3 = intensities(swt_kry, [q]; energies, kernel, kT)
        @test isapprox(res1.data, res3.data, rtol=1e-8)

        # Very high accuracy after 4 iterations
        swt_kry = SpinWaveTheoryKPM(sys; measure, niters=4)
        res3 = intensities(swt_kry, [q]; energies, kernel, kT)
        @test isapprox(res1.data, res3.data, rtol=1e-8)
    end
end
