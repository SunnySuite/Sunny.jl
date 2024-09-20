@testitem "Lanczos" begin
    using LinearAlgebra

    N = 200
    A = hermitianpart(randn(N, N))
    A = diagm(vcat(ones(N÷2), -ones(N÷2)))

    S = randn(N, N)
    S = S' * S
    S = S / eigmax(S)
    S = S + 0.0I

    v = randn(N)
    ev1 = eigvals(Sunny.lanczos_ref(A*S*A, S, v; niters=10))

    mulA!(w, v) = (w .= A * S * A * v)
    mulS!(w, v) = mul!(w, S, v)
    ev2 = eigvals(Sunny.lanczos(mulA!, mulS!, copy(v); niters=10))

    @test ev1 ≈ ev2

    ev3 = extrema(eigvals(Sunny.lanczos_ref(A*S*A, S, v; niters=100)))
    ev4 = extrema(eigvals((A*S)^2))
    @test isapprox(collect(ev3), collect(ev4); atol=1e-3)
end


@testitem "FCC KPM" begin
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


@testitem "AFM KPM" begin
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
        swt = SpinWaveTheory(sys; measure)
        swt_kpm = SpinWaveTheoryKPM(sys; measure, tol=1e-8)

        kernel = lorentzian(fwhm=0.1)
        energies = range(0, 6, 100)
        kT = 0.0
        res1 = intensities(swt, [q]; energies, kernel, kT)
        res2 = intensities(swt_kpm, [q]; energies, kernel, kT)

        # Note that KPM accuracy is currently limited by Gibbs ringing
        # introduced in the thermal occupancy (Heaviside step) function.
        @test isapprox(res1.data, res2.data, rtol=1e-5)
    end
end
