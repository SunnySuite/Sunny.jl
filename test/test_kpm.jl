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
        Sunny.eigbounds(swt, q, 20; extend=0.0)
        @test all(extrema(energies) .≈ Sunny.eigbounds(swt, q, 30; extend=0.0))
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
        extrema(energies)
        q_reshaped = Sunny.to_reshaped_rlu(sys, q)
        bounds = Sunny.eigbounds(swt, q_reshaped, 50; extend=0.0)
        @test all(extrema(energies) .≈ bounds)

        measure = ssf_perp(sys; formfactors=[1 => "Fe2"])
        σ = 0.05
        swt = SpinWaveTheory(sys; measure)
        swt_kpm = SpinWaveTheoryKPM(sys; measure, resolution=0.005, screening_factor=10)

        energies = range(0, 6, 100)
        kT = 0.2
        kernel = lorentzian(fwhm=2σ)
        formfactors = [FormFactor("Fe2")]

        res1 = intensities(swt, [q]; energies, formfactors, kernel, kT)
        res2 = intensities(swt_kpm, [q]; energies, formfactors, kernel, kT)

        @test isapprox(res1.data, res2.data, atol=1e-3)
    end
end
