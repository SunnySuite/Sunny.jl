@testitem "FCC KPM" begin
    using LinearAlgebra

    fcc = Sunny.fcc_crystal()
    S = 5/2
    g = 2
    J = 1.9
    K = 0.0129
    C = J + K
    J₁ = diagm([J, J, C])
    D = 25/24

    for mode in (:dipole, :SUN)
        sys = System(fcc, (1, 1, 1), [SpinInfo(1; S, g)], mode)
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
    a = 1
    latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    q = Sunny.Vec3(0.12, 0.23, 0.34)
    
    S = 1
    sys = System(cryst, (1, 1, 1), [SpinInfo(1; S, g=-1)], :dipole)
    set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
    set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
    set_field!(sys, [0, 0, h])
    
    sys_swt_dip = reshape_supercell(sys, [1 -1 0; 1 1 0; 0 0 1])
    
    c₂ = 1 - 1/(2S)
    θ = acos(h / (2S*(4J+D*c₂)))
    set_dipole!(sys_swt_dip, ( sin(θ), 0, cos(θ)), position_to_site(sys_swt_dip, (0,0,0)))
    set_dipole!(sys_swt_dip, (-sin(θ), 0, cos(θ)), position_to_site(sys_swt_dip, (1,0,0)))
    
    swt_dip = SpinWaveTheory(sys_swt_dip; measure=ssf_perp(sys_swt_dip))
    ϵq_num = dispersion(swt_dip, [q])
    
    energies, T = excitations(swt_dip, q)
    extrema(energies)
    
    q_reshaped = Sunny.to_reshaped_rlu(sys_swt_dip, q)
    bounds = Sunny.eigbounds(swt_dip, q_reshaped, 50; extend=0.0)
    
    @test all(extrema(energies) .≈ bounds)
    
    energies = collect(range(0, 6, 100))
    P = 1000
    kT = 0.01
    σ = 0.05
    broadening = (ω, x, σ) ->  (1/π) * (σ / ((x - ω)^2 + σ^2))
    kernel = nothing # "jackson"
    regularization_style = :cubic
    
    res1 = Sunny.kpm_intensities(swt_dip, [q], energies, P, kT, σ, broadening; kernel, regularization_style)
    res2 = intensities(swt_dip, [q]; energies, kernel=lorentzian(fwhm=2σ))
    
    @test isapprox(res1[1,:], res2.data[:,1], atol=1e-3)
end
