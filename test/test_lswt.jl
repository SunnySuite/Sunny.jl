@testitem "Kitchen sink" begin
    using LinearAlgebra

    # Pyrochlore with nonstandard, primitive lattice vectors
    latvecs = [[1, 1, 0] [1, 0, 1] [0, 1, 1]] / 2
    positions = [[5, 5, 1], [5, 1, 5], [1, 5, 5], [5, 5, 5]] / 8

    msg = "Cell is 1/4 the standard size for spacegroup 227. Consider `standardize`."
    cryst = @test_logs (:info, msg) Crystal(latvecs, positions)
    natoms = Sunny.natoms(cryst)

    moments = [1 => Moment(s=5/2, g=7.2)]
    sys = System(cryst, moments, :SUN; seed=0)

    A,B,C,D = 2.6, -1.3, 0.2, -5.7
    set_exchange!(sys, [A C -D; C A D; D -D B], Bond(1, 2, [0, 0, 0]))

    A,B,C,D,E,F,G,H,K = 2.6, -1.3, 0.2, -5.7, 8.2, 0.3, 2.5, -0.6, 1.3
    set_exchange!(sys, [A F+K E-H; F-K B D+G; E+H D-G C], Bond(1, 4, [1, 0, 0]))

    A,B,C,D = 2.6, -1.3, 0.2, -5.7
    set_exchange!(sys, [A D D; D B C; D C B], Bond(4, 4, [1, 1, 0]))

    O = stevens_matrices(spin_label(sys, 3))
    c1, c2, c3 = 2.6, -1.3, 0.2, -5.7
    Λ = c1 * (O[2,-2] - 2O[2,-1] - 2O[2,1]) +
        c2 * (-7O[4,-3] + 2O[4,-2] + O[4,-1] + O[4,1] + 7O[4,3]) +
        c3 * (O[4,0] + 5O[4,4])
    set_onsite_coupling!(sys, Λ, 3)

    A = [1 3 1; -1 1 0; 0 0 1]
    sys = reshape_supercell(sys, A)

    # minimize_energy!(sys; maxiters=1000)
    # println(round(reinterpret(reshape, ComplexF64, sys.coherents); digits=3))
    ground_state = ComplexF64[0.452 + 0.069im; 0.586 + 0.166im; 0.31 + 0.351im; 0.035 + 0.39im; -0.136 + 0.145im; -0.03 - 0.076im;;;;; -0.026 - 0.064im; 0.247 - 0.209im; -0.444 + 0.307im; 0.501 - 0.122im; -0.478 - 0.012im; 0.32 - 0.046im;;;;; -0.078 - 0.451im; -0.029 - 0.609im; 0.234 - 0.406im; 0.358 - 0.158im; 0.181 + 0.083im; -0.063 + 0.053im;;;;; -0.061 + 0.033im; -0.235 - 0.222im; 0.354 + 0.407im; -0.177 - 0.485im; 0.042 + 0.476im; -0.082 - 0.312im;;;;; -0.782 + 0.332im; 0.184 - 0.453im; 0.058 + 0.138im; 0.121 + 0.049im; -0.007 + 0.003im; -0.006 + 0.015im;;;;; 0.057 - 0.005im; -0.023 - 0.013im; -0.101 + 0.009im; -0.247 - 0.257im; -0.631 - 0.193im; -0.617 + 0.209im;;;;; -0.695 - 0.489im; 0.482 - 0.083im; -0.087 + 0.122im; 0.022 + 0.129im; -0.006 - 0.004im; -0.016 + 0.003im;;;;; 0.04 + 0.041im; -0.004 - 0.026im; -0.07 - 0.072im; 0.044 - 0.354im; -0.248 - 0.611im; -0.551 - 0.347im;;;;; 0.347 + 0.149im; -0.216 + 0.511im; -0.386 - 0.36im; 0.434 - 0.153im; -0.024 + 0.24im; 0.026 - 0.021im;;;;; 0.061 + 0.082im; 0.096 + 0.219im; 0.361 + 0.203im; 0.487 - 0.18im; 0.228 - 0.554im; -0.282 - 0.229im;;;;; -0.205 + 0.317im; -0.468 - 0.298im; 0.419 - 0.32im; 0.078 + 0.454im; -0.233 - 0.064im; 0.016 + 0.029im;;;;; 0.101 + 0.018im; 0.22 + 0.095im; 0.402 - 0.097im; 0.233 - 0.464im; -0.211 - 0.561im; -0.363 + 0.025im;;;;; 0.059 - 0.582im; 0.442 - 0.355im; 0.51 + 0.063im; 0.158 + 0.213im; -0.009 + 0.024im; 0.03 - 0.026im;;;;; 0.582 - 0.058im; 0.436 + 0.363im; 0.04 + 0.512im; -0.177 + 0.197im; -0.025 - 0.004im; 0.031 + 0.024im;;;;; 0.561 + 0.165im; 0.268 + 0.5im; -0.156 + 0.49im; -0.238 + 0.116im; -0.022 - 0.013im; 0.02 + 0.034im;;;;; -0.337 + 0.478im; -0.56 + 0.092im; -0.413 - 0.306im; -0.033 - 0.263im; 0.019 - 0.017im; -0.039 + 0.008im]

    sys.coherents .= reinterpret(reshape, Sunny.CVec{6}, ground_state)
    minimize_energy!(sys)
    @test energy_per_site(sys) ≈ -328.38255

    # Verify that this is a local minimum of energy
    @test norm(Sunny.proj.(Sunny.energy_grad_coherents(sys), sys.coherents)) < 1e-7

    # Test energies at an arbitrary wave vector
    qs = [[0.24331089495721447, 0.2818361515716459, 0.21954858411037714]]
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys; apply_g=false))

    res = intensities_bands(swt, qs; kT=100.0)
    # println(round.(vec(res.disp); digits=12))
    # println(round.(vec(res.data); digits=12))
    disps_golden = [1394.440092579925, 1393.728009951747, 1393.008551251224, 1392.919524974098, 1279.239919068637, 1279.094568472202, 1278.224518515366, 1277.69176148292, 1194.366336255056, 1193.750083625346, 1191.583519659758, 1189.794451340792, 1131.422439587906, 1131.202770074575, 1065.242927850642, 1065.09589244662, 1026.649340922954, 1024.028348558092, 1022.830406299696, 1020.767349639404, 945.202397530689, 944.795817851584, 835.545028394171, 832.001588695236, 827.939501409146, 827.307586947863, 821.216582166429, 820.430993567084, 820.29454877608, 818.594570998006, 810.207001090189, 808.553158273681, 766.524411070996, 766.51610275022, 766.513825852464, 766.508655558188, 758.579854167682, 754.683765885699, 750.572578891218, 750.471006252543, 665.954573008179, 662.421047653194, 651.465562550036, 651.417940124412, 581.258189152574, 568.105209800088, 559.053702296449, 558.493005822973, 552.043762736839, 550.131096070953, 539.733572947827, 530.698033192904, 499.661483510112, 494.928560823174, 435.233706061902, 427.70227706432, 408.128705853668, 399.856401749667, 370.0693430633, 369.845327686247, 365.04951424025, 363.639416669404, 354.648012591371, 346.609483926993, 341.989165167298, 339.373361067981, 318.363717384318, 276.219249203163, 263.161053831818, 257.409506246762, 230.539454193868, 229.778324172696, 203.971681278995, 197.504237153905, 193.879371534689, 189.866421874996, 189.815806967662, 167.944134431612, 154.923566498395, 146.219538847758] 
    data_golden = [0.0003866506, 0.0, 0.007231543496, 0.0, 0.008665025993, 0.0, 0.015340573883, 0.0, 0.0, 0.054200622367, 0.073310127326, 0.0, 0.00527607845, 0.0, 0.0, 0.026709096225, 0.0, 0.062334999507, 0.112031767918, 0.0, 0.031132103909, 0.0, 0.115596282255, 0.0, 0.0, 0.039004932814, 0.0, 0.031161016058, 0.029121117571, 0.0, 0.0, 0.004672396604, 0.000763285324, 0.000940997393, 0.0, 0.0, 0.008777727636, 0.0, 0.033281383926, 0.0, 0.066167396885, 0.0, 0.0, 0.060045694629, 0.335530855821, 0.007634795842, 0.0, 0.0, 0.0, 0.068726684036, 0.02950016123, 0.0, 0.0, 0.610003683362, 0.383607563841, 0.0, 0.873304893237, 0.0, 0.0, 0.297421427631, 0.358428673038, 0.0, 0.0, 0.986474626382, 0.0, 1.661999771965, 0.0, 0.215308474961, 0.195882715465, 0.0, 0.0, 0.348754438732, 0.884436719117, 0.0, 0.012121723144, 0.0, 0.177850646384, 0.402799537928, 0.0, 0.0]
    @test isapprox(res.disp, disps_golden; atol=1e-8)
    @test isapprox(res.data, data_golden; atol=1e-8)

    # Test first 5 output matrices
    formfactors = [1 => FormFactor("Fe2")]
    measure = ssf_custom((q, ssf) -> ssf, sys; apply_g=false, formfactors)
    swt = SpinWaveTheory(sys; measure)
    res = intensities_bands(swt, qs)
    data_flat = reinterpret(ComplexF64, res.data[1:5])
    # println(round.(data_flat; digits=12))
    data_golden = [0.003075023211 + 0.0im, 0.001813252796 - 0.000195741551im, 0.001874141877 + 0.000343251171im, 0.001813252796 + 0.000195741551im, 0.001081683041 + 0.0im, 0.001083277833 + 0.000321704428im, 0.001874141877 - 0.000343251171im, 0.001083277833 - 0.000321704428im, 0.001180553411 + 0.0im, 0.0 + 0.0im, 0.0 - 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 - 0.0im, 0.0 - 0.0im, 0.0 + 0.0im, 0.00087176192 + 0.0im, -0.000457598014 - 0.000845175127im, -0.000504075741 + 0.000799582738im, -0.000457598014 + 0.000845175127im, 0.001059597714 + 0.0im, -0.000510602006 - 0.000908412874im, -0.000504075741 - 0.000799582738im, -0.000510602006 + 0.000908412874im, 0.001024849661 + 0.0im, 0.0 + 0.0im, -0.0 - 0.0im, -0.0 + 0.0im, -0.0 + 0.0im, 0.0 + 0.0im, -0.0 - 0.0im, -0.0 - 0.0im, -0.0 + 0.0im, 0.0 + 0.0im, 0.000300068597 + 0.0im, 0.000826500187 + 0.000634405423im, 0.000960221539 - 0.000440668562im, 0.000826500187 - 0.000634405423im, 0.003617748774 + 0.0im, 0.001713144131 - 0.003243866269im, 0.000960221539 + 0.000440668562im, 0.001713144131 + 0.003243866269im, 0.003719863381 + 0.0im] 

    @test isapprox(data_flat, data_golden; atol=1e-9)
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
    moments = [1 => Moment(s=1, g=2)]
    sys = System(cryst, moments, :SUN; seed=0)
    set_exchange!(sys, J,  Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 1, [0, 0, 1]))
    set_onsite_coupling!(sys, S -> D * S[3]^2, 1)

    # Reshape to sheared supercell and minimize energy
    A = [1 1 1; -1 1 0; 0 0 1]
    sys = reshape_supercell(sys, A)
    randomize_spins!(sys)
    @test minimize_energy!(sys).converged

    q = rand(Float64, 3)
    swt = SpinWaveTheory(sys; measure=nothing)
    ωk_num = dispersion(swt, [q])

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
    ωk1, ωk2 = single_ion_analytical_disp(q)
    ωk3, ωk4 = single_ion_analytical_disp(q + [0.5, 0.5, 0.5])
    ωk_ref = sort([ωk1, ωk2, ωk3, ωk4]; rev=true)

    @test ωk_num ≈ ωk_ref
end


@testitem "Intensities" begin
    using LinearAlgebra

    # Crystal
    a = 8.289
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [[0, 0, 0]]
    fcc = Crystal(latvecs, positions, 225)

    units = Units(:meV, :angstrom)
    moments = [1 => Moment(s=5/2, g=2)]
    J = 22.06 * units.K
    K = 0.15  * units.K
    C = J + K
    J₁ = diagm([J, J, C])
    D = 25/24

    function compute(mode)
        sys = System(fcc, moments, mode)
        set_exchange!(sys, J₁, Bond(1, 2, [0, 0, 0]))
        set_onsite_coupling!(sys, S -> D * (S[1]^4 + S[2]^4 + S[3]^4), 1)
        set_dipole!(sys, (1, 1, 1), position_to_site(sys, (0, 0, 0)))
        set_dipole!(sys, (1, -1, -1), position_to_site(sys, (1/2, 1/2, 0)))
        set_dipole!(sys, (-1, -1, 1), position_to_site(sys, (1/2, 0, 1/2)))
        set_dipole!(sys, (-1, 1, -1), position_to_site(sys, (0, 1/2, 1/2)))
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        q = [0.8, 0.6, 0.1]
        res = intensities_bands(swt, [q])

        return filter(>(1e-12), abs.(res.data))
    end

    reference = [18.78918915723918, 19.679676835786527, 16.76890645406461]
    @test compute(:SUN) ≈ reference
    @test compute(:dipole) ≈ reference
end


@testitem "Scalar biquadratic" begin
    # Cubic crystal
    a = 2.0
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)

    function test_biquad(mode, q, s)
        # System
        sys = System(cryst, [1 => Moment(; s, g=2)], mode; dims=(2, 2, 2))
        α = -0.4π
        J = 1.0
        JL, JQ = J * cos(α), J * sin(α) / s^2
        set_pair_coupling!(sys, (Si, Sj) -> Si'*JL*Sj + JQ*(Si'*Sj)^2, Bond(1, 1, [1, 0, 0]))

        # Initialize Néel order
        sys = reshape_supercell(sys, [1 1 1; -1 1 0; 0 0 1])
        set_dipole!(sys, ( 1, 0, 0), position_to_site(sys, (0, 0, 0)))
        set_dipole!(sys, (-1, 0, 0), position_to_site(sys, (0, 1, 0)))

        # Numerical result
        swt = SpinWaveTheory(sys; measure=nothing)
        disp = dispersion(swt, [q])

        # Analytical result
        γq = 2 * (cos(2π*q[1]) + cos(2π*q[2]) + cos(2π*q[3]))
        disp_ref = J * (s*cos(α) - (2*s-2+1/s) * sin(α)) * √(36 - γq^2)

        @test disp[end-1] ≈ disp[end] ≈ disp_ref
    end

    q = [0.12, 0.23, 0.34]
    for mode in (:SUN, :dipole), s in (1, 3/2)
        test_biquad(mode, q, s)
    end
end


@testitem "General biquadratic" begin
    using LinearAlgebra

    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    msg = "Nonstandard tetragonal cell for spacegroup 99. Consider `standardize`."
    cryst = @test_logs (:info, msg) Crystal(latvecs, [[0,0,0], [0.4,0,0]]; types=["A", "B"])

    sys = System(cryst, [1 => Moment(s=1, g=2), 2 => Moment(s=2, g=2)], :dipole)
    set_pair_coupling!(sys, (S1, S2) -> +(S1'*diagm([2,-1,-1])*S1)*(S2'*diagm([2,-1,-1])*S2), Bond(1, 2, [0,0,0]))

    θ = randn()
    set_dipole!(sys, [1, 0, 0], (1, 1, 1, 1))
    set_dipole!(sys, [0, cos(θ), sin(θ)], (1, 1, 1, 2))
    energy(sys)
    @test energy(sys) ≈ -3

    swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
    res = intensities_bands(swt, [[0,0,0]])
    @test res.disp[1] ≈ 9
    @test res.data[1] ≈ 1
end


@testitem "Canted AFM" begin

    function test_canted_afm(s)
        J, D, h = 1.0, 0.54, 0.76
        a = 1
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
        positions = [[0, 0, 0]]
        cryst = Crystal(latvecs, positions)
        q = [0.12, 0.23, 0.34]

        sys = System(cryst, [1 => Moment(; s, g=-1)], :dipole)
        sys = reshape_supercell(sys, [1 -1 0; 1 1 0; 0 0 1])
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
        set_field!(sys, [0, 0, h])

        # Numerical
        c₂ = 1 - 1/2s
        θ = acos(h / (2s*(4J+D*c₂)))
        set_dipole!(sys, ( sin(θ), 0, cos(θ)), position_to_site(sys, (0,0,0)))
        set_dipole!(sys, (-sin(θ), 0, cos(θ)), position_to_site(sys, (1,0,0)))
        swt_dip = SpinWaveTheory(sys; measure=nothing)
        ϵq_num = dispersion(swt_dip, [q])

        # Analytical
        Jq = 2J*(cos(2π*q[1])+cos(2π*q[2]))
        ωq₊ = real(√Complex(4J*s*(4J*s+2D*s*c₂*sin(θ)^2) + cos(2θ)*(Jq*s)^2 + 2s*Jq*(4J*s*cos(θ)^2 + c₂*D*s*sin(θ)^2)))
        ωq₋ = real(√Complex(4J*s*(4J*s+2D*s*c₂*sin(θ)^2) + cos(2θ)*(Jq*s)^2 - 2s*Jq*(4J*s*cos(θ)^2 + c₂*D*s*sin(θ)^2)))
        ϵq_ana = [ωq₊, ωq₋]

        ϵq_num ≈ ϵq_ana
    end

    @test test_canted_afm(1)
    @test test_canted_afm(2)
end


@testitem "Local Stevens expansion" begin
    using LinearAlgebra
    a = 1
    latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
    positions = [[0, 0, 0]]
    # P1 point group for most general single-ion anisotropy
    cryst = Crystal(latvecs, positions, 1)

    s = 3
    sys_dip = System(cryst, [1 => Moment(; s, g=-1)], :dipole)
    sys_SUN = System(cryst, [1 => Moment(; s, g=-1)], :SUN)

    # The strengths of single-ion anisotropy (must be negative to favor the dipolar ordering under consideration)
    Ds = -rand(3)
    h  = 0.1*rand()
    M = normalize(rand(3))
    SM = M' * spin_matrices(s)
    aniso = Ds[1]*SM^2 + Ds[2]*SM^4 + Ds[3]*SM^6

    set_onsite_coupling!(sys_dip, aniso, 1)
    set_onsite_coupling!(sys_SUN, aniso, 1)

    set_field!(sys_dip, h*M)
    set_field!(sys_SUN, h*M)

    set_dipole!(sys_dip, M, (1,1,1,1))
    set_dipole!(sys_SUN, M, (1,1,1,1))

    energy(sys_dip)
    energy(sys_SUN)

    swt_dip = SpinWaveTheory(sys_dip; measure=nothing)
    swt_SUN = SpinWaveTheory(sys_SUN; measure=nothing)

    q = rand(3)
    disp_dip = dispersion(swt_dip, [q])
    disp_SUN = dispersion(swt_SUN, [q])

    @test only(disp_dip) ≈ disp_SUN[end-1]
end


@testitem "Dipole-dipole" begin
    latvecs = lattice_vectors(10, 10, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])

    for mode in (:dipole, :SUN)
        sys = System(cryst, [1 => Moment(s=1, g=1)], mode)
        enable_dipole_dipole!(sys, 1.0; demag=0)

        polarize_spins!(sys, (0,0,1))
        @test energy_per_site(sys) ≈ -0.1913132980155851

        swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
        qs = [[0, 0, 0], [0, 0, 1/2], [0, 1/2, 1/2], [0, 0, 0]]
        res = intensities_bands(swt, qs)
        disp_ref = [0.5689399140467553, 0.23914164251944922, 0.23914164251948083, 0.5689399140467553]
        @test isapprox(res.disp[end, :], disp_ref; atol=1e-7)
        @test res.data[end, :] ≈ [2/3, 1, 201/202, 2/3]
    end

    begin
        units = Units(:meV, :angstrom)
        cryst = Sunny.bcc_crystal()
        sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, seed=2)
        enable_dipole_dipole!(sys, units.vacuum_permeability; demag=0)
        polarize_spins!(sys, (1,2,3)) # arbitrary direction

        R = hcat([1,1,-1], [-1,1,1], [1,-1,1]) / 2
        sys_reshape = reshape_supercell(sys, R)
        @test energy_per_site(sys_reshape) ≈ energy_per_site(sys) ≈ -0.89944235377

        swt1 = SpinWaveTheory(sys; measure=nothing)
        swt2 = SpinWaveTheory(sys_reshape; measure=nothing)
        q = [0.5, -0.1, 0.3]
        disp1 = dispersion(swt1, [q])
        disp2 = dispersion(swt2, [q])

        @test disp1 ≈ [1.3236778041378718, 0.9206655245623444]
        @test disp2 ≈ [0.9206655245623366]
    end
end


@testitem "SW15-Langasite" begin
    # Ba3NbFe3Si2O14
    a = b = 8.539
    c = 5.2414
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    crystal = Crystal(latvecs, [[0.24964,0,0.5]], 150)
    sys = System(crystal, [1 => Moment(s=5/2, g=2)], :dipole; seed=5)
    set_exchange!(sys, 0.85,  Bond(3, 2, [1,1,0]))   # J1
    set_exchange!(sys, 0.24,  Bond(1, 3, [0,0,0]))   # J2
    set_exchange!(sys, 0.053, Bond(2, 3, [-1,-1,1])) # J3
    set_exchange!(sys, 0.017, Bond(1, 1, [0,0,1]))   # J4
    set_exchange!(sys, 0.24,  Bond(3, 2, [1,1,1]))   # J5

    for i in 1:3
        θ = -2π*(i-1)/3
        set_dipole!(sys, [cos(θ),sin(θ),0], (1,1,1,i))
    end

    sys = repeat_periodically_as_spiral(sys, (1, 1, 7); k=[0,0,1/7], axis=[0,0,1])

    measure = ssf_custom((q, ssf) -> ssf, sys; apply_g=false)
    swt = SpinWaveTheory(sys; measure)
    q = [0.41568,0.56382,0.76414]
    res = intensities_bands(swt, [q])

    SpinW_energies = [2.6267,2.6541,2.8177,2.8767,3.2458,3.3172,3.4727,3.7767,3.8202,3.8284,3.8749,3.9095,3.9422,3.9730,4.0113,4.0794,4.2785,4.4605,4.6736,4.7564,4.7865]
    SpinW_intensities = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2999830079, -0.2999830079im, 0,0.2999830079im, 0.2999830079, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3591387785, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5954018134, -0.5954018134im, 0,0.5954018134im, 0.5954018134, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.3708506016,1.3708506016im, 0, -1.3708506016im, 1.3708506016, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0511743697, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0734875342, 0.0 + 0.0734875342im, 0, 0.0 - 0.0734875342im, 0.0734875342, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0577275935, -0.0577275935im, 0,0.0577275935im, 0.0577275935, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.1733740706, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0338873034,0.0338873034im, 0, -0.0338873034im, 0.0338873034, 0, 0, 0, 0]

    @test isapprox(res.disp, reverse(SpinW_energies); atol=1e-3)
    @test isapprox(reinterpret(ComplexF64, res.data), SpinW_intensities; atol=1e-7)
end


@testitem "3Q Pyrochlore" begin
    tol = 1e-7
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=-1)], :dipole; dims=(7, 7, 1), seed=0)

    J₁ = -1
    J₃ = 1.6234898018587323
    set_exchange!(sys, J₁, Bond(1, 1, [1,0,0]))
    set_exchange!(sys, J₃, Bond(1, 1, [2,0,0]))

    D = 0.066*J₁
    set_onsite_coupling!(sys, S -> D*S[3]^2, 1)

    h = 2.0*J₃
    set_field!(sys, [0, 0, h])

    dipole_array = zeros(Float64,3,7,7)
    dipole_array[:,:,1]= [0.856769   0.811926  -0.485645   -0.950288   -0.429787  -0.0800935  0.217557;
                          0.35333   -0.2906     0.0685251  -0.297433   -0.854291  -0.004611   0.97169;
                          -0.375638   0.506288   0.871466    0.0921234   0.292353   0.996777   0.0921244]
    dipole_array[:,:,2]= [-0.2422     0.428786   0.365354  -0.0893441   0.0440378   4.84897e-7  -0.524944;
                          0.368806  -0.652928  -0.556337   0.136048   -0.0670579  -7.39508e-7   0.799351;
                         -0.897397  -0.624361   0.746328   0.986665    0.996777    1.0          0.292356]
    dipole_array[:,:,3]= [-0.978133   -0.664683  -0.0559526  0.130092  0.650549   0.954731   0.0360572;
                          0.186465   -0.645827  -0.860547   0.472888  0.753856   0.0549402  0.0716667;
                          0.0921244  -0.375638   0.506288   0.871466  0.0921234  0.292353   0.996777]
    dipole_array[:,:,4]= [-0.0731477  -0.474581  -0.179152  -0.245875  -0.226963   0.440493    0.73273;
                          -0.145399   -0.123781  -0.387346   0.348821   0.898546   0.0253494  -0.674256;
                          0.986665    0.871465   0.90436    0.904361  -0.375636  -0.897398    0.0921208]
    dipole_array[:,:,5]= [0.183479  0.299124  0.773231  0.425028   -0.717278  -0.779847   -0.122394;
                          -0.454847  0.594572  0.381818  0.0385228   0.478729  -0.0448744  -0.918648;
                          0.871464  0.74633   0.506289  0.90436     0.506291  -0.624359   -0.375639]
    dipole_array[:,:,6]= [-0.20424   -0.154295   0.351059   0.891644   0.34449   -0.664477   -0.657631;
                          -0.374726   0.848447   0.697806  -0.252717  -0.349109  -0.0382366  -0.557847;
                          0.90436    0.506291  -0.624359  -0.375639   0.871464   0.74633     0.506289]
    dipole_array[:,:,7]= [0.426641   -0.2224    -0.734377  -0.198296   0.327583   0.162493    0.302168;
                          0.0104859   0.364238   0.565321  -0.394152  -0.940321   0.00934933  0.38632;
                          0.90436     0.904361  -0.375636  -0.897398   0.0921208  0.986665    0.871465]
    for x in 1:7, y in 1:7
        set_dipole!(sys, dipole_array[:,x,y], (x,y,1,1))
    end

    q1 = [0.0391,0.0415,0.5530]
    q2 = [0.2360,0.7492,0.9596]
    q3 = [0.1131,0.7654,0.2810]
    q = [q1,q2,q3]
    formfactors = [1 => FormFactor("Cr4")]
    measure = ssf_custom((q, ssf) -> ssf, sys; formfactors)
    swt = SpinWaveTheory(sys; measure)
    res = intensities_bands(swt, q)
    disp_inds = [107, 89, 118, 140, 112, 16, 103, 75, 142, 18]
    int_inds = [9, 147, 131, 41, 15, 96, 48, 105, 129, 17]
    disp_ref = [8.464621970889235,2.965829202488746,6.539681848582543,2.524276472373584,7.536305768861917,7.21157510322424,9.267100207705882,5.603801899767303,2.2012141464553636,6.933800585478572]
    int_ref = [[1.2252940671236579e-5 + 0.0im 2.274979028773366e-5 - 2.9311941591805888e-6im 1.653212548779866e-5 + 1.1781816967947205e-5im; 2.274979028773366e-5 + 2.9311941591805888e-6im 4.294029257562512e-5 + 0.0im 2.7876377103396445e-5 + 2.582994102312879e-5im; 1.653212548779866e-5 - 1.1781816967947205e-5im 2.7876377103396445e-5 - 2.582994102312879e-5im 3.363456946935e-5 + 0.0im], [0.0023849068155162296 + 0.0im 0.00012123326788175223 + 0.0005597432044950685im 2.9024438353305692e-5 + 0.0019488022511109527im; 0.00012123326788175223 - 0.0005597432044950685im 0.00013753575531155366 + 0.0im 0.0004588638588955403 + 9.225242336242961e-5im; 2.9024438353305692e-5 - 0.0019488022511109527im 0.0004588638588955403 - 9.225242336242961e-5im 0.0015927970884407876 + 0.0im], [3.343948808906739e-5 + 0.0im -0.00020469198105989745 - 4.273264350562476e-5im 0.00014834426130997533 + 1.9072039037346926e-5im; -0.00020469198105989745 + 4.273264350562476e-5im 0.0013075823952430512 + 0.0im -0.0009324269345411454 + 7.282554609603465e-5im; 0.00014834426130997533 - 1.9072039037346926e-5im -0.0009324269345411454 - 7.282554609603465e-5im 0.0006689624696724299 + 0.0im], [0.000371215929258512 + 0.0im 0.00015333075423969357 + 0.0003656727495761781im -0.0005009238108982899 + 0.0003075677944476983im; 0.00015333075423969357 - 0.0003656727495761781im 0.0004235456174856762 + 0.0im 9.606844029120137e-5 + 0.0006204846586781491im; -0.0005009238108982899 - 0.0003075677944476983im 9.606844029120137e-5 - 0.0006204846586781491im 0.0009307860608149384 + 0.0im], [2.9077299015972345e-6 + 0.0im 7.621667422789561e-6 + 2.799505642406002e-6im -2.266990384026427e-6 - 2.148518183351213e-6im; 7.621667422789561e-6 - 2.799505642406002e-6im 2.2673029606106315e-5 + 0.0im -8.010728755108931e-6 - 3.4490269087400155e-6im; -2.266990384026427e-6 + 2.148518183351213e-6im -8.010728755108931e-6 + 3.4490269087400155e-6im 3.3549800413375375e-6 + 0.0im], [4.79789862126743e-5 + 0.0im -2.0507170782740867e-5 - 1.8427040499771927e-5im 1.1851577631841944e-7 + 3.659464580060244e-5im; -2.0507170782740867e-5 + 1.8427040499771927e-5im 1.5842349642893138e-5 + 0.0im -1.4105371891520816e-5 - 1.5595760044615497e-5im; 1.1851577631841944e-7 - 3.659464580060244e-5im -1.4105371891520816e-5 + 1.5595760044615497e-5im 2.7911847518508492e-5 + 0.0im], [0.004511188103971374 + 0.0im -0.0013855170428121722 - 0.003178210621811203im 0.0012245622270242976 + 0.0004875348349735645im; -0.0013855170428121722 + 0.003178210621811203im 0.0026646373317783553 + 0.0im -0.0007195754536780329 + 0.0007129893012763896im; 0.0012245622270242976 - 0.0004875348349735645im -0.0007195754536780329 - 0.0007129893012763896im 0.0003850965251566534 + 0.0im], [1.1287527510791432e-5 + 0.0im 9.784635886704085e-6 + 4.074063598919808e-5im 1.177917483079537e-5 + 1.7619012333705406e-5im; 9.784635886704085e-6 - 4.074063598919808e-5im 0.00015552905794128413 + 0.0im 7.380400215595438e-5 - 2.7242055744441945e-5im; 1.177917483079537e-5 - 1.7619012333705406e-5im 7.380400215595438e-5 + 2.7242055744441945e-5im 3.979423792148188e-5 + 0.0im], [9.457308666169535e-5 + 0.0im -1.7915569710789528e-6 + 2.743855511588096e-5im 2.3814745247182995e-5 - 0.00010529497917934995im; -1.7915569710789528e-6 - 2.743855511588096e-5im 7.994705575514363e-6 + 0.0im -3.100044278968906e-5 - 4.914720059098736e-6im; 2.3814745247182995e-5 + 0.00010529497917934995im -3.100044278968906e-5 + 4.914720059098736e-6im 0.00012322929432616498 + 0.0im], [2.4114075174233513e-6 + 0.0im -2.4607194255187506e-6 - 1.9339332806224738e-6im 2.7253619853740792e-6 + 6.205353429165704e-7im; -2.4607194255187506e-6 + 1.9339332806224738e-6im 4.062041755385664e-6 + 0.0im -3.278759427148752e-6 + 1.5524978029108779e-6im; 2.7253619853740792e-6 - 6.205353429165704e-7im -3.278759427148752e-6 - 1.5524978029108779e-6im 3.2398762990830975e-6 + 0.0im]]
    @test isapprox(res.disp[disp_inds], disp_ref; atol=1e-6)
    @test isapprox(res.data[int_inds], int_ref; atol=1e-7)
end


@testitem "Invariance to reshaping" begin
    # Diamond-cubic with antiferromagnetic exchange
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]], 227; choice="1")
    s = 3/2
    sys = System(cryst, [1 => Moment(; s, g=2)], :dipole; seed=0)
    set_exchange!(sys, 1.0, Bond(1, 3, [0,0,0]))
    randomize_spins!(sys)
    minimize_energy!(sys)
    @test energy_per_site(sys) ≈ -2s^2

    # Reshaped system
    shape = [0 1 1; 1 0 1; 1 1 0] / 2
    sys_prim = reshape_supercell(sys, shape)
    @test energy_per_site(sys_prim) ≈ -2s^2

    # Both systems should produce the same intensities
    formfactors = [1 => FormFactor("Co2")]
    swt1 = SpinWaveTheory(sys_prim; measure=ssf_perp(sys_prim; formfactors))
    swt2 = SpinWaveTheory(sys; measure=ssf_perp(sys; formfactors))
    kernel = lorentzian(fwhm=0.8)
    q = randn(3)
    energies = 0:0.01:6
    res1 = intensities(swt1, [q]; energies, kernel)
    res2 = intensities(swt2, [q]; energies, kernel)
    @test res1.data ≈ res2.data
end


@testitem "Invariance to spin rotation" begin
    using LinearAlgebra, Random

    function build_system(R, D1, D2, J, K1, K2, h, g)
        latvecs = lattice_vectors(1, 1, 1, 92, 93, 94)
        cryst = Crystal(latvecs, [[0,0,0], [0.4,0,0]]; types=["A", "B"])
        moments = [1 => Moment(s=1, g=2), 2 => Moment(s=2, g=R*g*R')]
        sys = System(cryst, moments, :dipole; seed=101)

        set_onsite_coupling!(sys, S -> S'*R*(D1+D1')*R'*S, 1)
        set_onsite_coupling!(sys, S -> S'*R*(D2+D2')*R'*S, 2)

        K1 = Sunny.tracelesspart(K1 + K1')
        K2 = Sunny.tracelesspart(K2 + K2')

        set_pair_coupling!(sys, (S1, S2) -> S1'*R*J*R'*S2 + (S1'*R*K1*R'*S1)*(S2'*R*K2*R'*S2), Bond(1, 2, [0,0,0]))
        set_field!(sys, R*h)

        return sys
    end

    Random.seed!(101)
    g  = randn(3,3)
    D1 = randn(3,3)
    D2 = randn(3,3)
    J  = randn(3,3)
    K1 = randn(3,3)
    K2 = randn(3,3)
    h  = randn(3)
    R1 = Sunny.Mat3(I)
    R2 = Sunny.axis_angle_to_matrix(randn(3), 0.46)
    sys1 = build_system(R1, D1, D2, J, K1, K2, h, g)
    sys2 = build_system(R2, D1, D2, J, K1, K2, h, g)

    randomize_spins!(sys1)
    minimize_energy!(sys1)

    for site in eachsite(sys1)
        sys2.dipoles[site] = R2 * R1' * sys1.dipoles[site]
    end
    @assert energy_per_site(sys1) ≈ energy_per_site(sys2)

    swt = SpinWaveTheory(sys1; measure=ssf_trace(sys1))
    res1 = intensities_bands(swt, [[0,0,0]])
    swt = SpinWaveTheory(sys2; measure=ssf_trace(sys2))
    res2 = intensities_bands(swt, [[0,0,0]])
    @assert res1.data ≈ res2.data
end


@testitem "Generalized interaction consistency" begin
    using LinearAlgebra

    function make_lswt_hamiltonian(sys, q)
        swt = SpinWaveTheory(sys; measure=nothing)
        L = Sunny.nbands(swt)
        H = zeros(ComplexF64, 2L, 2L)
        Sunny.swt_hamiltonian_SUN!(H, swt, Sunny.Vec3(q))
        return H
    end

    sys = System(Sunny.cubic_crystal(), [1 => Moment(s=1, g=1)], :SUN)

    q = [0.23, 0, 0]

    set_pair_coupling!(sys, (Si, Sj) -> -(Si'*Sj), Bond(1,1,[1,0,0]); extract_parts=true)
    H_conventional = make_lswt_hamiltonian(sys, q)

    set_pair_coupling!(sys, (Si, Sj) -> -(Si'*Sj), Bond(1,1,[1,0,0]); extract_parts=false)
    H_alternative = make_lswt_hamiltonian(sys, q)

    @test H_conventional ≈ H_alternative
end


@testitem "Equivalence of dense and sparse Hamiltonian constructions" begin
    using LinearAlgebra

    # Build System and SpinWaveTheory with exchange, field and single-site anisotropy
    function simple_swt(mode)
        cryst = Crystal(diagm([1, 1, 2]), [[0, 0, 0]], "P1")
        sys = System(cryst, [1 => Moment(s=1, g=1)], mode; dims=(8, 1, 1))

        K1 = diagm([2, -1, -1])
        K2 = diagm([-1, -1, 2])
        set_pair_coupling!(sys, (Si, Sj) -> -Si'*Sj + (Si'*K1*Si)*(Sj'*K2*Sj), Bond(1,1,[1,0,0]); extract_parts=true)
        set_onsite_coupling!(sys, S -> S[3]^2, 1)
        set_field!(sys, [0, 0, 0.1])

        randomize_spins!(sys)
        minimize_energy!(sys; maxiters=1_000)

        return SpinWaveTheory(sys; measure=nothing)
    end

    q = Sunny.Vec3(0.5, 0, 0)

    for mode = (:dipole, :SUN)
        # Construct Hamiltonian directly
        swt = simple_swt(mode)
        H1 = Sunny.dynamical_matrix(swt, q)

        # Construct Hamiltonian by sparse matrix-vector multiplies
        L = Sunny.nbands(swt)
        x = Matrix{ComplexF64}(I, 2L, 2L)
        H2 = transpose(Sunny.mul_dynamical_matrix(swt, x, fill(q, 2L)))

        @test H1 ≈ H2
    end
end


@testitem "Spin ladder intensities reference test" begin
    using LinearAlgebra

    function dimer_model(; J=1.0, J′=0.0, h=0.0, dims=(2,1,1), fast=false)
        cryst = Crystal(I(3), [[0,0,0]], 1)
        sys = System(cryst, [1 => Moment(s=3/2, g=1)], :SUN; dims)

        S = spin_matrices(1/2)
        S1, S2 = Sunny.to_product_space(S, S)
        onsite = J*(S1' * S2 - 0.25I(4)) - h*(S1[3] + S2[3])
        set_onsite_coupling!(sys, onsite, 1)

        S1i, S1j = Sunny.to_product_space(S1, S1)
        S2i, S2j = Sunny.to_product_space(S2, S2)
        scalar, bilin, biquad, tensordec = Sunny.decompose_general_coupling(J′*(S1i' * S1j + S2i' * S2j), 4, 4; extract_parts=fast)
        Sunny.set_pair_coupling_aux!(sys, scalar, Sunny.Mat3(I)*bilin, biquad, tensordec, Bond(1, 1, [-1,0,0]))

        return sys, cryst
    end

    # Set up dimer model and observables (0 and π channels manually specified)
    sys, cryst = dimer_model(; J=1.0, J′=0.2, h=0.0, dims=(2,1,1), fast=false)
    s = spin_matrices(1/2)
    S1, S2 = Sunny.to_product_space(s, s)
    observables0 = Hermitian.([
        1/√2*(S1[1] - S2[1]),
        1/√2*(S1[2] - S2[2]),
        1/√2*(S1[3] - S2[3]),
    ])
    observables = repeat(observables0, 1, size(eachsite(sys))...)
    corr_pairs = [(3,3), (2,2), (1,1)]
    combiner = (_, data) -> real(sum(data))
    measure = Sunny.MeasureSpec(observables, corr_pairs, combiner, [one(FormFactor)])

    # Set up SpinWaveTheory
    randomize_spins!(sys)
    minimize_energy!(sys)
    swt = SpinWaveTheory(sys; measure)

    qs = [[0,0,0], [0.5,0,0], [1,0,0]]
    energies = 0:0.5:5
    kernel = lorentzian(fwhm=0.3)
    res = intensities(swt, qs; energies, kernel)
    # println(round.(res.data; digits=12))
    data_ref = [0.042551644188 0.148531187027 0.042551644188; 0.123710785142 0.944407715529 0.123710785142; 1.079575108835 1.261286085876 1.079575108835; 0.492703854423 0.168505531387 0.492703854423; 0.087770506619 0.06066521855 0.087770506619; 0.034461978527 0.030825188879 0.034461978527; 0.018214262744 0.018585357642 0.018214262744; 0.011230027228 0.012410289203 0.011230027228; 0.007607320239 0.008868510563 0.007607320239; 0.005490942587 0.006651305827 0.005490942587; 0.004148615687 0.005172181047 0.004148615687]
    @test isapprox(res.data, data_ref; atol=1e-9)
end


@testitem "LSWT correction to classical energy" begin
    J = 1
    s = 1
    δE_afm1_ref = 0.488056/(2s) * (-2*J*s^2)

    # The results are taken from Phys. Rev. B 102, 220405(R) (2020) for the AFM1
    # phase on the FCC lattice
    function correction(mode)
        a = 1
        latvecs = lattice_vectors(a, a, a, 90, 90, 90)
        positions = [[0, 0, 0]]
        fcc = Crystal(latvecs, positions, 225)
        sys_afm1 = System(fcc, [1 => Moment(; s, g=1)], mode)
        set_exchange!(sys_afm1, J, Bond(1, 2, [0, 0, 0]))
        set_dipole!(sys_afm1, (0, 0,  1), position_to_site(sys_afm1, (0, 0, 0)))
        set_dipole!(sys_afm1, (0, 0, -1), position_to_site(sys_afm1, (1/2, 1/2, 0)))
        set_dipole!(sys_afm1, (0, 0, -1), position_to_site(sys_afm1, (1/2, 0, 1/2)))
        set_dipole!(sys_afm1, (0, 0,  1), position_to_site(sys_afm1, (0, 1/2, 1/2)))
        swt_afm1 = SpinWaveTheory(sys_afm1; measure=nothing)
        # Calculate at low accuracy for faster testing
        δE_afm1 = Sunny.energy_per_site_lswt_correction(swt_afm1; atol=5e-4)
        return isapprox(δE_afm1_ref, δE_afm1; atol=1e-3)
    end

    for mode in (:dipole, :SUN)
        @test correction(mode)
    end
end


@testitem "LSWT correction to the ordered moments (s maximized)" begin
    # Test example 1: The magnetization is maximized to `s`. Reference result
    # comes from Phys. Rev. B 79, 144416 (2009) Eq. (45) for the 120° order on
    # the triangular lattice.
    J = 1
    s = 1/2
    a = 1
    δS_ref = -0.261302

    function δS_triangular(mode)
        latvecs = lattice_vectors(a, a, 10a, 90, 90, 120)
        cryst = Crystal(latvecs, [[0, 0, 0]])
        sys = System(cryst, [1 => Moment(s=s, g=2)], mode)
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        polarize_spins!(sys, [0, 1, 0])
        sys = repeat_periodically_as_spiral(sys, (3, 3, 1); k=[2/3, -1/3, 0], axis=[0, 0, 1])
        swt = SpinWaveTheory(sys; measure=nothing)
        # Calculate first 3 digits for faster testing
        δS = Sunny.magnetization_lswt_correction(swt; atol=1e-3)[1]

        return isapprox(δS_ref, δS, atol=1e-3)
    end

    for mode in (:dipole, :SUN)
        @test δS_triangular(mode)
    end
end


@testitem "LSWT correction to the ordered moments (s not maximized)" begin
    using LinearAlgebra
    # Test example 2: The magnetization is smaller than `s` due to easy-plane
    # single-ion anisotropy The results are derived in the Supplemental
    # Information (Note 12) of Nature Comm. 12.1 (2021): 5331.
    a = b = 8.3193
    c = 5.3348
    lat_vecs = lattice_vectors(a, b, c, 90, 90, 90)
    types = ["Fe"]
    positions = [[0, 0, 0]]
    cryst = Crystal(lat_vecs, positions, 113; types)

    s = 1
    J₁  = 0.266
    J₁′ = 0.1J₁
    Δ = Δ′ = 1/3
    D = 1.42
    gab, gcc = 2.18, 1.93
    g = diagm([gab, gab, gcc])
    x = 1/2 - D/(8*(2J₁+J₁′))

    sys = System(cryst, [1 => Moment(; s, g)], :SUN; dims=(1, 1, 2), seed=0)
    set_exchange!(sys, diagm([J₁, J₁, J₁*Δ]),  Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, diagm([J₁′, J₁′, J₁′*Δ′]), Bond(1, 1, [0, 0, 1]))
    set_onsite_coupling!(sys, S -> D*S[3]^2, 1)

    randomize_spins!(sys)
    minimize_energy!(sys; maxiters=1000)
    swt = SpinWaveTheory(sys; measure=nothing)

    δS = Sunny.magnetization_lswt_correction(swt; atol=1e-2)[1]

    M_cl  = 2*√((1-x)*x)
    # Paper reported M_ref = 2.79, but actual result is closer to 2.78
    M_ref = 2.78
    @test isapprox(M_ref, (M_cl+δS)*√3*gab, atol=1e-2)
end
