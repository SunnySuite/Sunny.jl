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
    @test isapprox(res.disp, disps_golden; atol=1e-7)
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


@testitem "Vacant sites" begin
    # A vacant site must decouple from the spin wave calculation entirely: it
    # contributes only null bands, and the surviving sites behave as though it
    # were absent. `swt_data!` achieves this by zeroing the local frame, which is
    # load bearing for the Zeeman term alone -- the anisotropy and exchange terms
    # are already suppressed by the vanishing spin magnitude.
    cryst = Crystal(lattice_vectors(1, 1.3, 1.7, 90, 90, 90), [[0, 0, 0]])

    function build(mode, dims)
        sys = System(cryst, [1 => Moment(s=1, g=2)], mode; dims, seed=0)
        set_onsite_coupling!(sys, S -> -0.3*S[3]^2, 1)
        set_field!(sys, [0, 0, 0.7])
        randomize_spins!(sys)
        return sys
    end

    for mode in (:dipole, :SUN)
        # An isolated moment, and the same moment obtained by breaking a chain
        # with a vacancy. The lattice is orthorhombic so that the chain bond has
        # no symmetry equivalent wrapping onto the surviving site.
        sys_ref = build(mode, (1, 1, 1))
        minimize_energy!(sys_ref)
        sys = build(mode, (2, 1, 1))
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        sys = to_inhomogeneous(sys)
        set_vacancy_at!(sys, (2, 1, 1, 1))
        minimize_energy!(sys)
        @test 2energy_per_site(sys) ≈ energy_per_site(sys_ref)

        q = [[0.2, 0.3, 0.4]]
        disp_ref = dispersion(SpinWaveTheory(sys_ref; measure=nothing), q)
        disp = dispersion(SpinWaveTheory(sys; measure=nothing), q)
        L = length(disp_ref)
        @test disp[1:L] ≈ disp_ref
        @test all(<(1e-6), disp[L+1:end])
    end
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
    bond = Bond(1, 1, [1, 0, 0])

    q = [0.23, 0, 0]

    set_pair_coupling!(sys, (Si, Sj) -> -(Si'*Sj), bond; extract_parts=true)
    H_conventional = make_lswt_hamiltonian(sys, q)

    msg = "Overwriting coupling for $bond"
    @test_logs (:warn, msg) set_pair_coupling!(sys, (Si, Sj) -> -(Si'*Sj), bond; extract_parts=false)
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
        op = J′*(S1i' * S1j + S2i' * S2j)
        (; scalar, bilin, biquad, tensordec) = Sunny.decompose_general_coupling(op, 4, 4; extract_parts=fast)
        Sunny.set_pair_coupling_aux!(sys, scalar, Sunny.Mat3(I)*bilin, biquad, tensordec, Bond(1, 1, [-1,0,0]), nothing)

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
    # Layout is (nobs × d1 × d2 × d3 × nunits × nparts); here nparts=1.
    operators = repeat(reshape(observables0, 3, 1, 1, 1, 1, 1), 1, size(eachsite(sys))..., 1)
    corr_pairs = [(3,3), (2,2), (1,1)]
    combiner = (_, data) -> real(sum(data))
    formfactors = fill(one(FormFactor), 3, 1, 1)
    measure = Sunny.MeasureSpec(operators, corr_pairs, combiner, formfactors)

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
        # A fixed budget rather than an `rtol`, so that both modes integrate the same
        # function to the same accuracy; `rtol` alone stalls in :SUN mode, whose Nambu
        # space is larger and whose norm the relative test is measured against.
        δE_afm1 = Sunny.energy_per_site_lswt_correction(swt_afm1; maxevals=2000)
        return isapprox(δE_afm1_ref, δE_afm1; atol=1e-3)
    end

    for mode in (:dipole, :SUN)
        @test correction(mode)
    end

    # The onsite coupling contributes a constant at this same order, which
    # `energy_per_site_lswt_correction` now includes. It vanishes identically in
    # `:dipole` mode, where `rcs_factors` leaves the classical energy exact, so only
    # `:dipole_uncorrected` sees it; the reference above is therefore unaffected, having
    # no anisotropy at all. The two modes describe the same model with different
    # truncations of it, so their corrected energies need not agree, only the correction
    # must be present in one and absent in the other.
    function easy_plane(mode)
        cryst = Sunny.square_crystal(; c=3)
        sys = System(cryst, [1 => Moment(s=2, g=1)], mode)
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> 0.3*S[3]^2 + 0.05*S[3]^4, 1)
        sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
        set_dipole!(sys, [+1, 0, 0], (1, 1, 1, 1))
        set_dipole!(sys, [-1, 0, 0], (1, 1, 1, 2))
        minimize_energy!(sys)
        return SpinWaveTheory(sys; measure=nothing)
    end
    @test iszero(Sunny.anisotropy_correction(easy_plane(:dipole)).δE)
    @test Sunny.anisotropy_correction(easy_plane(:dipole_uncorrected)).δE ≈ 0.0542857 rtol=1e-5
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
        δS = Sunny.magnetization_lswt_correction(swt; rtol=1e-3)[1]

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
    # Information (Note 12) of https://doi.org/10.1038/s41467-021-25591-7.
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

    δS = Sunny.magnetization_lswt_correction(swt; rtol=1e-2)[1]

    M_cl  = 2*√((1-x)*x)
    # Paper reported M_ref = 2.79, but actual result is closer to 2.78
    M_ref = 2.78
    @test isapprox(M_ref, (M_cl+δS)*√3*gab, atol=1e-2)
end


@testitem "1/s corrections against exact diagonalization" begin
    using LinearAlgebra, SparseArrays

    # One generic cluster serves every exact-diagonalization check below. Three
    # sites give three bands, enough for off-diagonal elements of Σ̂; a triclinic
    # cell with generic positions leaves no site symmetry, so every Stevens word
    # is allowed; unequal spins exercise the per-site factors σᵢ/σⱼ; and a
    # generic field cants the moments out of collinearity, without which Σ̂ would
    # vanish identically. Every bond has zero offset, so the Hamiltonian is
    # 𝐪-independent and a single grid point integrates the self-energy exactly.
    # Ferromagnetic exchange is what keeps the Fock space affordable, the
    # anomalous mixing being ⟨n̂⟩ ≈ 0.004 here where a canted antiferromagnet of
    # any `s` would put it near 0.3.
    cryst = Crystal(lattice_vectors(1, 1.1, 1.2, 80, 90, 100),
                    [[0, 0, 0], [0.45, 0.05, 0.1], [0.1, 0.4, 0.05]], 1)
    Js = [diagm([-0.6, -0.6, -1.3]), diagm([-1.4, -0.7, -0.7]),
          -0.1*[1 0.3 -0.2; 0.25 1 0.15; -0.15 0.1 1]]
    B = [0.48, -0.32, 0.8]

    function cluster(ss; mode=:dipole, aniso=nothing)
        sys = System(cryst, [i => Moment(s=ss[i], g=1) for i in 1:3], mode)
        for (n, (i, j)) in enumerate([(1, 2), (2, 3), (1, 3)])
            set_exchange!(sys, Js[n], Bond(i, j, [0, 0, 0]))
        end
        set_field!(sys, B)
        isnothing(aniso) || for i in 1:3
            set_onsite_coupling!(sys, aniso(stevens_matrices(mode == :dipole ? ss[i] : Inf), ss[i], i), i)
        end
        # Deterministic, and converged to a torque of 1e-11 from this start
        polarize_spins!(sys, B)
        minimize_energy!(sys; jitter=0)
        return sys
    end

    # Boson operators on a truncated Fock space of `n` sites, labeled by the Nambu
    # index of a `BosonMonomial`: `a ≤ n` annihilates on site `a`, `a > n` creates.
    # Sparse, because the cluster below acts on 7³ states.
    function fock_ops(dims)
        n = length(dims)
        op(O, i) = reduce(kron, (k == i ? O : sparse(1.0I, dims[k], dims[k]) for k in 1:n))
        b(i) = spdiagm(1 => [√float(k) for k in 1:dims[i]-1])
        return a -> ComplexF64.(a <= n ? op(b(a), a) : op(b(a-n)', a-n))
    end

    # Normal-ordered quadratic Hamiltonian, read off from Sunny's own dynamical
    # matrix at 𝐪 = 0, together with that matrix. Any `terms2` are added to it,
    # which is how the anisotropy correction to H₂ enters below.
    function fock_quadratic(swt, bop, dim; terms2=nothing)
        L = Sunny.nbands(swt)
        H = zeros(ComplexF64, 2L, 2L)
        Sunny.dynamical_matrix!(H, swt, zero(Sunny.Vec3))
        isnothing(terms2) || Sunny.accum_quadratic!(H, terms2, zero(Sunny.Vec3))
        H2 = spzeros(ComplexF64, dim, dim)
        for i in 1:L, j in 1:L
            H2 .+= ((H[i, j] + H[L+j, L+i])/2) * bop(L+i)*bop(j) +
                   (H[i, L+j]/2) * bop(L+i)*bop(L+j) + (H[L+i, j]/2) * bop(i)*bop(j)
        end
        return (H, H2)
    end

    expand(bop, terms, dim) = sum(t -> t.c * prod(bop, t.as), terms; init=spzeros(ComplexF64, dim, dim))

    # Exact retarded Green function of the cluster, in the quasi-particle operators
    # y = T⁻¹x = τ₃T†τ₃x, from which the self-energy follows as
    # Σ̂ = ω - diag(ε) - (Gτ₃)⁻¹, returned for each of the `ωs`. A complex frequency
    # keeps every denominator away from a pole, so the comparison is independent of
    # broadening.
    function cluster_self_energy(H2, Hpert, bop, T0, ε, λ, ωs)
        L = size(T0, 1) ÷ 2
        τ₃ = Diagonal([ones(L); -ones(L)])
        Y = [sum(a -> (τ₃ * T0' * τ₃)[m, a] * bop(a), 1:2L) for m in 1:2L]
        (Es, ψs) = eigen(Hermitian(Matrix(H2 + λ*Hpert)))
        ΔE = Es .- Es[1]
        us = [ψs' * (Y[m]' * ψs[:, 1]) for m in 1:2L]   # ⟨0|y_m|j⟩
        vs = [ψs' * (Y[m] * ψs[:, 1]) for m in 1:2L]    # ⟨0|y_m†|j⟩
        return map(ωs) do ω
            G = [sum(@. us[m]*conj(us[m′])/(ω - ΔE) - conj(vs[m])*vs[m′]/(ω + ΔE))
                 for m in 1:2L, m′ in 1:2L]
            return (ω*I - Diagonal(ε) - inv(G * τ₃)) / λ^2
        end
    end


    # ---- The boson expansion against the exact cluster Hamiltonian ----

    # Compares the whole Holstein-Primakoff expansion, order by order in the boson
    # number, against exact diagonalization in the *spin* Hilbert space, so that the
    # truncation of the series is itself what is under test.
    function cluster_errors(ss; mode=:dipole, aniso=nothing)
        sys = cluster(ss; mode, aniso)
        # Regularization would otherwise leak into the quadratic coefficients
        swt = SpinWaveTheory(sys; measure=nothing, regularization=0)
        Ns = ntuple(i -> Int(2ss[i]+1), 3)
        dim = prod(Ns)
        bop = fock_ops(Ns)
        op(O, i) = reduce(kron, (k == i ? sparse(O) : sparse(1.0I, Ns[k], Ns[k]) for k in 1:3))
        S(a, i) = op(spin_matrices(ss[i])[a], i)

        # Exact cluster Hamiltonian, in the local frames that the boson expansion
        # uses and with the same Zeeman convention (+𝐁⋅𝐒). The anisotropy is rotated
        # by hand, so this is independent of the implementation, which instead works
        # from the Stevens coefficients that `swt_data!` stored.
        Rs = swt.data.local_rotations
        Hex = spzeros(ComplexF64, dim, dim)
        for i in 1:3
            Bi = Rs[i]' * (sys.gs[1, 1, 1, i]' * sys.extfield[1, 1, 1, i])
            Hex .+= sum(a -> Bi[a] * S(a, i), 1:3)
            isnothing(aniso) && continue
            A = Hermitian(Matrix(aniso(stevens_matrices(ss[i]), ss[i], i)))
            Hex .+= op(Matrix(Sunny.rotate_operator(A, Rs[i])), i)
        end
        # Each bond appears twice, once culled
        for int in sys.interactions_union, c in int.pair
            c.isculled && continue
            @assert iszero(c.bond.n)
            (i, j) = (c.bond.i, c.bond.j)
            J = Rs[i]' * c.bilin * Rs[j]
            Hex .+= sum(J[a, b] * S(a, i) * S(b, j) for a in 1:3, b in 1:3)
        end

        (; terms2, δE) = Sunny.anisotropy_correction(swt)
        H1 = expand(bop, Sunny.anisotropy_monomials(swt, Val{1}()), dim)
        H3 = expand(bop, Sunny.cubic_monomials(swt), dim)
        H4 = expand(bop, Sunny.quartic_monomials(swt), dim)
        # Quadratic Hamiltonian of LSWT, plus its own correction at order 1/s
        (_, H2) = fock_quadratic(swt, bop, dim; terms2)

        # Flattened Kronecker index of a boson occupation triple, the states of a
        # given total boson number, and the boson numbers of an index
        idx(ms) = 1 + sum(k -> ms[k] * prod(Ns[k+1:3]), 1:3)
        states(n) = [idx(ms) for ms in Iterators.product(0:n, 0:n, 0:n) if sum(ms) == n]
        ms(a) = ntuple(k -> mod(div(a-1, prod(Ns[k+1:3])), Ns[k]), 3)
        (vac, n1, n2) = (idx((0, 0, 0)), states(1), states(2))
        E0 = real(Hex[vac, vac])
        R = Hex - E0*I - H1 - H2 - H3 - H4

        # With no anisotropy the first four blocks below are *exact*, not merely
        # asymptotic in `s`: a matrix element between one and two bosons is purely
        # cubic, since H₄ conserves boson number modulo two, and every term of H₅ is
        # a transverse operator times the classical Sᶻ of its partner, hence
        # proportional to the transverse field that vanishes at a classical minimum.
        return (; classical = abs(E0 - (energy(sys) + 3δE)) / norm(Hex),
                  quadratic = norm(R[n1, n1]) / norm(Hex),
                  anomalous = norm(R[n2, [vac]]) / norm(Hex),
                  cubic = norm(R[n2, n1]) / norm(Hex),
                  # Unlike the blocks above this one does receive an H₆ contribution,
                  # so it is asymptotic even with no anisotropy
                  quartic = norm(R[n2, n2]) / norm(H4[n2, n2]),
                  # Largest residual over every element four bosons can reach, which
                  # is the only available measure once anisotropy makes each word
                  # truncated rather than exact
                  reachable = maximum(abs(R[a, b]) for a in axes(R, 1), b in axes(R, 2)
                                      if all(ms(a) .+ ms(b) .<= 4)) / norm(Hex),
                  coherent = max(abs(3δE), norm(H1)) / norm(Hex),
                  hermiticity = (norm(H3 - H3') + norm(H4 - H4')) / (norm(H3) + norm(H4)))
    end

    # Exchange and field only. Two `s` triples: one certifies every exact block, the
    # second supplies the 1/s scaling of the leftover H₆ piece, which would instead
    # approach unity were the quartic term wrong at O(s⁰).
    rs = map(((1.0, 3/2, 1.0), (2.0, 3.0, 2.0))) do ss
        err = cluster_errors(ss)
        @test err.classical < 1e-12
        @test err.quadratic < 1e-8
        @test err.anomalous < 1e-12
        @test err.cubic < 1e-8
        @test err.hermiticity < 1e-12
        return err.quartic
    end
    @test 0.4 < rs[2] / rs[1] < 0.6

    # Now an onsite anisotropy on every site, whose Stevens coefficients carry a
    # factor s^-k so that the classical energy landscape is held fixed as `s` varies,
    # making the rate at which the residual vanishes meaningful. Here the words are
    # truncated rather than exact — one order in 1/s is kept per word — so what is
    # left over is the next order. Requiring it to vanish faster than 1/s, the size
    # of the retained correction itself, pins every retained coefficient: an error at
    # the order kept would leave a residual of the same size as the correction, and
    # duplicating a term that LSWT already holds would leave one of order unity.
    # Mode :dipole renormalizes the stored Stevens coefficients and
    # :dipole_uncorrected does not, so both are checked.
    aniso(O, s, i) = ((0.3*O[2, 0] + 0.15*O[2, 1])/s^2 + (0.02*O[4, 2] - 0.01*O[4, -3])/s^4 +
                      0.004*O[6, i]/s^6)

    for mode in (:dipole, :dipole_uncorrected)
        # s ≥ 3 is required for the order-6 word to be nonzero
        es = [cluster_errors(f .* (3, 3, 3); mode, aniso) for f in (1, 4/3)]
        @test es[1].hermiticity < 1e-12
        @test es[1].reachable < 0.02
        @test es[2].reachable < 0.5 * es[1].reachable
        # The Stevens coefficients are renormalized in mode :dipole so that the
        # classical energy function is exact in a spin coherent state, to all orders
        # in 1/s. That makes the corrections to the energy and to the linear term
        # vanish identically, and leaves the anomalous coefficient A₂ as the only
        # correction to LSWT's H₂.
        @test mode == :dipole ? es[1].coherent < 1e-12 : es[1].coherent > 1e-3
    end


    # ---- Mean fields, self-energy and pair amplitudes, one Fock space ----

    sys = cluster((1, 3/2, 1))
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)
    nmax = 6
    dim = (nmax+1)^L
    bop = fock_ops(ntuple(_ -> nmax+1, L))
    (H, H2) = fock_quadratic(swt, bop, dim)
    ψ = eigen(Hermitian(Matrix(H2))).vectors[:, 1]
    expect(O) = dot(ψ, O, ψ)

    # The moments cant by up to 20°, which is what makes Σ̂ nonvanishing below
    ds = [normalize(sys.dipoles[1, 1, 1, i]) for i in 1:L]
    @test maximum(norm(ds[i] × ds[j]) for i in 1:L, j in 1:L) ≈ 0.35 atol=0.02

    # Truncation is harmless only if the vacuum has no weight in the top sector
    ψn = reshape(ψ, ntuple(_ -> nmax+1, L))
    @test maximum(i -> norm(selectdim(ψn, i, nmax+1)), 1:L) < 1e-4

    # Every Nambu correlation, including the anomalous ones
    ckeys = [(a, a′, (0, 0, 0)) for a in 1:2L for a′ in 1:2L]
    gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; rtol=1e-8)
    @test gs ≈ [expect(bop(a)*bop(a′)) for (a, a′, _) in ckeys] atol=1e-7

    terms4 = Sunny.quartic_monomials(swt)
    ckeys = Sunny.correlation_keys(L, terms4)
    gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; rtol=1e-8)
    (terms2, δE) = Sunny.hartree_fock_decoupling(terms4, Sunny.correlation_lookup(ckeys, gs, L))
    H4 = expand(bop, terms4, dim)
    H4mf = expand(bop, terms2, dim)

    # Wick's theorem is exact in a Gaussian state, so the constant subtracted by
    # the decoupling is precisely -⟨H₄⟩. This pins the three pairings, their sign,
    # and the fact that the constant is removed once rather than twice.
    @test -δE ≈ expect(H4) atol=1e-7

    # Defining property of the decoupling: the mean-field operator reproduces the
    # response of H₄ to every quadratic perturbation.
    Qs = [bop(a)*bop(a′) for a in 1:2L, a′ in 1:2L]
    @test [expect(H4*Q - Q*H4) for Q in Qs] ≈ [expect(H4mf*Q - Q*H4mf) for Q in Qs] atol=1e-6

    # The analogous decoupling of H₃, which Wick-contracts down to the linear
    # operator that tadpole relaxation must cancel. Only linear perturbations test
    # anything here, since a Gaussian state gives ⟨[H₃, Q]⟩ = 0 for quadratic Q; a
    # commutator of two linear operators is a c-number, so this pins the three
    # cubic pairings exactly rather than to integration accuracy.
    terms3 = Sunny.cubic_monomials(swt)
    ckeys = Sunny.correlation_keys(L, terms3)
    gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; rtol=1e-8)
    ℓ = Sunny.tadpole_vector(terms3, Sunny.correlation_lookup(ckeys, gs, L), L, 1e-6)
    H3raw = expand(bop, terms3, dim)
    H3mf = sum(a -> ℓ[a] * bop(a), 1:2L)
    @test [expect(H3raw*bop(a) - bop(a)*H3raw) for a in 1:2L] ≈
          [expect(H3mf*bop(a) - bop(a)*H3mf) for a in 1:2L] atol=1e-8
    # The linear part must be removed from the perturbation below: normal ordering
    # H₃ in the quasi-particle basis leaves that piece behind, and it is the
    # tadpole correction rather than a loop.
    H3 = H3raw - H3mf

    τ₃ = Diagonal([ones(L); -ones(L)])
    T0 = zeros(ComplexF64, 2L, 2L)
    ε = copy(Sunny.bogoliubov!(T0, H))

    # The whole Nambu matrix of the cubic self-energy, against the exact retarded
    # Green function. This fixes the off-diagonal elements, which mix bands and
    # correct spectral weights, and the anomalous blocks, which admix three magnons
    # into the ground state. Nothing cheaper constrains the off-diagonal ones: a
    # collinear structure has Σ̂ = 0 outright, and the branches of a spiral or an
    # umbrella live in momentum sectors that the cubic vertex cannot mix, so Σ̂ is
    # diagonal there too. They are also the one part of Σ̂ sensitive to the per-band
    # phase that `bogoliubov!` fixes independently; leaving it free corrupts them by
    # 70%, two orders of magnitude above this tolerance.
    ωs = [0.5 + 0.3im, -1.1 + 0.25im]
    Σed = cluster_self_energy(H2, H3, bop, T0, ε, 0.06, ωs)
    Σm = Sunny.cubic_self_energy(swt, [[0, 0, 0]], ωs; η=1e-10, grid=(1, 1, 1))
    @test [Σm[:, :, iω, 1] for iω in eachindex(ωs)] ≈ Σed atol=3e-4
    @test maximum(abs, Σm[:, :, 1, 1] - Diagonal(diag(Σm[:, :, 1, 1]))) >
          0.2 * maximum(abs, diag(Σm[:, :, 1, 1]))

    # Below the three-magnon threshold, unitarity requires τ₃Σ̂ to be Hermitian,
    # which is what makes the Dyson equation preserve spectral weight. The
    # broadening η is what breaks it, by an amount η ∂Σ/∂ω.
    Σh = Sunny.cubic_self_energy(swt, [[0, 0, 0]], [0.5]; η=1e-10, grid=(1, 1, 1))[:, :, 1, 1]
    @test τ₃ * Σh ≈ (τ₃ * Σh)' atol=1e-9

    # Second-order perturbation theory in H₃, evaluated exactly in the truncated
    # Fock space, against the on-shell form of `cubic_self_energy`, which takes no
    # frequencies. The vacuum shift constrains the source channel alone, and the
    # level shifts the diagonal of Σ̂ at the on-shell frequency. Three magnons are
    # created and destroyed in the former, and the unrestricted sum over their
    # bands supplies 3! orderings, cancelling one of the two factors of 3! that
    # relate the symmetrized vertex to Γ₂. Levels are identified by proximity to
    # the harmonic band energy, the Fock spectrum interleaving two-magnon states
    # among the one-magnon ones, so a band permutation would be caught here. The
    # residual is the O(h²) error of the differencing.
    function levels(λ)
        E = eigen(Hermitian(Matrix(H2 + λ*H3))).values
        ΔE = E .- E[1]
        return [E[1]; [ΔE[argmin(abs.(ΔE .- ε[n]))] for n in 1:L]]
    end
    shifts = ((levels(0.01) + levels(-0.01))/2 - levels(0)) / 0.01^2
    Σ = Sunny.cubic_self_energy(swt, [[0, 0, 0]]; η=1e-10, grid=(1, 1, 1))
    @test real(vec(Σ)) ≈ shifts[2:L+1] atol=1e-3
    @test maximum(abs, imag(Σ)) < 1e-8
    U3 = Sunny.vertex(swt, terms3, ntuple(_ -> zero(Sunny.Vec3), 3))
    @test shifts[1] ≈ -6 * sum(abs2(U3[L+n1, L+n2, L+n3]) / (ε[n1] + ε[n2] + ε[n3])
                           for n1 in 1:L, n2 in 1:L, n3 in 1:L) atol=1e-6

    # Both amplitudes for the observable to create a pair of magnons, against exact
    # matrix elements in the same Fock space. This is what pins their relative phase,
    # on which the interference of `corrected_channels` depends and nothing else does:
    # each squared amplitude is invariant under a phase on either one separately, so
    # the two diagonal blocks of the measure — the self-energy and the `direct` channel
    # — are blind to an error here. So are the sum rule and the rotation-invariance
    # checks elsewhere, both being homogeneous in the interference. A sign error was
    # found this way, `pair_amplitude` having returned the amplitude of +b†b where
    # Sᶻ = s - b†b carries a minus.
    #
    # The comparison needs a two-magnon state and the operator that reaches it. For a
    # normalized pair |ab⟩ = y_a†y_b†|0⟩/√(1+δ_ab), each of Sunny's amplitudes is
    # conj(⟨ab|·|·⟩) up to the combinatorial factor √(2/(1+δ_ab)) by which the
    # symmetrized state differs from Sunny's sum over ordered pairs. The readout is
    # the full S^{μν} tensor, so that every component is reachable and nothing is
    # protected by the scalar sum rule that a trace measure enjoys.
    swt2m = SpinWaveTheory(sys; measure=ssf_custom((q, ssf) -> ssf, sys; apply_g=false))
    Nobs = Sunny.num_observables(swt2m.measure)
    u2 = zeros(ComplexF64, 2L, Nobs)
    q2m = [0.23, -0.41, 0.17]
    qr2 = Sunny.to_reshaped_rlu(sys, q2m)
    qg2 = Sunny.orig_crystal(sys).recipvecs * q2m
    Sunny.set_swt_observable_vectors!(u2, swt2m, qr2, qg2)
    pref2 = zeros(ComplexF64, Nobs, L)
    Sunny.pair_amplitude_prefactors!(pref2, swt2m, qr2, qg2)

    # Observable A_ν(q) as a linear form in the Nambu vector. The index swap is forced,
    # not chosen: `intensities_bands` forms Avec[μ] = dot(u[:,μ], T[:,n]) for the left
    # amplitude ⟨0|A†|n⟩, so A = Σ_a u[ā] x_a. The `nambu_correlations` block above
    # fixes the state, so nothing here is free, and the transverse amplitude must come
    # out equal to Sunny's own w = T†u.
    Aodd = [sum(a -> u2[Sunny.nambu_conj(a, L), ν] * bop(a), 1:2L) for ν in 1:Nobs]
    # The longitudinal part, -Σᵢ prefᵢ b†ᵢbᵢ. `pref` carries the conjugated Fourier
    # phase, so it is conjugated back to describe the same operator as `Aodd`.
    Aeven = [-sum(i -> conj(pref2[ν, i]) * bop(L+i) * bop(i), 1:L) for ν in 1:Nobs]

    Y = [sum(a -> (τ₃ * T0' * τ₃)[m, a] * bop(a), 1:2L) for m in 1:2L]
    ψ1 = [Y[m]' * ψ for m in 1:L]
    @test [dot(ψ1[n], Aodd[ν] * ψ) for n in 1:L, ν in 1:Nobs] ≈ (T0' * u2)[1:L, :] atol=1e-7

    U3m = zeros(ComplexF64, 2L, 2L, 2L)
    T3m = [conj(T0[Sunny.nambu_conj(a, L), Sunny.nambu_conj(b, L)]) for a in 1:2L, b in 1:2L]
    Sunny.vertex!(U3m, terms3, ntuple(_ -> zero(Sunny.Vec3), 3), (T0, T0, T3m), similar(U3m))
    w2m = T0' * u2

    # Each amplitude separately, then both at once, which is the quantity the
    # interference actually depends on. A phase convention shared by both would cancel
    # from each separately and survive in the coherent sum, so that sum is compared
    # against first-order perturbation theory for ⟨ab|A_ν|0⟩: the even part reaches the
    # pair directly, the odd part through a one-magnon intermediate.
    #
    # The same matrix elements, accumulated into the measure that `corrected_channels`
    # actually consumes, certify it in its own storage convention rather than as
    # rederived here. Linear bin splitting preserves the zeroth and first moments of
    # the measure exactly, so those two moments pin it without replicating the binning.
    # The 12 block is the interference, which no published calculation constrains.
    ρ2m = Matrix{ComplexF64}[]
    onshell2m = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
    Sunny.accum_pair_measure!(ρ2m, swt2m, terms3, qr2, [zero(Sunny.Vec3)];
                              source_freqs=onshell2m, bin_width=1e-3, pref=pref2)
    (M0, M1) = (sum(ρ2m), sum(((i, r),) -> (i - 1) * 1e-3 * r, enumerate(ρ2m)))
    (E0, E1) = (zeros(ComplexF64, L+Nobs, L+Nobs), zeros(ComplexF64, L+Nobs, L+Nobs))
    for a in 1:L, b in 1:L
        ψpair = (Y[a]' * (Y[b]' * ψ)) / sqrt(1 + (a == b))
        c = sqrt(2 / (1 + (a == b)))
        x = ε[a] + ε[b]
        # Magnon-mediated route: the pair is reached from one magnon through H₃.
        # Fixes the √18 and the external leg.
        med = [dot(ψpair, H3 * ψ1[m]) for m in 1:L]
        @test med ≈ c * [conj(√18 * U3m[a, b, Sunny.nambu_conj(m, L)]) for m in 1:L] atol=1e-6
        # Direct route, including the sign of Sᶻ = s - b†b
        dir = [dot(ψpair, Aeven[ν] * ψ) for ν in 1:Nobs]
        @test dir ≈ c * [conj(Sunny.pair_amplitude(pref2, T0, T0, a, b, ν, L)) for ν in 1:Nobs] atol=1e-6
        # Both routes at once
        mediated(ν) = sum(m -> conj(w2m[m, ν]) * √18 * U3m[a, b, Sunny.nambu_conj(m, L)] / (x - ε[m]), 1:L)
        both = (dir + [sum(m -> med[m] * dot(ψ1[m], Aodd[ν] * ψ) / (x - ε[m]), 1:L) for ν in 1:Nobs]) ./ c
        @test both ≈ [conj(Sunny.pair_amplitude(pref2, T0, T0, a, b, ν, L) + mediated(ν)) for ν in 1:Nobs] atol=1e-6
        y = conj([med; dir]) ./ c
        E0 .+= y * y'
        E1 .+= x .* (y * y')
    end
    @test M0 ≈ E0 atol=1e-6
    @test M1 ≈ E1 atol=1e-5

    # ---- Onsite anisotropy in closed form ----

    # For the words of at most two bosons the correction has a closed form: the
    # leading word that LSWT already holds, times ℓ = -binomial(k, 2)/2s, which is the
    # leading deviation of `rcs_factors` from unity and so vanishes in mode :dipole,
    # plus a further 1/4s on the anomalous A₂ alone. That last piece is present in
    # both modes because LSWT reads A₂ off the classical energy, in effect using a
    # boson coherent state, whose amplitude on two spin deviations exceeds a spin
    # coherent state's by 1/√(1 - 1/2s). Unlike the tests above this is exact rather
    # than asymptotic, so it pins both scalars including their signs. One Stevens
    # order at a time is required, since ℓ depends on k.
    function closed_form_errors(mode, k, s)
        sys = System(Crystal(lattice_vectors(1, 1.3, 1.7, 88, 95, 100), [[0, 0, 0]], 1),
                     [1 => Moment(; s, g=1)], mode)
        O = stevens_matrices(mode == :dipole ? s : Inf)
        set_onsite_coupling!(sys, (0.3*O[k, 0] + 0.2*O[k, 1] - 0.1*O[k, -2])/s^k, 1)
        polarize_spins!(sys, [0.3, 0.5, 0.8])
        minimize_energy!(sys; jitter=0)
        swt = SpinWaveTheory(sys; measure=nothing, regularization=0)

        # LSWT's own coefficients: A₁ of the word b†b and A₂ of the anomalous b†b†
        H = zeros(ComplexF64, 2, 2)
        Sunny.dynamical_matrix!(H, swt, zero(Sunny.Vec3))
        (; terms2, δE) = Sunny.anisotropy_correction(swt)
        coef(as) = sum(t.c for t in terms2 if t.as == as; init=0.0im)

        ℓ = mode == :dipole ? 0 : -binomial(k, 2)/2s
        E = energy_per_site(sys)
        # Errors are measured against a common scale, because in mode :dipole the
        # first two comparisons are between exact zeros
        rel(a, b) = abs(a - b) / max(abs(E), abs(H[1, 1]), abs(H[1, 2]))
        return (; energy = rel(δE, ℓ*E),
                  diagonal = rel(coef((2, 1)), ℓ*H[1, 1]),
                  anomalous = rel(2coef((2, 2)), (ℓ + 1/4s)*H[1, 2]))
    end

    for mode in (:dipole, :dipole_uncorrected), k in (2, 4, 6)
        err = closed_form_errors(mode, k, 4.5)
        @test err.energy < 1e-12
        @test err.diagonal < 1e-12
        @test err.anomalous < 1e-12
    end


    # ---- The vertex at every slot count, with nonzero bond offsets ----

    # Canted square-lattice antiferromagnet, whose bonds connect distinct cells, so
    # that the Fourier phases the cluster above cannot reach are exercised. The
    # tetragonal anisotropy is diagonal in the global frame but not in either local
    # frame, so it contributes to every vertex.
    sys = System(Sunny.square_crystal(; c=3), [1 => Moment(s=2.0, g=1)], :dipole)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    O = stevens_matrices(2)
    set_onsite_coupling!(sys, 0.1*O[4, 0] + 0.05*O[4, 4], 1)
    set_field!(sys, [1.5, 0, 0])
    sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
    set_dipole!(sys, [0, 0, +1], (1, 1, 1, 1))
    set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
    minimize_energy!(sys)
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)

    # Reference contraction: accumulate the symmetrized monomial coefficients into
    # an explicit Nambu-space tensor, then transform every slot at once
    function reference_vertex(terms, qs, Ts)
        K = length(qs)
        W = zeros(ComplexF64, ntuple(_ -> 2L, K))
        perms = Sunny.slot_permutations(K)
        for t in terms, p in perms
            c = t.c / length(perms)
            for slot in 1:K
                c *= cis(2π * dot(qs[slot], t.ns[p[slot]]))
            end
            W[CartesianIndex(ntuple(slot -> t.as[p[slot]], K))] += c
        end
        U = zero(W)
        for n in CartesianIndices(U), a in CartesianIndices(W)
            U[n] += W[a] * prod(slot -> Ts[slot][a[slot], n[slot]], 1:K)
        end
        return U
    end

    q(x, y) = Sunny.Vec3(x, y, 0)
    cases = ((Sunny.cubic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.54, -0.22))),
             (Sunny.quartic_monomials(swt), (q(0.13, 0.29), q(0.41, -0.07), q(-0.22, 0.35), q(-0.32, -0.57))))
    for (terms, qs) in cases
        K = length(qs)
        Ts = Sunny.bogoliubov_matrices(swt, qs)
        buf() = zeros(ComplexF64, ntuple(_ -> 2L, K))
        U = Sunny.vertex!(buf(), terms, qs, Ts)
        @test U ≈ reference_vertex(terms, qs, Ts)

        # Permuting the (momentum, band) slots together leaves the vertex invariant
        for p in Sunny.slot_permutations(K)
            Up = Sunny.vertex!(buf(), terms, ntuple(t -> qs[p[t]], K), ntuple(t -> Ts[p[t]], K))
            @test permutedims(Up, invperm(collect(p))) ≈ U
        end

        # Hermiticity of the underlying operator implies U(-𝐤; n̄) = conj(U(𝐤; n)),
        # where n̄ = n + L exchanges creation and annihilation. Only the magnitudes
        # are compared, because `bogoliubov!` fixes the phase of each band
        # independently at 𝐤 and -𝐤.
        Um = Sunny.vertex(swt, terms, ntuple(t -> -qs[t], K))
        bar(n) = mod1(n + L, 2L)
        @test abs.(U) ≈ [abs(Um[CartesianIndex(map(bar, Tuple(i)))]) for i in CartesianIndices(U)]
    end
end

@testitem "1/s corrections on a lattice" begin
    using LinearAlgebra

    # Square-lattice antiferromagnet, optionally canted by a field. Sunny's Zeeman
    # coupling is +𝐁⋅𝐒, so the moments cant away from the field, with cos θ = -B/8s.
    # At B = 0 the structure is collinear Néel.
    function canted_square(s, B; mode=:dipole)
        sys = System(Sunny.square_crystal(; c=3), [1 => Moment(; s, g=1)], mode)
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        set_field!(sys, [0, 0, B])
        sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
        θ = acos(-B / 8s)
        set_dipole!(sys, [+sin(θ), 0, cos(θ)], (1, 1, 1, 1))
        set_dipole!(sys, [-sin(θ), 0, cos(θ)], (1, 1, 1, 2))
        @assert energy_per_site(sys) ≈ -2s^2 - B^2/16
        return sys
    end

    # Default accuracy for the momentum integrals, enough for the checks below that
    # compare against a reference value. Those that instead compare two corrections
    # computed from the same quadrature, or extract a ratio from them, loosen it
    # individually.
    rtol = 1e-6

    # ---- Oguchi's Z_c, on a supercell and on a multi-atom basis ----

    # For a collinear structure the cubic vertex vanishes, so the mean field is the
    # entire O(1/s) shift of the dispersion, and it is a uniform rescaling by Oguchi's
    # Z_c = 1 + ζ/2s [Prog. Theor. Phys. 13, 148 (1960)]. Both the uniformity in 𝐪 and
    # the 1/s scaling are strong tests of the four-boson coefficients; the numerical
    # values of ζ test their overall normalization. Two values of `s` suffice for the
    # scaling, and ⟨H₄⟩ being of order s⁰ makes the energy correction s-independent.
    # The two lattices differ in how the sublattices arise: a reshaped supercell for
    # the square, and the two atoms of the chemical cell for the honeycomb, whose
    # ζ = 1 - ⟨√(1-|γ_𝐪|²)⟩ with γ_𝐪 = (1 + e^{iq₁} + e^{iq₂})/3 was obtained by
    # direct quadrature.
    function neel_honeycomb(s)
        sys = System(Sunny.hexagonal_crystal(; c=3), [1 => Moment(; s, g=1)], :dipole)
        set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
        set_dipole!(sys, [0, 0, +1], (1, 1, 1, 1))
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
        @assert energy_per_site(sys) ≈ -3s^2/2
        return sys
    end

    qs = [[0.1, 0, 0], [0.3, 0.17, 0], [0.5, 0, 0], [0.13, -0.4, 0.22]]
    for (model, ζ, δE₀, ss) in ((s -> canted_square(s, 0), 0.1579474, 0.01247369, (1/2, 2)),
                                (neel_honeycomb, 0.20984170, 0.01651258, (1/2,)))
        for s in ss
            swt = SpinWaveTheory(model(s); measure=nothing)
            (; terms2, δE) = Sunny.hartree_fock_correction(swt; rtol)
            Zc = Sunny.corrected_dispersion(swt, qs, terms2) ./ dispersion(swt, qs)
            @test maximum(abs, Zc .- Zc[1]) < 1e-6
            @test 2s * (Zc[1] - 1) ≈ ζ atol=1e-6
            @test δE ≈ δE₀ atol=1e-7
        end
    end

    # ---- What a collinear structure switches off ----

    # The onsite correlation ⟨b†ᵢbᵢ⟩ must reproduce Sunny's independent calculation of
    # the moment reduction, and the commutator ⟨bᵢb†ᵢ⟩ - ⟨b†ᵢbᵢ⟩ must come out to one.
    # Self-consistency is inert for this collinear antiferromagnet, where the mean field
    # merely rescales H₂ and so leaves the Bogoliubov transformation, hence the mean
    # fields themselves, unchanged. And every cubic monomial carries a transverse
    # component of the exchange in the local frame, so the cubic vertex, the self-energy
    # and the tadpole all vanish identically, the last because the cubic monomials are
    # proportional to a transverse effective field that vanishes at a classical minimum.
    swt = SpinWaveTheory(canted_square(1/2, 0); measure=nothing)
    L = Sunny.nbands(swt)
    gs = Sunny.nambu_correlations(swt, [(L+1, 1, (0, 0, 0)), (1, L+1, (0, 0, 0))],
                                  Sunny.BosonMonomial{2}[]; rtol=1e-4)
    @test real(gs[1]) ≈ -Sunny.magnetization_lswt_correction(swt; rtol=1e-4)[1] atol=1e-5
    @test gs[2] - gs[1] ≈ 1 atol=1e-10
    @test maximum(abs, Sunny.cubic_self_energy(swt, [[0.3, 0.1, 0]]; η=0.01, grid=(6, 6, 1))) < 1e-12
    Zcs = map(maxiters -> Sunny.corrected_dispersion(swt, [[0.3, 0.1, 0]],
                  Sunny.hartree_fock_correction(swt; maxiters, tol=1e-9, rtol=1e-4).terms2), (1, 100))
    @test Zcs[1] ≈ Zcs[2] atol=1e-7
    tad = Sunny.tadpole_correction(swt; rtol)
    @test all(t -> abs(t.c) < 1e-12, tad.terms2)
    @test abs(tad.δE) < 1e-12
    @test tad.dipoles ≈ [[1, 0, 0], [-1, 0, 0]] / 2

    # ---- Two independent routes to the tadpole, canted ----

    # Thermodynamic consistency, and the observable correction, both on the canted
    # structure and both of the form "two routes agree to leading order, so their
    # difference falls off like 1/s".
    #
    # The correction to the uniform magnetization can be read from `corrected_dipoles`,
    # which combines the tadpole tilt with the reduction of the moment magnitude, or
    # obtained as ∂/∂B of the zero-point energy. Those routes share no machinery, and the
    # second is s-independent. The order in which the two corrections compose is itself a
    # 1/s² choice: shortening along the tilted axis, as `corrected_dipoles` does, and
    # shortening along the original one differ by 1.6% at s = 1 and converge at the same
    # rate, the former slightly faster.
    #
    # The tilt also corrects the amplitude for creating one magnon, at relative order
    # 1/s. Here `observable_corrections` displaces the boson as b → b + v within the
    # untilted frame, whereas re-expressing the transverse spin components about the
    # tilted axis rotates the observable vectors. The rotation is exact in the tilt
    # angle while the displacement is linear in it, so their difference must be smaller
    # than the correction itself by one more power of the tilt. Since the local frames
    # make `v` complex, this is sensitive to the conjugations and to the Nambu labeling.
    function canted_routes(s)
        # Differencing the zero-point energy loses accuracy, so 1e-5 is as tight as the
        # magnetization can be read anyway; the tilt comparison is a ratio of two
        # corrections from the same quadrature, so it is looser still
        B = 3s
        sys = canted_square(s, B)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        tad = Sunny.tadpole_correction(swt; rtol=1e-5)
        dipoles = [sys.dipoles[1, 1, 1, i] for i in 1:2]
        # `corrected_dipoles` applies the tilt and the shortening together, which is what
        # makes it comparable to ∂/∂B below; either alone would be an incomplete 1/s
        corrected = Sunny.corrected_dipoles(swt; rtol=1e-5)
        δmz = sum(i -> (corrected[i] - dipoles[i])[3], 1:2) / 2
        zp(B′) = Sunny.energy_per_site_lswt_correction(
            SpinWaveTheory(canted_square(s, B′); measure=nothing); rtol=1e-5)
        δzp = (zp(B + 1e-4) - zp(B - 1e-4)) / 2e-4

        δc = Sunny.observable_corrections(swt; v=tad.v, rtol=1e-4) -
             Sunny.observable_corrections(swt; rtol=1e-4)
        (err, mag) = (0.0, 0.0)
        for i in 1:2
            R = swt.data.local_rotations[i]
            σ = √2 * swt.data.sqrtS[i]
            Rot = Sunny.rotation_between_vectors([0, 0, 1], R' * tad.dipoles[i])
            for μ in 1:3
                (O, O′) = (swt.data.observables[μ, i], Rot' * swt.data.observables[μ, i])
                δcp = ((O′[1] - O[1]) - im*(O′[2] - O[2])) / 2
                δcm = ((O′[1] - O[1]) + im*(O′[2] - O[2])) / 2
                err = max(err, abs(σ*δcp - δc[i, μ]), abs(σ*δcm - δc[2+i, μ]))
                mag = max(mag, abs(δc[i, μ]), abs(δc[2+i, μ]))
            end
        end
        return (δmz, δzp, err / mag)
    end
    rs = map(canted_routes, (1, 2))
    @test rs[1][1] ≈ 0.0485824 atol=1e-6
    @test rs[1][2] ≈ rs[2][2] atol=1e-8
    @test rs[1][1] - rs[1][2] ≈ 2 * (rs[2][1] - rs[2][2]) rtol=0.03
    @test rs[1][3] < 0.005
    @test rs[1][3] / rs[2][3] ≈ 2 rtol=0.02

    # The tilt half of `corrected_dipoles` needs the cubic vertex, which is unavailable in
    # :SUN mode, whereas the shortening half is not. Rather than error, it must then
    # return the dipoles shortened and untilted, which is what it also returns whenever
    # the tilt is zero. Both quantities are keyword-required, an `rtol` or a `maxevals`
    # being the only control over the cubature.
    let
        sys = canted_square(1, 3; mode=:SUN)
        swt = SpinWaveTheory(sys; measure=nothing)
        dipoles = [sys.dipoles[1, 1, 1, i] for i in 1:2]
        δS = Sunny.magnetization_lswt_correction(swt; rtol=1e-5)
        corrected = Sunny.corrected_dipoles(swt; rtol=1e-5)
        @test corrected ≈ [(1 + δS[i]/norm(dipoles[i])) * dipoles[i] for i in 1:2] rtol=1e-6
        # The tilt that :dipole mode does supply is real, so the two must differ
        @test !isapprox(corrected, Sunny.corrected_dipoles(
            SpinWaveTheory(canted_square(1, 3); measure=nothing); rtol=1e-5), rtol=1e-3)
        @test_throws ErrorException Sunny.corrected_dipoles(swt)
        @test_throws ErrorException Sunny.energy_per_site_lswt_correction(swt)
    end

    # ---- The tadpole sourced by an anisotropy alone ----

    # An onsite anisotropy sources the tadpole on its own, because the quantum
    # correction δE₀ that `anisotropy_correction` makes to the classical energy depends
    # on the direction of the moment. The linear monomial that `tadpole_correction` adds
    # to its source ℓ must therefore be the gradient of δE₀. Displacing the boson by v
    # tilts the moment to (Sˣ, Sʸ) = σ(Re v, Im v) in the local frame, and a linear term
    # c b† + h.c. contributes 2 Re(c v̄), so the gradient is 2(Re c, Im c)/σ. The
    # identity relates the words of one operator in a rotated frame, so it holds
    # separately at every order in 1/s; both sides are read off at the same order,
    # making it exact rather than asymptotic.
    function anisotropic_site(n)
        cryst = Crystal(lattice_vectors(1, 1.1, 1.3, 88, 92, 95), [[0, 0, 0]], 1)
        sys = System(cryst, [1 => Moment(s=3, g=1)], :dipole_uncorrected)
        O = stevens_matrices(Inf)
        set_onsite_coupling!(sys, 0.3*O[2, 0] - 0.2*O[2, -1] + 0.02*O[4, 2], 1)
        set_dipole!(sys, n, (1, 1, 1, 1))
        return SpinWaveTheory(sys; measure=nothing)
    end
    swt = anisotropic_site([0.3, 0.5, 0.8])
    Rloc = swt.data.local_rotations[1]
    δE₀(n) = Sunny.anisotropy_correction(anisotropic_site(n)).δE
    # Tilting n by t along a transverse axis of the local frame moves that component of
    # the dipole by s t, to first order
    grad = [(δE₀(normalize(Rloc[:, 3] + 1e-5*Rloc[:, k])) - δE₀(normalize(Rloc[:, 3] - 1e-5*Rloc[:, k]))) / (2e-5 * 3)
            for k in 1:2]
    terms1 = Sunny.anisotropy_monomials(swt, Val{1}())
    cb = only(t.c for t in terms1 if t.as == (2,))  # coefficient of b†
    @test only(t.c for t in terms1 if t.as == (1,)) ≈ conj(cb)
    @test grad ≈ 2 * [real(cb), imag(cb)] / √6 rtol=1e-6

    # ---- The mean-field operator, and its two consumers ----

    # On the canted structure, whose bonds connect distinct cells, `accum_quadratic!`
    # must agree with the phase convention of `vertex!`, which the exact-diagonalization
    # tests check independently. Contracting a quadratic monomial list at momenta
    # (𝐪, -𝐪) gives U₂ with Σ U₂[n₁,n₂] y_𝐪[n₁] y_{-𝐪}[n₂], whereas `accum_quadratic!`
    # produces H with (1/2) x†_𝐪 H x_𝐪 = (1/2) y†_𝐪 T†HT y_𝐪. Using y_𝐪[m]† = y_{-𝐪}[m̄],
    # the two agree once the latter is symmetrized over its slots. The Bogoliubov matrix
    # at -𝐪 is built from the one at 𝐪, via T_{-𝐪}[a,m] = conj(T_𝐪[ā,m̄]), because
    # `bogoliubov!` fixes the phase of each band independently and only a consistent
    # pair of matrices can be compared.
    swt = SpinWaveTheory(canted_square(1, 3); measure=nothing)
    L = Sunny.nbands(swt)
    terms2 = Sunny.hartree_fock_correction(swt; rtol).terms2
    bar(m) = mod1(m + L, 2L)
    function mean_field_matrix(q)
        H = zeros(ComplexF64, 2L, 2L)
        Sunny.accum_quadratic!(H, terms2, q)
        return H
    end
    let q = Sunny.Vec3(0.23, -0.41, 0.17)
        Hq = zeros(ComplexF64, 2L, 2L)
        T = zeros(ComplexF64, 2L, 2L)
        Sunny.dynamical_matrix!(Hq, swt, q)
        Sunny.bogoliubov!(T, Hq)
        Tm = [conj(T[bar(a), bar(m)]) for a in 1:2L, m in 1:2L]
        U = Sunny.vertex!(zeros(ComplexF64, 2L, 2L), terms2, (q, -q), (T, Tm))
        A = T' * mean_field_matrix(q) * T / 2
        B = Tm' * mean_field_matrix(-q) * Tm / 2
        @test U ≈ [(A[bar(n2), n1] + B[bar(n1), n2])/2 for n1 in 1:2L, n2 in 1:2L]
    end

    # Iterating the mean fields to self-consistency must reach a fixed point that does
    # not depend on how strongly the iteration is damped. Both paths share the same
    # quadrature, so its accuracy is irrelevant here and a loose `rtol` keeps the
    # twenty-odd iterations cheap.
    ress = map((0.0, 0.5)) do damping
        r = Sunny.hartree_fock_correction(swt; maxiters=100, tol=1e-9, damping, rtol=1e-3)
        return (r.δE, Sunny.corrected_dispersion(swt, [[0.3, 0.1, 0]], r.terms2))
    end
    @test ress[1][1] ≈ ress[2][1] atol=1e-8
    @test ress[1][2] ≈ ress[2][2] atol=1e-7

    # `static_self_energy` is the term linear in the correction of the shift that
    # `corrected_dispersion` obtains by rediagonalizing.
    scaled(λ) = [Sunny.BosonMonomial(λ*t.c, t.as, t.ns) for t in terms2]
    q = [[0.3, 0.1, 0]]
    δdisp = (Sunny.corrected_dispersion(swt, q, scaled(1e-4)) - Sunny.corrected_dispersion(swt, q, scaled(-1e-4))) / 2e-4
    @test δdisp ≈ Sunny.static_self_energy(swt, q, terms2) atol=1e-6

    # ---- Goldstone modes survive every correction ----

    # Near a protected zero mode each correction separately diverges like 1/ε, the
    # Bogoliubov matrix doing so, which makes their cancellation stringent; it is the
    # product ε × shift that must vanish. Two symmetries are gated, and they switch on
    # different corrections.
    #
    # Rotation about the field axis leaves the canted structure a zero mode at 𝐪 = 0,
    # and there all three of mean field, tadpole and cubic self-energy contribute. What
    # remains is the discretization error of the self-energy integral, falling off
    # like 1/nk.
    t2 = [terms2; Sunny.tadpole_correction(swt; rtol).terms2]
    δ = Sunny.static_self_energy(swt, [[0, 0, 0]], t2)[2]
    @test δ > 1e3
    rs = map(nk -> (δ + real(Sunny.cubic_self_energy(swt, [[0, 0, 0]]; η=0.005, grid=(nk, nk, 1))[2])) / δ, (16, 32))
    @test rs[1] ≈ 2 * rs[2] rtol=0.01
    @test rs[2] < 0.02

    # The second symmetry is broken only by an onsite anisotropy. Every term below is
    # invariant under rotation about ẑ, so the in-plane moment of this easy-plane
    # ferromagnet leaves a zero mode at 𝐪 = 0, while the cubic and linear vertices
    # vanish identically by that same symmetry. That leaves `anisotropy_correction` to
    # cancel the mean field on its own -- which it can only do because `anisotropy_words`
    # keeps exactly one order in 1/s per word. Since a perturbation of the quadratic
    # form opens a gap like its square root, expanding the anisotropy exactly instead,
    # thereby injecting a partial set of O(1/s²) terms, would gap the mode at O(1/s).
    for mode in (:dipole, :dipole_uncorrected)
        s = 2.0
        sys = System(Sunny.square_crystal(; c=3), [1 => Moment(; s, g=1)], mode)
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        O = stevens_matrices(mode == :dipole ? s : Inf)
        set_onsite_coupling!(sys, (0.5*O[2, 0] - 0.2*O[4, 0]/s^2) / (3s^2), 1)
        set_dipole!(sys, [1, 0, 0], (1, 1, 1, 1))
        sw = SpinWaveTheory(sys; measure=nothing)
        @test maximum(abs(t.c) for t in Sunny.cubic_monomials(sw)) < 1e-12
        @test maximum(abs(t.c) for t in Sunny.anisotropy_monomials(sw, Val{1}()); init=0.0) < 1e-12
        ε = dispersion(sw, [[0, 0, 0]])[1]
        mf = Sunny.hartree_fock_correction(sw; rtol).terms2
        resid(t2) = Sunny.static_self_energy(sw, [[0, 0, 0]], t2)[1] * ε
        @test abs(resid([mf; Sunny.anisotropy_correction(sw).terms2])) < 1e-7 < abs(resid(mf))
    end
end

@testitem "1/s corrections to intensities" begin
    using LinearAlgebra

    # Triangular-lattice antiferromagnet in a three-site cell, small enough for the
    # cubic self-energy to be affordable. Isotropic and in zero field, this is the
    # s = 1/2 model of arXiv:0901.4803 in its 120° state. Given an easy-plane
    # anisotropy and a tilted field it instead cants into a state with three
    # inequivalent sublattices, which is what makes the off-diagonal elements of Σ̂
    # large. The state is built explicitly rather than minimized, so that the
    # chirality of the former and the degenerate direction of the latter are fixed.
    cryst = Sunny.triangular_crystal(; a=1.0, c=10.0)
    Q = cryst.recipvecs * [1/3, 1/3, 0]
    ds = [[0.120272, -0.449974, -0.181818], [-0.465766, -0.002090, -0.181818],
          [0.112161, 0.452064, -0.181818]]
    function triangular(; Δ=1.0, field=nothing, φ=0, g=2)
        sys = System(cryst, [1 => Moment(s=1/2, g=g)], :dipole)
        set_exchange!(sys, diagm([1.0, 1.0, Δ]), Bond(1, 1, [1, 0, 0]))
        sys = reshape_supercell(sys, [2 -1 0; 1 1 0; 0 0 1])
        R = [cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1]
        if isnothing(field)
            for site in eachsite(sys)
                θ = dot(Q, global_positions(sys)[site])
                set_dipole!(sys, [cos(θ), sin(θ), 0], site)
            end
            @assert energy_per_site(sys) ≈ -1.5 * (1/2)^2
        else
            set_field!(sys, R * field)
            for i in 1:3
                set_dipole!(sys, R * ds[i], (1, 1, 1, i))
            end
            minimize_energy!(sys)
            @assert energy_per_site(sys) ≈ -0.51131313 atol=1e-8
        end
        return sys
    end

    sys = triangular()
    swt = SpinWaveTheory(sys; measure=nothing)
    L = Sunny.nbands(swt)

    # ---- Harmonic results, and the folding convention ----

    # The energy per site is -0.5388 J, and the magnon energy at the M point of the
    # original lattice is 2Js, the lowest of the three folded bands. Both Γ and K fold
    # onto 𝐪 = 0, where the harmonic dispersion vanishes, so all three bands are
    # gapless there. Over the whole zone the dispersion is Eq. (11) of arXiv:1306.1231,
    # written in the reciprocal lattice units of the original one-site cell; the
    # three-site cell folds 𝐪 together with 𝐪 ± 𝐊, so each wavevector gates all three
    # bands at once, and with them the folding convention.
    @test energy_per_site(sys) + Sunny.energy_per_site_lswt_correction(swt; rtol=1e-4) ≈ -0.53881 atol=1e-4
    q = [[1/2, 0, 0]]
    @test dispersion(swt, q)[:] ≈ [√2.5, √2.5, 1] atol=1e-6
    γ(q) = (cos(2π*q[1]) + cos(2π*q[2]) + cos(2π*(q[1] + q[2]))) / 3
    ε11(q) = 3 * (1/2) * sqrt(max(0, (1 - γ(q)) * (1 + 2γ(q))))
    K = [1/3, 1/3, 0]
    @test all([[h, k, 0] for h in 0.1:0.4:0.9, k in 0.1:0.4:0.9]) do q
        isapprox(sort(dispersion(swt, [q])[:]), sort([ε11(q + n*K) for n in -1:1]); atol=1e-6)
    end

    # ---- The cubic self-energy on a lattice ----

    # Being a symmetric energy minimum, the 120° structure cannot be tilted by
    # zero-point fluctuations. Unlike the collinear case the cubic monomials are
    # individually nonzero, so their cancellation here tests their relative phases.
    tad = Sunny.tadpole_correction(swt; rtol=1e-3)
    @test maximum(t -> abs(t.c), tad.terms2) < 1e-6
    @test abs(tad.δE) < 1e-12

    terms2 = [Sunny.hartree_fock_correction(swt; rtol=1e-3).terms2; tad.terms2]
    δ = Sunny.static_self_energy(swt, q, terms2)[:]
    Σs = map(nk -> Sunny.cubic_self_energy(swt, q; η=0.02, grid=(nk, nk, 1))[:], (24, 48))

    # The self-energy converges like 1/nk, so a Richardson step gives the O(1/s) magnon
    # energy at the M point. It falls 27% below the harmonic value, most of that coming
    # from the cubic self-energy rather than the mean-field shift.
    εs = [1 + δ[3] + real(Σ[3]) for Σ in Σs]
    @test 2εs[2] - εs[1] ≈ 0.7316 atol=2e-3

    # Because the harmonic dispersion vanishes at Γ and K, the lower edge of the
    # two-magnon continuum touches the one-magnon branch at every wavevector, and a
    # magnon acquires a width only where the branch lies strictly inside the continuum.
    # There the width survives η → 0, as it does for the top of the band here, where
    # 2Γ/ε extrapolates to about 0.2, of the order of the maximum ~0.3 that the
    # reference reports. The M-point magnon instead sits on the boundary, and its
    # apparent width is entirely the Lorentzian tail of the regularization, falling off
    # like η.
    Σ = Sunny.cubic_self_energy(swt, q; η=0.01, grid=(48, 48, 1))[:]
    @test imag(Σ[3]) ≈ imag(Σs[2][3]) / 2 rtol=0.01
    @test imag(Σ[1]) / imag(Σs[2][1]) > 0.85

    # The loop grid must keep both internal lines, 𝐩 and 𝐪-𝐩, off the zone centre,
    # where the cubic vertex diverges. Offsetting by half a step does that only for 𝐩:
    # 𝐪 = [0, 1/4, 0] is [1/4, 1/4, 0] in the reshaped cell, so at nk = 26 the
    # reflected grid hits the zone centre exactly, and one point out of 26² then
    # dominates the integral. Nothing about this 𝐪 is singular, so grids on either side
    # of it must agree; a fixed half-step offset instead gave +1.07 - 0.23im for the
    # first band, wrong even in sign, and -2.82 for the second.
    Σgrid = [Sunny.cubic_self_energy(swt, [[0, 1/4, 0]]; η=0.02, grid=(nk, nk, 1))[:] for nk in (24, 26)]
    @test all(Σ -> isapprox(Σ, [-0.655 - 0.086im, -0.996 - 0.040im, -0.996 - 0.040im]; atol=0.012), Σgrid)

    # Binning the decay measure in the pair energy is a choice of quadrature, not a
    # change of interface, so it must reproduce the frequency loop it replaces. The
    # error is second order in the bin width relative to the regulator Γ, which here is
    # carried by the imaginary part of the frequencies. The reference is written out
    # term by term, straight from the formula at the head of SelfEnergy.jl, so that
    # nothing but `foreach_cubic_line` is shared with the implementation under test.
    terms3 = Sunny.cubic_monomials(swt)
    ps = Sunny.loop_wavevectors((12, 12, 1))
    ε = dispersion(swt, q)[:]
    onshell = [(ε[m] + ε[m′])/2 for m in 1:L, m′ in 1:L]
    ωs = range(0, 2, 6) .+ im*0.06
    k = Sunny.to_reshaped_rlu(sys, q[1])
    Σloop = let Σ = zeros(ComplexF64, L, L, length(ωs))
        Sunny.foreach_cubic_line(swt, terms3, k, ps, L) do a, _b, u, x, _T1, _T2
            for m′ in 1:L, m in 1:L
                R = 18 * conj(u[m]) * u[m′]
                # The source channel a > L is frozen at its on-shell frequency, hence
                # is ω-independent and enters with the opposite sign
                for iω in eachindex(ωs)
                    Σ[m, m′, iω] += a > L ? -R / (onshell[m, m′] - x) : R / (ωs[iω] - x)
                end
            end
        end
        Σ ./ length(ps)
    end
    Σbin = let ρ = Matrix{ComplexF64}[]
        Σsrc = Sunny.accum_pair_measure!(ρ, swt, terms3, k, ps; source_freqs=onshell, bin_width=0.06/16)
        Sunny.pair_self_energy!(zeros(ComplexF64, L, L, length(ωs)), ρ, Σsrc, ωs, 0.06/16)
    end
    @test maximum(abs, Σbin - Σloop) / maximum(abs, Σloop) < 1e-3

    # ---- The corrected spectral function ----

    # The corrected structure factor is a spectral function in its own right, not
    # merely one to the order worked to. Because the Dyson equation is solved in the
    # particle block, with the source channel of the self-energy frozen on shell, its
    # denominator has imaginary part at least the regulator η, so the intensity is
    # non-negative and no taller than a resolution-limited peak of the same weight; and
    # because that denominator grows as ωI, the transverse weight of each 𝐪 is exactly
    # the static weight of the corrected observables. All three properties are violated
    # at s = 1/2 by inverting the full Nambu denominator instead: near 𝐪 = [0.476, 0, 0]
    # a mirror pole is pushed up through ω = 0, giving intensities of -2.8 and +3.6
    # against a bound of 2.3, and at 𝐪 = [0.375, 0.125, 0] the weight comes out 27% low.
    # The momentum-space integrals need no great accuracy here: the identities hold for
    # any self-energy and any observable amplitudes. The options below are those that
    # `tol` selects, so the reference weight is built from the same mean fields as the
    # spectrum.
    opts = (; rtol=0.01, maxevals=100_000)
    swt2 = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
    Ncells = Sunny.nsites(sys) / Sunny.natoms(cryst)
    qs2 = [[0.476, 0, 0], [0.375, 0.125, 0]]

    # Static weight of the one-magnon bands, Σ_{n ≤ L} w[n, μ] conj(w[n, ν]) with
    # w = T'ũ, contracted through the measure's combiner exactly as
    # `intensities_corrected` contracts the spectral function. Passing `drop` subtracts
    # |δA|², leaving the cross term 2Re(A conj(δA)) that is the correction of relative
    # order 1/s; the sum rule below wants that, whereas the spectral weight wants the
    # full corrected amplitude. This is used again for the square lattice further down.
    function static_weights(swt, qs, δc; drop=false)
        sys = swt.sys
        cryst = Sunny.orig_crystal(sys)
        L = Sunny.nbands(swt)
        T = zeros(ComplexF64, 2L, 2L)
        H = zeros(ComplexF64, 2L, 2L)
        u = zeros(ComplexF64, 2L, Sunny.num_observables(swt.measure))
        δu = zero(u)
        Ncells = Sunny.nsites(sys) / Sunny.natoms(cryst)
        return map(qs) do q
            q_reshaped = Sunny.to_reshaped_rlu(sys, q)
            q_global = cryst.recipvecs * q
            Sunny.excitations!(T, H, swt, q)
            Sunny.set_swt_observable_vectors!(u, swt, q_reshaped, q_global)
            fill!(δu, 0)
            isnothing(δc) || Sunny.accum_observable_corrections!(δu, swt, q_reshaped, q_global, δc)
            w = T' * (u + δu)
            δw = T' * δu
            corr = map(swt.measure.corr_pairs) do (μ, ν)
                c = dot(view(w, 1:L, μ), view(w, 1:L, ν))
                drop && (c -= dot(view(δw, 1:L, μ), view(δw, 1:L, ν)))
                c / Ncells
            end
            real(swt.measure.combiner(q_global, corr))
        end
    end
    δc = Sunny.observable_corrections(swt2; v=Sunny.tadpole_correction(swt2; opts...).v, opts...)
    refs = static_weights(swt2, qs2, δc)

    # A Lorentzian tail needs range rather than resolution, so the window is wide and
    # the step is a fraction of η. The residual 0.17% is the truncated tail.
    η = 0.06
    energies = range(-20, 24, 1501)
    chans = Sunny.corrected_channels(swt2, qs2; energies, η, tol=opts.rtol,
                                     loop_grid=(12, 12, 1), mean_field_maxevals=opts.maxevals)
    # The two weight identities are properties of the transverse spectral function, so
    # that channel is taken on its own; the pair channel created directly by the
    # observable carries weight of its own, and their interference sums to zero only
    # over the whole zone.
    (; transverse) = chans
    @test all(≥(0), transverse + chans.cross + chans.direct)
    @test all(vec(maximum(transverse; dims=1)) .< refs ./ (π * η))
    @test vec(sum(transverse; dims=1)) * step(energies) ≈ refs rtol=5e-3

    # The interference is invisible to a trace measure integrated over the zone: 𝐒·𝐒 is
    # a scalar, so its expansion in bosons has no term linking an odd number of them to
    # an even one, and the cancellation is exact for every ω once Σ_𝐪 restores momentum
    # conservation. It is the only check here that constrains the *relative phase* of
    # the two routes to a pair, a phase that each channel alone is free of. The window
    # must reach below ω = 0, or a pair at x ≈ 0 contributes just half of its
    # Lorentzian.
    qs4 = vec([[i, j, 0] ./ 3 for i in 0:2, j in 0:2])
    chans4 = Sunny.corrected_channels(swt2, qs4; energies=range(-4, 12, 161), η=0.3,
                                      loop_grid=(3, 3, 1), mean_field_maxevals=opts.maxevals)
    @test abs(sum(chans4.cross)) < 1e-3 * sum(chans4.direct)

    # Reference for the `direct` channel of `corrected_channels`: the same sum over pairs
    # of magnons, but broadening each pair individually instead of binning its energy
    # first, and contracting the three observable amplitudes β by `contract` rather than
    # through a `measure`. The binning is therefore the only approximation separating the
    # two, so this doubles as the gate on `bin_width`, which `corrected_channels` derives
    # from `η` rather than exposing. Too slow for a converged grid, but an identity holds
    # grid by grid.
    function direct_unbinned(contract, swt, qs, energies, η, grid)
        (; sys) = swt
        L = Sunny.nbands(swt)
        Ncells = Sunny.nsites(sys) / Sunny.natoms(Sunny.orig_crystal(sys))
        ref = zeros(length(energies), length(qs))
        for (iq, q) in enumerate(qs)
            q_reshaped = Sunny.to_reshaped_rlu(sys, Sunny.Vec3(q))
            pref = zeros(ComplexF64, 3, L)
            q_global = Sunny.orig_crystal(sys).recipvecs * Sunny.Vec3(q)
            Sunny.pair_amplitude_prefactors!(pref, swt, q_reshaped, q_global)
            ps = Sunny.loop_wavevectors(grid, q_reshaped)
            Sunny.foreach_magnon_pair(swt, q_reshaped, ps) do _p, T1, T2, ε1, ε2
                for b in 1:L, a in 1:L
                    β = ntuple(μ -> Sunny.pair_amplitude(pref, T1, T2, a, b, μ, L), 3)
                    x = ε1[a] + ε2[b]
                    for (iω, ω) in enumerate(energies)
                        ref[iω, iq] += contract(β) * (η/π) / ((ω - x)^2 + η^2)
                    end
                end
            end
            view(ref, :, iq) ./= Ncells * length(ps)
        end
        return ref
    end

    # Orientation of the μν pair in the channels of `corrected_channels`. Every combiner
    # used above (`ssf_trace`, `ssf_perp`) puts zero weight on off-diagonal `corr_pairs`
    # and on imaginary parts, so a μν transpose — which for a Hermitian S is a complex
    # conjugation — is invisible to all of them, as are the rotation and sum-rule checks,
    # the latter being homogeneous in the interference. This measure sees it, Im Sˣʸ being
    # genuinely nonzero on the 120° structure, so `direct` is fixed in sign and not just
    # in magnitude. A transpose would deviate by 2 rather than by the 3e-4 of the binning.
    let
        swt5 = SpinWaveTheory(sys; measure=ssf_custom((q, ssf) -> imag(ssf[1, 2]), sys; apply_g=false))
        (η, grid) = (0.15, (8, 8, 1))
        qs = [[0.3, 0.2, 0], [1/6, 1/6, 0]]
        energies = range(0, 4, 61)
        direct = Sunny.corrected_channels(swt5, qs; energies, η, loop_grid=grid).direct
        ref = direct_unbinned(β -> imag(β[1] * conj(β[2])), swt5, qs, energies, η, grid)
        scale = maximum(abs, ref)
        @test scale > 1e-2   # the measure is not trivially zero
        @test maximum(abs, direct - ref) < 1e-3 * scale
    end

    # ---- Gauge invariance, where Σ̂ is genuinely off-diagonal ----

    # With three inequivalent sublattices the off-diagonal elements of Σ̂ are a third of
    # the diagonal, whereas for the 120° structure above they vanish identically: its
    # branches sit in momentum sectors that the cubic vertex cannot connect, as do those
    # of the umbrella that a field along z produces. Only the exact-diagonalization
    # cluster elsewhere in this file constrains them otherwise. They are also the one
    # part of Σ̂ sensitive to the per-band phase that `bogoliubov!` fixes independently,
    # since a gauge conjugation Σ̂ → DΣ̂D† with D diagonal and unitary leaves the
    # eigenvalues of the Dyson denominator alone, and with them the poles, and leaves
    # τ₃Σ̂ Hermitian; only the contraction against the observable amplitudes T'u sees it.
    # Rotation about z is an exact symmetry of this model, fixing both diagm([1, 1, Δ])
    # and the trace structure factor while moving every local frame, so the corrected
    # intensities must be invariant under it. Breaking that gauge violates this by 5%.
    # The loop grid is coarse because an invariance holds grid by grid and needs no
    # converged integral.
    swts = map(φ -> let sys3 = triangular(; Δ=0.6, field=[0.7, 0, 1.2], φ, g=1)
                        SpinWaveTheory(sys3; measure=ssf_trace(sys3; apply_g=false))
                    end, (0.0, 0.9))
    qs3 = [[0.23, 0.11, 0], [0.37, 0.09, 0]]
    Σ3 = Sunny.cubic_self_energy(swts[1], qs3[1:1], dispersion(swts[1], qs3[1:1])[1:1];
                                 η=0.05, grid=(6, 6, 1))[1:L, 1:L, 1, 1]
    @test maximum(abs, Σ3 - Diagonal(diag(Σ3))) > 0.2 * maximum(abs, diag(Σ3))
    datas = map(swt3 -> Sunny.intensities_corrected(swt3, qs3; energies=range(0.2, 2.0, 25),
                                                    η=0.15, loop_grid=(4, 4, 1)).data, swts)
    @test datas[1] ≈ datas[2] rtol=1e-6

    # ---- Quantum sum rule on the square lattice ----

    # A field-polarized ferromagnet conserves Sᶻ. Its Bogoliubov transformation is
    # trivial, so the two-magnon channel must carry no weight at all, and nothing else
    # in the 1/s expansion is nonzero for this state either: every mean field vanishes
    # in the empty vacuum, and a collinear structure has no cubic vertex. So the
    # corrected intensities must reduce to those of linear spin wave theory, which is
    # itself exact here, the polarized state and its one-magnon excitations being
    # eigenstates.
    square = Sunny.square_crystal(; c=3)
    let
        sys = System(square, [1 => Moment(s=1, g=2)], :dipole)
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        set_field!(sys, [0, 0, 0.5])
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 1))
        @assert energy_per_site(sys) ≈ -2 - 1.0
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        qs = [[0.3, 0.2, 0]]
        (energies, η) = (range(0, 10, 101), 0.2)
        chans = Sunny.corrected_channels(swt, qs; energies, η, loop_grid=(4, 4, 1), mean_field_maxevals=1000)
        @test maximum(abs, chans.direct) < 1e-25
        @test chans.transverse + chans.cross + chans.direct ≈
              intensities(swt, qs; energies, kernel=lorentzian(fwhm=2η)).data atol=1e-12
    end

    # Easy-axis Néel order on the square lattice, whose gap makes every momentum
    # integral below converge exponentially.
    function square_afm(; field)
        sys = System(square, [1 => Moment(s=1.0, g=1)], :dipole_uncorrected)
        set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> -0.5*S[3]^2, 1)
        set_field!(sys, field)
        sys = reshape_supercell(sys, [1 1 0; 1 -1 0; 0 0 1])
        set_dipole!(sys, [0, 0, +1], (1, 1, 1, 1))
        set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
        minimize_energy!(sys)
        return sys
    end

    # The fast path of `corrected_channels`. A collinear structure has no cubic vertex,
    # but its monomials cancel on merging rather than being absent, so the gate must be on
    # magnitude: the coefficients of the 26 monomials here sum to 1e-12 of the quadratic
    # Hamiltonian, against 5e-2 once a field cants the same model. Its observable
    # consequence is that `cross` comes out identically zero rather than at the 1e-24 of
    # the round-off it replaces, the rest of the result being untouched — checked to be
    # bit-identical, which the two tests below cannot see.
    let
        # Largest cubic coefficient relative to the quadratic Hamiltonian, the
        # dimensionless quantity the gate thresholds
        function vertex_scale(swt)
            L = Sunny.nbands(swt)
            H = zeros(ComplexF64, 2L, 2L)
            Sunny.dynamical_matrix!(H, swt, Sunny.Vec3(0, 0, 0))
            terms3 = Sunny.cubic_monomials(swt)
            @test length(terms3) == 26
            return maximum(abs(t.c) for t in terms3) / norm(H)
        end

        sys = square_afm(; field=[0, 0, 0])
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys))
        @test vertex_scale(swt) < 1e-11
        canted = SpinWaveTheory(square_afm(; field=[1.5, 0, 0]); measure=nothing)
        @test vertex_scale(canted) > 1e-2

        chans = Sunny.corrected_channels(swt, [[0.3, 0.2, 0]]; energies=range(0, 12, 121),
                                         η=0.2, loop_grid=(8, 8, 1), mean_field_maxevals=1000)
        @test iszero(chans.cross) && !iszero(chans.direct)
    end

    # Weights of the three channels into which the quantum sum rule decomposes, all per
    # site and in units where a trace measure is used, so that no local frame projection
    # survives.
    #
    # Because Sᶻ = s - b†b is exact in the local frame, the two-magnon spectrum
    # saturates the longitudinal sum rule ⟨(δSᶻ)²⟩ = n(1+n) + |Δ|², where n = ⟨b†b⟩ and
    # Δ = ⟨bb⟩ follow from Wick's theorem. The transverse channel obeys an identity
    # sharper still: because the truncated S⁺ = σ(b - b†bb/4s) makes S⁻S⁺ = n̂(2s+1-n̂)
    # exact, and because completeness turns a sum of one-magnon weights over bands and
    # wavevectors into the static expectation value ⟨Â†Â⟩, the one-magnon bands must
    # carry ⟨(Sˣ)² + (Sʸ)²⟩ = s + 2s⟨n̂⟩ - ⟨n̂²⟩. Linear spin wave theory produces only
    # the first two terms; the -⟨n̂²⟩ is supplied entirely by `observable_corrections`,
    # so this fixes both the sign and the magnitude of that correction.
    function channel_weights(sys; nq=16, rtol=1e-6)
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        L = Sunny.nbands(swt)

        # Onsite ⟨b†b⟩ and ⟨bb⟩, from which ⟨n̂²⟩ = ⟨n̂⟩² + ⟨n̂⟩(1+⟨n̂⟩) + |Δ|²
        ckeys = [[(L+i, i, (0, 0, 0)) for i in 1:L]; [(i, i, (0, 0, 0)) for i in 1:L]]
        gs = Sunny.nambu_correlations(swt, ckeys, Sunny.BosonMonomial{2}[]; rtol)
        ss = [swt.data.sqrtS[i]^2 for i in 1:L]
        n = real.(gs[1:L])
        n2 = @. n^2 + n * (1 + n) + abs2(gs[L+1:2L])

        # A uniform grid offset by half a step cancels the phase factors of all
        # correlations at distances below `nq`, leaving only the onsite ones above.
        qs = vec([[(a - 0.5)/nq, (b - 0.5)/nq, 0] for a in 1:nq, b in 1:nq])
        δc = Sunny.observable_corrections(swt; rtol)
        harm = sum(static_weights(swt, qs, nothing)) / length(qs)
        transverse = sum(static_weights(swt, qs, δc; drop=true)) / length(qs)

        # Elastic weight of the ordered moment, shortened by zero-point fluctuations
        δm = Sunny.magnetization_lswt_correction(swt; rtol)
        elastic = sum(i -> (ss[i] + δm[i])^2, 1:L) / L

        # Two-magnon continuum, integrated over energy. Its wavevector average converges
        # quickly enough to use a coarser grid, which matters because each point requires
        # its own momentum-space integral. For a one-atom chemical cell this average is
        # the longitudinal weight per site.
        qs2 = vec([[(a - 0.5)/3, (b - 0.5)/3, 0] for a in 1:3, b in 1:3])
        energies = range(-2, 16, 181)
        direct = Sunny.corrected_channels(swt, qs2; energies, η=0.2, loop_grid=(12, 12, 1),
                                          mean_field_maxevals=1000).direct
        longitudinal = sum(direct) * step(energies) / length(qs2)

        return (; harm, transverse, elastic, longitudinal,
                harm_ref = sum(@. ss + 2ss*n) / L,
                transverse_ref = sum(@. ss + 2ss*n - n2) / L,
                casimir = sum(@. ss * (ss + 1)) / L)
    end

    for field in ([0, 0, 0], [1.5, 0, 0])
        # Canting makes the onsite ⟨bb⟩ nonzero, exercising the anomalous contraction
        w = channel_weights(square_afm(; field))
        @test abs(w.harm / w.harm_ref - 1) < 1e-7
        @test abs(w.transverse / w.transverse_ref - 1) < 1e-7

        # The quantum sum rule, and the point of the whole exercise. Because 𝐒⋅𝐒 is a
        # Casimir, the elastic weight of the ordered moment, the one-magnon bands and the
        # two-magnon continuum must together carry exactly s(s+1), and each is produced by
        # a different part of this module. Linear spin wave theory saturates the rule only
        # through O(s): using its uncorrected one-magnon weights instead overshoots by
        # ⟨n̂²⟩, which is the entire O(s⁰) content of the rule, and some 3% of s(s+1) here.
        # The residual error is that of the energy integral above.
        @test abs((w.elastic + w.transverse + w.longitudinal) / w.casimir - 1) < 1e-3
        @test (w.elastic + w.harm + w.longitudinal) / w.casimir - 1 > 0.02
    end

    # The sum rules above constrain the total weight of the continuum, which converges
    # exponentially here, but not its distribution in energy. This gates the shape, and
    # with it the two discretizations that produce it. Measured against a 64×64 grid, the
    # error falls as 3.5e-2, 4.5e-3, 2.8e-4 at 8, 16 and 32 points per direction, so the
    # first figure below is a convergence rate; the second isolates the binning of the
    # pair energy, which is an order of magnitude smaller than the grid error it is paired
    # with and so never the limiting approximation.
    let
        sys = square_afm(; field=[1.5, 0, 0])
        swt = SpinWaveTheory(sys; measure=ssf_trace(sys; apply_g=false))
        qs = [[0.3, 0.2, 0]]
        (energies, η) = (range(-2, 16, 181), 0.2)
        direct(grid) = Sunny.corrected_channels(swt, qs; energies, η, loop_grid=grid,
                                                mean_field_maxevals=1000).direct
        ref = direct((32, 32, 1))
        scale = maximum(abs, ref)
        @test 1e-3 < maximum(abs, direct((16, 16, 1)) - ref) / scale < 1e-2
        unbinned = direct_unbinned(β -> sum(abs2, β), swt, qs, energies, η, (32, 32, 1))
        @test maximum(abs, unbinned - ref) < 3e-4 * scale
    end
end
