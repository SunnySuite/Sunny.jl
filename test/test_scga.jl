@testitem "diamond_lattice" begin
    using LinearAlgebra

    # test against JuliaSCGA (S. Gao)
    dia = [0.7046602277469309, 0.8230846832863896, 0.23309034250417973, 0.40975668535137943, 0.8474163786642979, 0.8230846832694241, 0.723491683211756, 0.5939752161027589, 0.6506966347286152, 0.8012263819500781, 0.23309034265153963, 0.5939752161512792, 0.7779185770415442, 0.9619923476121188, 0.28363795492206234, 0.4097566857133061, 0.6506966347914551, 0.9619923474164132, 0.8576708385848646, 0.45534457001475764, 0.8474163779976567, 0.8012263817424181, 0.2836379554342603, 0.4553445702621035, 0.9352400940660723]
    a = 8.5031 # (Å)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]], 227; choice="1")
    s = 3/2
    sys = System(cryst, [1 => Moment(; s, g=1)], :dipole; seed=0)
    # sys = reshape_supercell(sys, primitive_cell(cryst)) FIXME
    set_exchange!(sys, -1/s^2, Bond(1, 3, [0, 0, 0]))
    set_exchange!(sys, 1/(2s)^2, Bond(1, 2, [0, 0, 0]))
    measure = ssf_perp(sys)
    scga = Sunny.SCGA(sys; measure, Nq=8, sublattice_resolved=false)
    kT = 15*meV_per_K
    γ = s^2 * Sunny.natoms(cryst)
    grid = q_space_grid(cryst, [1, 0, 0], 0:0.9:4, [0, 1, 0], 0:0.9:4; orthogonalize=true)
    res = intensities_static(scga, grid; kT)
    @test isapprox(vec(res.data)/γ, dia; rtol=1e-8)
end

@testitem "square_lattice" begin
    # test against JuliaSCGA (S. Gao)
    tol = 1e-5
    sq = [0.3578208506624103, 0.3850668587212228, 0.4211633125595451, 0.3930558742921638, 0.3586715762816389, 0.3775147329089877, 0.4188431376589759, 0.4009835864744404, 0.36119385315852837, 0.3850668587212228, 0.4063051066815142, 0.4182331099724885,
          0.36751592492720886, 0.32829927921801794, 0.3490454303292574, 0.4062747596231383, 0.41626772462487194, 0.38789161791308147, 0.4211633125595451, 0.4182331099724885, 0.3710074695764304, 0.2924827325221696, 0.253584799014828, 0.2732638980307635,
          0.34406882609120737, 0.409719043169102, 0.4213848627804445, 0.3930558742921638, 0.36751592492720886, 0.2924827325221696, 0.21946942640699554, 0.18901809131848027, 0.2041521166894061, 0.26469888441402106, 0.3466451094806684, 0.3903985131607334,
          0.3586715762816389, 0.32829927921801794, 0.25358479901482806, 0.18901809131848027, 0.16311271576100106, 0.17594081183630647, 0.2284557904965832, 0.3060310966855855, 0.35526653426084637, 0.3775147329089877, 0.3490454303292574, 0.2732638980307635,
          0.20415211668940608, 0.17594081183630647, 0.18993269166531307, 0.2466353934125701, 0.32714470142342356, 0.37442377142308597, 0.4188431376589759, 0.4062747596231383, 0.34406882609120737, 0.26469888441402106, 0.2284557904965832, 0.2466353934125701,
          0.3154328042516787, 0.3918720231667212, 0.4179005987958521, 0.4009835864744404, 0.41626772462487194, 0.409719043169102, 0.3466451094806684, 0.3060310966855855, 0.32714470142342356, 0.3918720231667212, 0.42107414035570556, 0.40321330759757756,
          0.36119385315852837, 0.38789161791308147, 0.4213848627804445, 0.3903985131607334, 0.35526653426084637, 0.37442377142308597, 0.4179005987958521, 0.40321330759757756, 0.3645203324273335];
    a = 1 # (Å)
    latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, [1 => Moment(; s=1, g=1)], :dipole; seed=0)
    set_exchange!(sys, -1, Bond(1, 1, [1,0,0]))
    set_exchange!(sys, 0.5, Bond(1, 1, [1,1,0]))
    set_exchange!(sys, 0.25, Bond(1, 1, [2,0,0]))
    measure = ssf_perp(sys)
    scga = Sunny.SCGA(sys; measure, Nq=8, sublattice_resolved=false)
    kT = 27.5*meV_per_K
    grid = q_space_grid(cryst, [1, 0, 0], -1:0.12:0, [0, 1, 0], -1:0.12:0; orthogonalize=true)
    res = Sunny.intensities_static(scga, grid; kT)
    @test isapprox(vec(res.data)/2, sq; rtol=1e-5)
end


@testitem "MgCr2O4" begin
    # test against Conlon and Chalker
    latvecs    = lattice_vectors(8.3342, 8.3342, 8.3342, 90, 90, 90)
    positions  = [[0.1250, 0.1250, 0.1250],
                [0.5000, 0.5000, 0.5000],
                [0.2607, 0.2607, 0.2607]]
    types      = ["Mg","Cr","O"]
    spacegroup = 227 # Space Group Number
    xtal_mgcro = Crystal(latvecs, positions, spacegroup; types)
    cryst = subcrystal(xtal_mgcro,"Cr")
    sys = System(cryst, [1 => Moment(; s=3/2, g=1)], :dipole) # Same on MgCr2O4 crystal
    J1      = 3.27  # value of J1 in meV from Bai's PRL paper
    J_mgcro = [1.00,0.0815,0.1050,0.0085]*J1; # further neighbor pyrochlore relevant for MgCr2O4
    set_exchange!(sys, J_mgcro[1], Bond(1, 2, [0,0,0]))  # J1
    set_exchange!(sys, J_mgcro[2], Bond(1, 7, [0,0,0]))  # J2
    set_exchange!(sys, J_mgcro[3], Bond(1, 3, [1,0,0]))  # J3a -- Careful here!
    set_exchange!(sys, J_mgcro[4], Bond(1, 3, [0,0,0])); # J3b -- And here!
    # values from Conlon + Chalker
    MgCrO=[ 5.7767400e-02   7.4075451e-02   1.5971443e-01   5.5634087e-01   1.3136038e+00   5.5634087e-01   1.5971443e-01   7.4075451e-02   5.7767400e-02   7.4075451e-02   1.5971443e-01   5.5634087e-01   1.3136038e+00   5.5634087e-01   1.5971443e-01   7.4075451e-02   5.7767400e-02   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   5.7767400e-02   7.4075451e-02   1.5971443e-01   5.5634087e-01   1.3136038e+00   5.5634087e-01   1.5971443e-01   7.4075451e-02   5.7767400e-02   7.4075451e-02   1.5971443e-01   5.5634087e-01   1.3136038e+00   5.5634087e-01   1.5971443e-01   7.4075451e-02   5.7767400e-02   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   1.8741195e+00   3.2690074e+00   1.8741195e+00   1.3136038e+00   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.2377330e+00   2.8136620e+00   2.4168819e+00   1.8741195e+00   2.4168819e+00   2.8136620e+00   1.2377330e+00   5.5634087e-01   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   6.1505984e-01   1.7143609e+00   2.8136620e+00   3.2690074e+00   2.8136620e+00   1.7143609e+00   6.1505984e-01   1.5971443e-01   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   1.6242284e-01   6.1505984e-01   1.2377330e+00   1.8741195e+00   1.2377330e+00   6.1505984e-01   1.6242284e-01   7.4075451e-02   5.7767400e-02   7.4075451e-02   1.5971443e-01   5.5634087e-01   1.3136038e+00   5.5634087e-01   1.5971443e-01   7.4075451e-02   5.7767400e-02   7.4075451e-02   1.5971443e-01   5.5634087e-01   1.3136038e+00   5.5634087e-01   1.5971443e-01   7.4075451e-02   5.7767400e-02]
    measure = ssf_custom((q, ssf) -> real(sum(ssf)), sys)
    scga = Sunny.SCGA(sys; measure, Nq=8, quantum_sum_rule=true, sublattice_resolved=false)
    kT = 20*meV_per_K
    grid = q_space_grid(cryst, [1, 0, 0], range(-4, 4, 17), [0, 1, 0], range(-4, 4, 17); orthogonalize=true)
    res = Sunny.intensities_static(scga, grid; kT)
    S = 3/2
    γ=S*(S+1)*length(cryst.positions)
    @test isapprox(vec(res.data)/γ, (3/4)*vec(MgCrO); rtol=1e-3)
    # factor 3/4 comes from the fact that C+C solve for a single spin component and
    # have a four site unit cell.
end

@testitem "Arbitrary Anisotropy" begin
    # test against JuliaSCGA (S. Gao)
    arb=[0.7586033244771696, 0.7615797236771995, 0.7067376202481463, 0.6713652267703079, 0.7037910523546276, 0.7542251266797739, 0.7483886564702621, 0.7289040202049699, 0.7594392143289115, 0.7113287543280767, 0.6214180238797973, 0.5790787072798553, 0.6192096683970346, 0.7067376202481463, 0.7574706875426616, 0.7617944185477006, 0.704270642625653, 0.6205167495949255,
    0.5219040226511134, 0.4819720409900168, 0.5213611759220986, 0.6202962239550798, 0.7083751254473052, 0.7350050888005822, 0.6689048865527003, 0.5778569918362061, 0.48158130944975375, 0.4442454028677924, 0.48203511793705517, 0.5800585042005479, 0.6760080974704061, 0.7093415594909234, 0.700944539713289, 0.6175127864675868, 0.5205893697530085, 0.48172879530096047,
    0.5219040226511135, 0.6218544838522482, 0.7110417061581527, 0.738269833048121, 0.75085044977932, 0.704270642625653, 0.6188902360706803, 0.5792753532751317, 0.6214319191875379, 0.7113287543280766, 0.763936847880277, 0.7690196643410564, 0.7450119108505591, 0.7544416648799901, 0.7062484206003281, 0.6745924555605828, 0.7100415800949798, 0.7633912278503215,
    0.7586033244771694, 0.7392096662818042, 0.7256714709772182, 0.7586014794304904, 0.7324922896074505, 0.7074913262740532, 0.7368083541302771, 0.768006438327862, 0.7387669446218719, 0.7113287543280766]
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]],1)
    sys = System(cryst, [1 => Moment(; s=1, g=1)], :dipole; seed=0)
    set_exchange!(sys, -1, Bond(1, 1, [1,0,0]))
    set_exchange!(sys, -1, Bond(1, 1, [0,1,0]))
    set_exchange!(sys, 0.5, Bond(1, 1, [1,1,0]))
    set_exchange!(sys, 0.5, Bond(1, 1, [-1,1,0]))
    set_exchange!(sys, 0.25, Bond(1, 1, [2,0,0]))
    set_exchange!(sys, 0.25, Bond(1, 1, [0,2,0]))
    to_inhomogeneous(sys)
    anis = [0.23 0.56 0.34;
            0.56 0.12 0.45;
            0.34 0.45 0.67]
    anis = 0.5.*(anis + anis')
    set_onsite_coupling!(sys, S -> S'*anis*S, 1)
    k_grid = [0.125:0.125:1;]
    measure = ssf_perp(sys)
    scga = Sunny.SCGA(sys; measure, Nq=8, sublattice_resolved=false)
    kT = 55*meV_per_K
    grid = q_space_grid(cryst, [1, 0, 0], 0.125:0.125:1, [0, 1, 0], 0.125:0.125:1; orthogonalize=true)
    res = Sunny.intensities_static(scga, grid; kT)
    @test isapprox(vec(res.data), arb; rtol=1e-6)
end

@testitem "Ferrimagnetic chain sum rule" begin
    tol = 1e-2
    latvecs    = lattice_vectors(3, 5, 8, 90, 90, 90)
    positions  = [[0,0,0],
                [0.5, 0,0]]
    types = ["Ni2","Fe3"]
    cryst = Crystal(latvecs, positions;types)
    s1 = 1
    s2 = 5/2
    moments = [1 => Moment(; s=s1, g=1),2 => Moment(; s=s2, g=1)]
    sys = System(cryst, moments, :dipole)
    J1 = 1
    set_exchange!(sys, J1, Bond(1, 2, [0,0,0]))
    measure = ssf_trace(sys; apply_g=false)
    scga = Sunny.SCGA(sys; measure, sublattice_resolved=true, Nq=60)
    kT = 22.5*meV_per_K
    qarray = range(0,2,60)
    qs1 = [[qx, 0, 0] for qx in qarray]
    res = Sunny.intensities_static(scga, qs1; kT, λs_init=[7.67033234814451, 1.2272531757030904])
    sum_rule = s1^2 + s2^2
    @test abs(sum(res.data)/length(qs1)-sum_rule )/sum_rule < tol
end

@testitem "Ferrimagnetic chain" begin
    tol = 5e-2
    latvecs = lattice_vectors(3, 5, 8, 90, 90, 90)
    positions = [[0, 0, 0], [0.5, 0, 0]]
    types = ["Ni2","Fe3"]
    cryst = Crystal(latvecs, positions; types)
    moments = [1 => Moment(; s=1, g=1),2 => Moment(; s=5/2, g=1)]
    sys = System(cryst, moments, :dipole)
    J1 = 1
    set_exchange!(sys, J1, Bond(1, 2, [0, 0, 0]))
    measure = ssf_trace(sys; apply_g=false)
    scga = Sunny.SCGA(sys; measure, sublattice_resolved=true, Nq=40)
    qs = q_space_path(cryst,[[0, 0, 0], [2, 0,0 ]], 17)
    kT = 17.5 * meV_per_K
    res_SCGA = Sunny.intensities_static(scga, qs; kT, λs_init = [6.742957842952556, 1.078873254872408])
    golden_data = [5.290737305656931, 4.7748982402676265, 4.08328246178467, 4.042960776695628, 4.590499089153959, 6.013520972997832, 8.995831669892262, 13.95353903683679, 17.263753113653635, 13.95353903683679, 8.995831669892262, 6.013520972997832, 4.590499089153959, 4.042960776695629, 4.0832824617846715, 4.7748982402676265, 5.290737305656931]
    @test sum(abs.(golden_data - res_SCGA.data)./ golden_data)/length(golden_data) < tol
end

@testitem "SCGA Kitchen sink" begin
    using LinearAlgebra
    tol = 1e-9
    a = 5.
    c = 17.
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    positions = [[2/3, 1/3, 0.5], [1/3, 2/3, 0.25]]
    types = ["Fe", "Mn"]
    cryst = Crystal(latvecs, positions, 148; types)
    moments = [1 => Moment(; s=1, g=3.4), 7 => Moment(; s=2, g=2)]
    sys = System(cryst, moments, :dipole; seed=2)
    J1 = -2.0*diagm([1.,1.,1.2])
    J2 = 5.25*[1 0.05 0; 0.05 1 0; 0 0 0]
    J3 = 0.75*diagm([1.,0.25,1.23]) + 0.23dmvec([0.23,0.87,0.43])
    J4 = 0.25*diagm([0.34,0.12,1.16]) + 0.23dmvec([0.98,0.65,0.353])
    J5 = 0.23*[1 0 0; 0 1 0.34; 0 0.34 0]
    J6 = -0.1*diagm([1.,1.,1.8]) + 0.086dmvec([0,0,1])
    J7 = 0.098*diagm([1.,1.,1.8]) + 0.021dmvec([0.2,0.34,0.65])
    J8 =  -0.23*diagm([1.,0.8,1.]) + 0.021dmvec([0.94,0.24,0.15])
    J9 = 1*diagm([1.21,1.11,0.76]) + 0.3dmvec([0.7,0.61,0.62])
    J10 = -0.265diagm([1,1.,0.45])
    D1 = -0.27
    D2 = 0.12
    set_exchange!(sys,J1,Bond(7, 8, [0, 0, 0])) # Mn dimer XXZ
    set_exchange!(sys,J2,Bond(1, 2, [0, 0, 0])) # Fe HC 1nn XYZ + PsD
    set_exchange!(sys,J3,Bond(1, 7, [0, 0, 0])) # Mn-Fe  XYZ + DMI [GHI]
    set_exchange!(sys,J4,Bond(1, 8, [0, 0, 0])) # Mn-Fe  XYZ + DMI [GHI]
    set_exchange!(sys,J7,Bond(1, 1, [1, 0, 0])) # Fe HC 2nn XYZ + DMI [G H I]
    set_exchange!(sys,J8,Bond(7, 7, [1, 0, 0])) # Mn in-plane XYZ + DMI [G H I]
    set_onsite_coupling!(sys, S -> D1*S[3]^2, 1)
    set_onsite_coupling!(sys, S -> D2*S[3]^2, 7)
    kT = 80.5*meV_per_K
    measure = ssf_perp(sys;)
    scga = Sunny.SCGA(sys; measure, sublattice_resolved=true, Nq=10)
    qs =  [[0, 0, 0], [0, 1/2, 1/2], [0.06, 0.49, 0.59]]
    λs_init = vcat(fill(22.47960976417936, 6), fill(5.706769019158553, 6))
    res = Sunny.intensities_static(scga, qs; kT, λs_init)
    golden_data = [37.87398338375552, 23.916275650938417, 23.911706638987532]
    @test isapprox(res.data, golden_data; atol=1e-5)
end
