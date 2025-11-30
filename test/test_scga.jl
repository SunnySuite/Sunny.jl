@testitem "Square lattice" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    s = 1
    sys = System(cryst, [1 => Moment(; s, g=1)], :dipole; seed=0)
    set_exchange!(sys, -1, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, 0.5, Bond(1, 1, [1, 1, 0]))
    set_exchange!(sys, 0.25, Bond(1, 1, [2, 0, 0]))
    measure = ssf_perp(sys)
    kT = 27.5*meV_per_K
    scga = SCGA(sys; measure, kT, dq=1/8)
    path = q_space_path(cryst, [[-1, -1, 0], [-0.04, -1, 0]], 9)
    res = Sunny.intensities_static(scga, path)
    ref_jscga = [0.3578208506624103, 0.3850668587212228, 0.4211633125595451, 0.3930558742921638, 0.3586715762816389,
                 0.3775147329089877, 0.4188431376589759, 0.4009835864744404, 0.36119385315852837]
    ref_jscga *= s^2 * Sunny.natoms(cryst)
    @test isapprox(vec(res.data)/2, ref_jscga; rtol=1e-5)
end

@testitem "Diamond lattice" begin
    a = 8.5031 # (Å)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]], 227; choice="1")
    s = 3/2
    sys = System(cryst, [1 => Moment(; s, g=1)], :dipole; seed=0)
    sys = reshape_supercell(sys, primitive_cell(cryst))
    set_exchange!(sys, -1/s^2, Bond(1, 3, [0, 0, 0]))
    set_exchange!(sys, 1/(2s)^2, Bond(1, 2, [0, 0, 0]))
    measure = ssf_perp(sys)
    kT = 15*meV_per_K
    scga = SCGA(sys; measure, kT, dq=1/8)
    grid = q_space_grid(cryst, [1, 0, 0], range(0, 3.6, 5), [0, 1, 0], range(0, 0.9, 2))
    res = intensities_static(scga, grid)
    ref_jscga = [0.7046602277469309, 0.8230846832863896, 0.23309034250417973, 0.40975668535137943, 0.8474163786642979,
                 0.8230846832694241, 0.723491683211756, 0.5939752161027589, 0.6506966347286152, 0.8012263819500781]
    ref_jscga *= s^2 * Sunny.natoms(cryst)
    @test isapprox(vec(res.data), ref_jscga; rtol=1e-8)
end

@testitem "MgCr2O4" begin
    using LinearAlgebra
    # Reproduce calculation in Bai et al., Phys. Rev. Lett. 122, 097201 (2019).
    # See also Conlon and Chalker, Phys. Rev. B 81, 224413 (2010).
    latvecs = lattice_vectors(8.3342, 8.3342, 8.3342, 90, 90, 90)
    positions = [[1/2, 1/2, 1/2]]
    cryst = Crystal(latvecs, positions, 227)
    s = 3/2
    sys = System(cryst, [1 => Moment(; s, g=1)], :dipole)
    sys = reshape_supercell(sys, primitive_cell(cryst))
    set_spin_rescaling_for_static_sum_rule!(sys)
    J1 = 3.27                                              # In meV, from Bai's PRL
    J_mgcro = [1.0, 0.0815, 0.1050, 0.0085] * J1           # Further exchanges for MgCr2O4
    set_exchange!(sys, J_mgcro[1], Bond(1, 2, [0, 0, 0]))  # J1
    set_exchange!(sys, J_mgcro[2], Bond(1, 7, [0, 0, 0]))  # J2
    set_exchange!(sys, J_mgcro[3], Bond(1, 3, [1, 0, 0]))  # J3a
    set_exchange!(sys, J_mgcro[4], Bond(1, 3, [0, 0, 0]))  # J3b
    measure = ssf_custom((q, ssf) -> real(ssf[1,1]), sys)
    kT = 20*meV_per_K
    scga = SCGA(sys; measure, kT, dq=1/4)
    grid = q_space_grid(cryst, [1, 0, 0], range(-1.5, 1.5, 4), [0, 1, 0], range(-1.5, 1.5, 4))
    res = Sunny.intensities_static(scga, grid)
    ref = [2.4168819 1.237733 1.237733 2.4168819;
           1.237733 0.16242284 0.16242284 1.237733;
           1.237733 0.16242284 0.16242284 1.237733;
           2.4168819 1.237733 1.237733 2.4168819]
    ref *= s*(s+1) / det(primitive_cell(cryst)) # Differing normalization convention
    @test isapprox(vec(res.data), vec(ref); rtol=1e-3)
end

@testitem "Arbitrary Anisotropy" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0 ,0]], 1)
    sys = System(cryst, [1 => Moment(; s=1, g=1)], :dipole; seed=0)
    set_exchange!(sys, -1, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, 0.5, Bond(1, 1, [1, 1, 0]))
    set_exchange!(sys, 0.5, Bond(1, 1, [-1, 1, 0]))
    set_exchange!(sys, 0.25, Bond(1, 1, [2, 0, 0]))
    set_exchange!(sys, 0.25, Bond(1, 1, [0, 2, 0]))
    to_inhomogeneous(sys)
    anis = [0.23 0.56 0.34;
            0.56 0.12 0.45;
            0.34 0.45 0.67]
    anis = 0.5.*(anis + anis')
    set_onsite_coupling!(sys, S -> S'*anis*S, 1)
    measure = ssf_perp(sys)
    kT = 55*meV_per_K
    scga = SCGA(sys; measure, kT, dq=1/8)
    path = q_space_path(cryst, [[0.125, 0.625, 0], [1, 0.625, 0]], 8)
    res = Sunny.intensities_static(scga, path)
    ref_jscga = [0.700944539713289, 0.6175127864675868, 0.5205893697530085, 0.48172879530096047,
                 0.5219040226511135, 0.6218544838522482, 0.7110417061581527, 0.738269833048121]
    @test isapprox(vec(res.data), ref_jscga; rtol=1e-6)
end

@testitem "Ferrimagnetic chain sum rule" begin
    latvecs = lattice_vectors(3, 5, 8, 90, 90, 90)
    positions = [[0, 0, 0], [0.5, 0, 0]]
    types = ["Ni2", "Fe3"]
    cryst = Crystal(latvecs, positions;types)
    s1 = 1
    s2 = 5/2
    moments = [1 => Moment(; s=s1, g=1), 2 => Moment(; s=s2, g=1)]
    sys = System(cryst, moments, :dipole)
    J1 = 1
    set_exchange!(sys, J1, Bond(1, 2, [0, 0, 0]))
    measure = ssf_trace(sys; apply_g=false)
    kT = 22.5*meV_per_K
    scga = SCGA(sys; measure, kT, dq=1/100)
    qs = [[qx, 0, 0] for qx in range(0, 2, 100)]
    res = Sunny.intensities_static(scga, qs)
    @test isapprox(sum(res.data)/length(qs), s1^2 + s2^2; rtol=1e-2)
end

@testitem "Ferrimagnetic chain" begin
    latvecs = lattice_vectors(3, 5, 8, 90, 90, 90)
    positions = [[0, 0, 0], [0.5, 0, 0]]
    types = ["Ni2", "Fe3"]
    cryst = Crystal(latvecs, positions; types)
    moments = [1 => Moment(; s=1, g=1), 2 => Moment(; s=5/2, g=1)]
    sys = System(cryst, moments, :dipole)
    J1 = 1
    set_exchange!(sys, J1, Bond(1, 2, [0, 0, 0]))
    measure = ssf_trace(sys; apply_g=false)
    kT = 17.5 * meV_per_K
    scga = SCGA(sys; measure, kT, dq=1/40)
    qs = q_space_path(cryst, [[0, 0, 0], [2, 0, 0]], 5)
    res = Sunny.intensities_static(scga, qs)
    # println(round.(res.data; digits=10))
    golden_data = [5.2798224086, 4.8642933106, 16.3317444291, 4.8642933106, 5.2798224086]
    @test isapprox(golden_data, res.data; rtol=1e-7)
end

@testitem "SCGA Kitchen sink" begin
    using LinearAlgebra
    a = 5.
    c = 17.
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    positions = [[2/3, 1/3, 0.5], [1/3, 2/3, 0.25]]
    types = ["Fe", "Mn"]
    cryst = Crystal(latvecs, positions, 148; types)
    moments = [1 => Moment(; s=1, g=3.4), 7 => Moment(; s=2, g=2)]
    sys = System(cryst, moments, :dipole)
    sys = reshape_supercell(sys, primitive_cell(cryst)) # Works either way!
    J1 = -2.0*diagm([1, 1, 1.2])
    J2 = 5.25*[1 0.05 0; 0.05 1 0; 0 0 0]
    J3 = 0.75*diagm([1, 0.25, 1.23]) + 0.23dmvec([0.23, 0.87, 0.43])
    J4 = 0.25*diagm([0.34, 0.12, 1.16]) + 0.23dmvec([0.98, 0.65, 0.353])
    J5 = 0.23*[1 0 0; 0 1 0.34; 0 0.34 0]
    J6 = -0.1*diagm([1, 1, 1.8]) + 0.086dmvec([0, 0, 1])
    J7 = 0.098*diagm([1, 1, 1.8]) + 0.021dmvec([0.2, 0.34, 0.65])
    J8 =  -0.23*diagm([1, 0.8, 1]) + 0.021dmvec([0.94, 0.24, 0.15])
    J9 = 1*diagm([1.21, 1.11, 0.76]) + 0.3dmvec([0.7, 0.61, 0.62])
    J10 = -0.265diagm([1, 1, 0.45])
    D1 = -0.27
    D2 = 0.12
    set_exchange!(sys, J1, Bond(7, 8, [0, 0, 0])) # Mn dimer XXZ
    set_exchange!(sys, J2, Bond(1, 2, [0, 0, 0])) # Fe HC 1nn XYZ + PsD
    set_exchange!(sys, J3, Bond(1, 7, [0, 0, 0])) # Mn-Fe  XYZ + DMI [GHI]
    set_exchange!(sys, J4, Bond(1, 8, [0, 0, 0])) # Mn-Fe  XYZ + DMI [GHI]
    set_exchange!(sys, J7, Bond(1, 1, [1, 0, 0])) # Fe HC 2nn XYZ + DMI [G H I]
    set_exchange!(sys, J8, Bond(7, 7, [1, 0, 0])) # Mn in-plane XYZ + DMI [G H I]
    set_onsite_coupling!(sys, S -> D1*S[3]^2, 1)
    set_onsite_coupling!(sys, S -> D2*S[3]^2, 7)
    measure = ssf_perp(sys)
    kT = 80.5*meV_per_K
    scga = SCGA(sys; measure, kT, dq=1/4)
    qs = [[0, 0, 0], [0, 1/2, 1/2], [0.06, 0.49, 0.59]]
    res = Sunny.intensities_static(scga, qs)
    # println(round.(res.data; digits=10))
    golden_data = [129.5389081127, 111.0369328807, 112.1694261822]
    @test isapprox(res.data, golden_data; rtol=1e-8)
end
