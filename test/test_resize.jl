@testitem "Stripe order" begin
    using LinearAlgebra, IOCapture

    latvecs = lattice_vectors(1, 1, 1.1, 90, 90, 120)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, (4, 6, 2), [SpinInfo(1, S=1)], :dipole)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, 1.2, Bond(1, 1, [0, 0, 1]))
    set_exchange!(sys, diagm([0.2, 0.3, 0.4]), Bond(1, 1, [1, 0, 1]))

    for idx in Sunny.all_sites(sys)
        row = idx[2]
        dir = [1, 0, 2mod(row, 2) - 1]
        Sunny.setspin!(sys, Sunny.dipolarspin(sys, idx, dir), idx)
    end

    capt = IOCapture.capture() do
        print_wrapped_intensities(sys)
    end
    @test capt.output == """
        Dominant wavevectors for spin sublattices:

            [0, 0, 0]              50.00% weight
            [0, 1/2, 0]            50.00%
        """

    capt = IOCapture.capture() do
        suggest_magnetic_supercell([[0,1/2,0]], sys.latsize)
    end
    @test capt.output == """
        Suggested magnetic supercell in multiples of lattice vectors:

            [1 0 0; 0 2 0; 0 0 1]

        for wavevectors [[0, 1/2, 0]].
        """

    A1 = [1 0 0; 0 2 0; 0 0 1]
    A2 = [1 0 0; 1 2 0; 0 0 1]
    newsys1 = reshape_geometry(sys, A1)
    newsys2 = reshape_geometry(sys, A2)

    @test energy(sys) / prod(sys.latsize) ≈ 2.55
    
    newsys = reshape_geometry(sys, A1)
    @test energy(newsys) / prod(newsys.latsize) ≈ 2.55
    newsys = reshape_geometry(sys, A2)
    @test energy(newsys) / prod(newsys.latsize) ≈ 2.55
    newsys = reshape_geometry(sys, A1)
    @test energy(newsys) / prod(newsys.latsize) ≈ 2.55
end

@testitem "Equivalent reshaping" begin
    using LinearAlgebra

    latvecs = lattice_vectors(1, 1, 1, 90, 90, 120)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, (3, 3, 3), [SpinInfo(1, S=1)], :dipole)
    randomize_spins!(sys)

    # Reshape to sheared volume
    sys2 = reshape_geometry(sys, [3 0 0; 2 3 0; 0 0 3])
    # Reshape back to original volume
    sys3 = reshape_geometry(sys2, diagm([3,3,3]))
    @test sys2.dipoles != sys.dipoles
    @test sys3.dipoles == sys.dipoles

    # Two equivalent ways of sizing up sys
    sys2 = repeat_periodically(sys, (2, 2, 1))
    sys3 = resize_periodically(sys, (6, 6, 3))
    @test sys2.dipoles == sys3.dipoles
end


@testitem "Interactions after reshaping" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, (3, 3, 3), [SpinInfo(1, S=1)], :dipole)
    randomize_spins!(sys)
    
    # Commensurate shear that is specially designed to preserve the periodicity of
    # the system volume
    sys2 = reshape_geometry(sys, [3 3 0; 0 3 0; 0 0 3])
    
    # Users always specify a bond using atom indices of the original unit cell,
    # but `sys2.interactions_union` is internally reindexed.
    set_exchange!(sys,  1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys2, 1.0, Bond(1, 1, [1, 0, 0]))
    
    @test energy(sys) ≈ energy(sys2)
end


@testitem "FeI2 equivalent energy" begin
    using LinearAlgebra

    a = b = 4.05012
    c = 6.75214
    lat_vecs = lattice_vectors(a, b, c, 90, 90, 120)
    cryst = Crystal(lat_vecs, [[0,0,0]], 164)
    sys = System(cryst, (4,4,4), [SpinInfo(1,S=1)], :SUN, seed=0)

    J1pm   = -0.236 
    J1pmpm = -0.161
    J1zpm  = -0.261
    J2pm   = 0.026
    J3pm   = 0.166
    J′0pm  = 0.037
    J′1pm  = 0.013
    J′2apm = 0.068

    J1zz   = -0.236
    J2zz   = 0.113
    J3zz   = 0.211
    J′0zz  = -0.036
    J′1zz  = 0.051
    J′2azz = 0.073

    J1xx = J1pm + J1pmpm 
    J1yy = J1pm - J1pmpm
    J1yz = J1zpm

    set_exchange!(sys, [J1xx   0.0    0.0;
                        0.0    J1yy   J1yz;
                        0.0    J1yz   J1zz], Bond(1,1,[1,0,0]))
    set_exchange!(sys, diagm([J2pm, J2pm, J2zz]), Bond(1,1,[1,2,0]))
    set_exchange!(sys, diagm([J3pm, J3pm, J3zz]), Bond(1,1,[2,0,0]))
    set_exchange!(sys, diagm([J′0pm, J′0pm, J′0zz]), Bond(1,1,[0,0,1]))
    set_exchange!(sys, diagm([J′1pm, J′1pm, J′1zz]), Bond(1,1,[1,0,1]))
    set_exchange!(sys, diagm([J′2apm, J′2apm, J′2azz]), Bond(1,1,[1,2,1]))

    D = 2.165
    set_anisotropy!(sys, -D*𝒮[3]^2, 1)

    # periodic ground state for FeI2
    s = Sunny.Vec3(1, 0, 0)
    for j = 1:4, k = 1:4
        sys.dipoles[:, j, k, 1] = circshift([s, s, -s, -s], mod(j+k, 4))
    end

    sys_supercell = reshape_geometry(sys, [2 0 1; -1 1 0; -1 -1 1])
    sys2 = resize_periodically(sys_supercell, (4,4,4))

    E0 = energy(sys) / length(sys.dipoles)
    E1 = energy(sys_supercell) / length(sys_supercell.dipoles)
    E2 = energy(sys2) / length(sys2.dipoles)
    @test E0 ≈ E1 ≈ E2
end
