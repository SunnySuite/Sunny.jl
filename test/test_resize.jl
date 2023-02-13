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
        print_dominant_wavevectors(sys)
    end
    @test capt.output ==
    """
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
