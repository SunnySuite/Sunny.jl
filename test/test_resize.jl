@testitem "Stripe order" begin
    using LinearAlgebra, IOCapture

    latvecs = lattice_vectors(1, 1, 1.1, 90, 90, 120)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(4, 6, 2))
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, 1.2, Bond(1, 1, [0, 0, 1]))
    set_exchange!(sys, diagm([0.2, 0.3, 0.4]), Bond(1, 1, [1, 0, 1]))

    for site in eachsite(sys)
        row = site[2]
        dir = [1, 0, 2mod(row, 2) - 1]
        set_dipole!(sys, dir, site)
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
        suggest_magnetic_supercell([[0,1/2,0]])
    end
    @test capt.output == """
        Possible magnetic supercell in multiples of lattice vectors:

            [1 0 0; 0 2 0; 0 0 1]

        for the rationalized wavevectors:

            [[0, 1/2, 0]]
        """
    capt = IOCapture.capture() do
        suggest_magnetic_supercell([[1/2,0,0], [0,1/2,0]])
    end
    @test capt.output == """
        Possible magnetic supercell in multiples of lattice vectors:

            [2 0 0; 0 2 0; 0 0 1]

        for the rationalized wavevectors:

            [[1/2, 0, 0], [0, 1/2, 0]]
        """

    A1 = [1 0 0; 0 2 0; 0 0 1]
    A2 = [1 0 0; 1 2 0; 0 0 1]
    newsys1 = reshape_supercell(sys, A1)
    newsys2 = reshape_supercell(sys, A2)

    @test energy_per_site(sys) ≈ 2.55

    newsys = reshape_supercell(sys, A1)
    @test energy_per_site(newsys) ≈ 2.55
    newsys = reshape_supercell(sys, A2)
    @test energy_per_site(newsys) ≈ 2.55
    newsys = reshape_supercell(sys, A1)
    @test energy_per_site(newsys) ≈ 2.55
end


@testitem "Validate reshaping" begin
    using LinearAlgebra
    
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    positions = [[0, 0, 0]/2, [1, 1, 1]/2]
    cryst = Crystal(latvecs, positions)
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
    
    primcell = primitive_cell(cryst)
    reshape_supercell(sys, Diagonal([2, 3, 4])) # Fine
    reshape_supercell(sys, primcell) # Fine
    reshape_supercell(sys, primcell * Diagonal([2, 3, 4])) # Fine
    msg = "Elements of `primitive_cell(cryst) \\ shape` must be integer. Calculated [3.5 0.5 -0.5; 1.0 3.0 -1.0; 0.5 -0.5 2.5]."
    @test_throws msg reshape_supercell(sys, Diagonal([2, 3, 4]) * primcell)
    
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
    reshape_supercell(sys, Diagonal([2.0, 3, 4]))
    msg = "Elements of `shape` must be integer. Received [2.5 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0]."
    @test_throws msg reshape_supercell(sys, Diagonal([2.5, 3, 4]))
end


@testitem "Equivalent reshaping" begin
    using LinearAlgebra

    latvecs = lattice_vectors(1, 1, 1, 90, 90, 120)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(3, 3, 3))
    randomize_spins!(sys)

    # Reshape to sheared volume
    sys2 = reshape_supercell(sys, [3 0 0; 2 3 0; 0 0 3])
    # Reshape back to original volume
    sys3 = reshape_supercell(sys2, diagm([3,3,3]))
    @test sys2.dipoles != sys.dipoles
    @test sys3.dipoles == sys.dipoles

    # Two equivalent ways of sizing up sys
    sys2 = repeat_periodically(sys, (2, 2, 1))
    sys3 = resize_supercell(sys, (6, 6, 3))
    @test sys2.dipoles == sys3.dipoles
end


@testitem "Interactions after reshaping" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(3, 3, 3))
    randomize_spins!(sys)
    
    # Commensurate shear that is specially designed to preserve the periodicity of
    # the system volume
    sys2 = reshape_supercell(sys, [3 3 0; 0 3 0; 0 0 3])
    
    # Users always specify a bond using atom indices of the original unit cell,
    # but `sys2.interactions_union` is internally reindexed.
    set_exchange!(sys,  1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys2, 1.0, Bond(1, 1, [1, 0, 0]))
    
    @test energy(sys) ≈ energy(sys2)
end


@testitem "Reshape then inhomogeneous" begin
    using LinearAlgebra

    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst = Crystal(latvecs, [[0,0,0]])
    dims = (6, 6, 2)
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims)
    polarize_spins!(sys, (0,0,1))

    J = 1.5
    bond = Bond(1, 1, [1, 0, 0])
    set_exchange!(sys, J, bond)

    E0 = J * prod(dims) * (6/2)
    @test energy(sys) == E0

    A = [dims[1] dims[2]÷2 0
         0       dims[2]   0
         0       0         dims[3]]
    @test det(A) ≈ prod(dims)
    sys2 = reshape_supercell(sys, A)
    @test Sunny.natoms(sys2.crystal) == 2
    @test energy(sys2) ≈ E0

    for sys′ in to_inhomogeneous.((sys, sys2))
        @test energy(sys′) ≈ E0

        # Check that we can overwrite all the bonds in inhomogeneous mode
        for (site1, site2, offset) in symmetry_equivalent_bonds(sys′, bond)
            set_exchange_at!(sys′, 2J, site1, site2; offset)
        end
        @test energy(sys′) ≈ 2*E0
    end
end


@testitem "Remove periodicity" begin
    a = 1
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])

    bond = Bond(1, 1, (1, 0, 0))
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(3, 3, 3), seed=0)
    set_exchange!(sys, 1, bond)
    @test energy(sys) == 81

    # Remove periodic boundaries in certain directions
    sys2 = to_inhomogeneous(sys)
    remove_periodicity!(sys2, (true, false, false))
    @test energy(sys2) == 72
    remove_periodicity!(sys2, (true, false, true))
    @test energy(sys2) == 63
end


@testitem "FeI2 equivalent energy" begin
    using LinearAlgebra

    a = b = 4.05012
    c = 6.75214
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    cryst = Crystal(latvecs, [[0,0,0]], 164)
    sys = System(cryst, [1 => Moment(s=1, g=2)], :SUN; dims=(4, 4, 4), seed=0)

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
    set_onsite_coupling!(sys, S -> -D*S[3]^2, 1)

    # periodic ground state for FeI2
    s = Sunny.Vec3(1, 0, 0)
    for j = 1:4, k = 1:4
        sys.dipoles[:, j, k, 1] = circshift([s, s, -s, -s], mod(j+k, 4))
    end

    sys_supercell = reshape_supercell(sys, [2 0 1; -1 1 0; -1 -1 1])
    sys2 = resize_supercell(sys_supercell, (4,4,4))

    E0 = energy_per_site(sys)
    E1 = energy_per_site(sys_supercell)
    E2 = energy_per_site(sys2)
    @test E0 ≈ E1 ≈ E2
end
