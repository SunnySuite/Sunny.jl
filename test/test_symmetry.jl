@testitem "Crystal Construction" begin
    include("shared.jl")

    cell_type(cryst::Crystal) = Sunny.cell_type(cryst.latvecs)
    lattice_params(cryst::Crystal) = Sunny.lattice_params(cryst.latvecs)

    ### Test construction of diamond lattice

    # Spglib inferred symmetry
    latvecs = [1 1 0; 0 1 1; 1 0 1]' / 2
    positions = [[1, 1, 1], [-1, -1, -1]] / 8
    cryst = Crystal(latvecs, positions)
    ref_bonds = reference_bonds(cryst, 2.)
    dist1 = [Sunny.global_distance(cryst, b) for b in ref_bonds]

    # Using explicit symops
    latvecs = Sunny.Mat3(latvecs)
    positions = [Sunny.Vec3(1, 1, 1) / 8]
    types = [""]
    cryst = Sunny.crystal_from_symops(latvecs, positions, types, cryst.symops, cryst.spacegroup)
    ref_bonds = reference_bonds(cryst, 2.)
    dist2 = [Sunny.global_distance(cryst, b) for b in ref_bonds]

    # Using Hall number
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90) # must switch to standard cubic unit cell
    positions = [Sunny.Vec3(1, 1, 1) / 4]
    cryst = Sunny.crystal_from_hall_number(latvecs, positions, types, 525)
    ref_bonds = reference_bonds(cryst, 2.)
    dist3 = [Sunny.global_distance(cryst, b) for b in ref_bonds]

    # Using international symbol
    positions = [[1, 1, 1] / 4]
    # cryst = Crystal(latvecs, positions, "F d -3 m") # Ambiguous!
    cryst = Crystal(latvecs, positions, "F d -3 m"; setting="1")
    ref_bonds = reference_bonds(cryst, 2.)
    dist4 = [Sunny.global_distance(cryst, b) for b in ref_bonds]

    @test dist1 ≈ dist2 ≈ dist3 ≈ dist4



    ### FCC lattice, primitive vs. standard unit cell

    latvecs = [1 1 0; 0 1 1; 1 0 1]' / 2
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)

    latvecs = [1 0 0; 0 1 0; 0 0 1]'
    positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    cryst′ = Crystal(latvecs, positions)

    @test cryst.sitesyms[1] == cryst′.sitesyms[1]

    # Calculate interaction table
    ref_bonds = reference_bonds(cryst, 2.)
    b = ref_bonds[2]
    basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, b)
    J = basis' * randn(length(basis))
    (bs, Js) = Sunny.all_symmetry_related_couplings_for_atom(cryst, b.i, b, J)
    @test length(Js) == Sunny.coordination_number(cryst, b.i, b)


    ### Triangular lattice, primitive unit cell

    c = 10
    latvecs = [1 0 0;  -1/2 √3/2 0;  0 0 c]'
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    @test cell_type(cryst) == Sunny.hexagonal
    @test Sunny.natoms(cryst) == 1
    @test Sunny.cell_volume(cryst) ≈ c * √3 / 2 
    @test all(lattice_params(cryst) .≈ (1., 1., c, 90., 90., 120.))

    ### Kagome lattice

    latvecs = [1 0 0;  -1/2 √3/2 0;  0 0 c]'
    positions = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
    cryst = Crystal(latvecs, positions)
    @test cell_type(cryst) == Sunny.hexagonal
    @test Sunny.natoms(cryst) == 3
    @test Sunny.cell_volume(cryst) ≈ c * √3 / 2 
    @test all(lattice_params(cryst) .≈ (1., 1., c, 90., 90., 120.))


    ### Arbitrary monoclinic

    mono_lat_params = (6, 7, 8, 90, 90, 40)
    latvecs = lattice_vectors(mono_lat_params...)
    positions = [[0,0,0]]
    # cryst = Crystal(latvecs, positions, "C 2/c")
    cryst = Crystal(latvecs, positions, "C 2/c", setting="c1")
    @test cell_type(cryst) == Sunny.monoclinic
    @test Sunny.natoms(cryst) == 4
    @test all(lattice_params(cryst) .≈ mono_lat_params)


    ### Arbitrary trigonal

    latvecs = lattice_vectors(5, 5, 6, 90, 90, 120)
    positions = [[0,0,0]]
    cryst1 = Crystal(latvecs, positions, "P -3")
    @test Sunny.natoms(cryst1) == 1
    @test cell_type(cryst1) == Sunny.hexagonal
    cryst2 = Crystal(latvecs, positions, "R -3")
    @test Sunny.natoms(cryst2) == 3
    cryst3 = Crystal(latvecs, positions, 147) # spacegroup number
    @test cell_type(cryst1) == cell_type(cryst2) == cell_type(cryst3) == Sunny.hexagonal


    ### Arbitrary triclinic

    latvecs = lattice_vectors(6, 7, 8, 70, 80, 90)
    positions = [[0,0,0]]
    cryst1 = Crystal(latvecs, positions, "P 1")
    @test Sunny.natoms(cryst1) == 1
    cryst2 = Crystal(latvecs, positions) # Infers 'P -1'
    @test Sunny.natoms(cryst1) == Sunny.natoms(cryst2) == 1
    @test cell_type(cryst1) == cell_type(cryst2) == Sunny.triclinic

    ### Orthorhombic test, found by Ovi Garlea

    latvecs = lattice_vectors(13.261, 7.718, 6.278, 90.0, 90.0, 90.0);
    types = ["Yb1","Yb2"];
    positions = [[0,0,0], [0.266,0.25,0.02]]; # Locations of atoms as multiples of lattice vectors
    capt = IOCapture.capture() do
        Crystal(latvecs, positions, 62; types, symprec=1e-4)
    end
    @test capt.output == """
        The spacegroup '62' allows for multiple settings!
        Returning a list of the possible crystals:
            1. "P n m a", setting="", with  8 atoms
            2. "P m n b", setting="ba-c", with 12 atoms
            3. "P b n m", setting="cab", with 12 atoms
            4. "P c m n", setting="-cba", with  8 atoms
            5. "P m c n", setting="bca", with 12 atoms
            6. "P n a m", setting="a-cb", with 12 atoms
        
        Note: To disambiguate, you may pass a named parameter, setting="...".
        
        """
    @test length(capt.value) == 6
    cryst = Crystal(latvecs, positions, 62; types, symprec=1e-4, setting="-cba")
    @test count(==(1), cryst.classes) == 4
    @test count(==(2), cryst.classes) == 4    
end


@testitem "Allowed anisotropy" begin
    using LinearAlgebra
    
    # Test some inferred anisotropy matrices
    let
        N = 7
        k = 6
        i = 1
        cryst = Sunny.diamond_crystal()
        O = Sunny.StevensMatrices(N)

        # print_site(cryst, i)
        Λ = O[6,0]-21O[6,4]
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)

        R = [normalize([1 1 -2]); normalize([-1 1 0]); normalize([1 1 1])]
        # print_site(cryst, i; R)
        Λ = O[6,0]-(35/√8)*O[6,3]+(77/8)*O[6,6]
        Λ′ = rotate_operator(Λ, R)
        @test Sunny.is_anisotropy_valid(cryst, i, Λ′)

        latvecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)
        cryst = Crystal(latvecs, [[0., 0., 0.]])
        # print_site(cryst, i)
        Λ = randn()*(O[6,0]-21O[6,4]) + randn()*(O[6,2]+(16/5)*O[6,4]+(11/5)*O[6,6])
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)
    end

    # Test validity of symmetry inferred anisotropies 
    let
        latvecs = [1 0 0; 0 1 0; 0 0 10]'
        # All atoms are a distance of 0.1 from the origin, and are arranged at
        # angles 0, π/2, π, 3π/2, counterclockwise.
        positions = [[0.1, 0, 0], [0, 0.1, 0], [-0.1, 0, 0], [0, -0.1, 0]]
        cryst = Crystal(latvecs, positions)

        for mode in (:dipole, :SUN)
            sys = System(cryst, (1,1,1), [SpinInfo(1, S=2)], mode)
            randomize_spins!(sys)

            # Most general allowed anisotropy for this crystal
            c = randn(9)
            O = stevens_operators(sys, 1)
            Λ = sum(c .* [O[2,0], O[2,2], O[4,0], O[4,2], O[4,4], O[6,0], O[6,2], O[6,4], O[6,6]])
            set_onsite!(sys, Λ, 1)

            E1 = energy(sys)

            # By circularly shifting the atom index, we are effectively rotating
            # site positions by π/2 clockwise (see comment above `positions`).
            # Rotate spin state correspondingly
            R = Sunny.Mat3([0 1 0; -1 0 0; 0 0 1])
            sys.dipoles .= circshift(sys.dipoles, (0,0,0,1))
            sys.dipoles .= [R*d for d in sys.dipoles]

            # If coherents are present, perform same operation
            if mode == :SUN
                U = Sunny.unitary_for_rotation(R; N=sys.Ns[1])
                sys.coherents .= circshift(sys.coherents, (0,0,0,1))
                sys.coherents .= [U*z for z in sys.coherents]
            end

            # Verify energy is invariant
            E2 = energy(sys)
            @test E1 ≈ E2
        end
    end
end


@testitem "Renormalization" begin
    latvecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)
    cryst = Crystal(latvecs, [[0., 0., 0.]])
    
    
    # Dipole system with renormalized anisotropy
    sys0 = System(cryst, (1,1,1), [SpinInfo(1, S=3)], :dipole)
    randomize_spins!(sys0)

    i = 1
    O = stevens_operators(sys0, i)
    Λ = randn()*(O[2,0]+3O[2,2]) +
        randn()*(O[4,0]-5O[4,2]) + randn()*(O[4,0]+5O[4,4]) +
        randn()*(O[6,0]-21O[6,4]) + randn()*(O[6,0]+(105/16)O[6,2]+(231/16)O[6,6])
    set_onsite!(sys0, Λ, i)
    E0 = energy(sys0)
    
    # Corresponding SU(N) system
    sys = System(cryst, (1,1,1), [SpinInfo(1, S=3)], :SUN)
    for site in all_sites(sys)
        polarize_spin!(sys, sys0.dipoles[site], site)
    end
    set_onsite!(sys, Λ, i)
    E = energy(sys)
    
    @test E ≈ E0    
end
