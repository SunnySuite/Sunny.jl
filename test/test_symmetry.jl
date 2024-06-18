@testitem "Crystal Construction" begin
    using IOCapture

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
    bond = ref_bonds[2]
    basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, bond)
    bs = Sunny.all_symmetry_related_bonds_for_atom(cryst, bond.i, bond)
    @test length(bs) == Sunny.coordination_number(cryst, bond.i, bond)

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


@testitem "Standardize" begin
    using LinearAlgebra

    function test_standardize(latvecs, r)
        cryst = Crystal(latvecs, [r])
        cryst2 = standardize(cryst; idealize=false)
        @test cryst2.latvecs * cryst2.positions[1] ≈ cryst.latvecs * r
        cryst3 = standardize(cryst)
        @test norm(cryst3.positions[1]) < 1e-12    
    end
    
    r = [0.1, 0.2, 0.3]
    test_standardize([1 0 1; 1 1 0; 0 1 1], r)
    test_standardize(lattice_vectors(1, 1, 1, 90, 90, 60), r)
end


@testitem "Allowed anisotropy" begin
    using LinearAlgebra
    
    # Test some inferred anisotropy matrices
    let
        # This test should also work for S = Inf, but there is a false negative
        # within second call to `is_anisotropy_valid`. TODO: Create minimized
        # test for `isapprox` bug and report to DynamicPolynomials repo. Cf.
        # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/created_by/kbarros

        S = 3
        k = 6
        i = 1
        cryst = Sunny.diamond_crystal()
        O = stevens_matrices(S)

        # print_site(cryst, i)
        Λ = O[6,0]-21O[6,4]
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)

        R = [normalize([1 1 -2]); normalize([-1 1 0]); normalize([1 1 1])]
        # print_site(cryst, i; R)
        Λ = O[6,0]-(35/√8)*O[6,3]+(77/8)*O[6,6]
        Λ′ = rotate_operator(Λ, R)
        @test Sunny.is_anisotropy_valid(cryst, i, Λ′)

        latvecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)

        warnstr = "Found a nonconventional tetragonal unit cell. Consider using `lattice_vectors(a, a, c, 90, 90, 90)`"
        cryst = @test_warn warnstr Crystal(latvecs, [[0, 0, 0]])
    
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
            sys = System(cryst, (1,1,1), [SpinInfo(1, S=2, g=2)], mode)
            randomize_spins!(sys)

            # Most general allowed anisotropy for this crystal
            c = randn(9)
            O = stevens_matrices(spin_label(sys, 1))
            Λ = sum(c .* [O[2,0], O[2,2], O[4,0], O[4,2], O[4,4], O[6,0], O[6,2], O[6,4], O[6,6]])
            set_onsite_coupling!(sys, Λ, 1)

            E1 = energy(sys)

            # By circularly shifting the atom index, we are effectively rotating
            # site positions by π/2 clockwise (see comment above `positions`).
            # Rotate spin state correspondingly
            R = Sunny.Mat3([0 1 0; -1 0 0; 0 0 1])
            sys.dipoles .= circshift(sys.dipoles, (0,0,0,1))
            sys.dipoles .= [R*d for d in sys.dipoles]

            # If coherents are present, perform same operation
            if mode == :SUN
                U = Sunny.unitary_irrep_for_rotation(R; N=sys.Ns[1])
                sys.coherents .= circshift(sys.coherents, (0,0,0,1))
                sys.coherents .= [U*z for z in sys.coherents]
            end

            # Verify energy is invariant
            E2 = energy(sys)
            @test E1 ≈ E2
        end
    end
end

@testitem "Symmetry table" begin
    using LinearAlgebra
    import IOCapture

    # Pyrochlore
    cryst = Crystal(Sunny.Mat3(I), [[0, 0, 0]], 227, setting="2")
    @test cryst.positions ≈ [
        [0, 0, 0], [1/4, 1/4, 0], [1/2, 1/2, 0], [3/4, 3/4, 0], [1/4, 0, 1/4], [0, 1/4, 1/4], [3/4, 1/2, 1/4], [1/2, 3/4, 1/4], [1/2, 0, 1/2], [3/4, 1/4, 1/2], [0, 1/2, 1/2], [1/4, 3/4, 1/2], [3/4, 0, 3/4], [1/2, 1/4, 3/4], [1/4, 1/2, 3/4], [0, 3/4, 3/4],
    ]
    capt = IOCapture.capture() do
        print_symmetry_table(cryst, 0.8)
    end
    @test capt.output == """
        Atom 1
        Position [0, 0, 0], multiplicity 16
        Allowed g-tensor: [A B B
                           B A B
                           B B A]
        Allowed anisotropy in Stevens operators:
            c₁*(𝒪[2,-2]+2𝒪[2,-1]+2𝒪[2,1]) +
            c₂*(-7𝒪[4,-3]-2𝒪[4,-2]+𝒪[4,-1]+𝒪[4,1]+7𝒪[4,3]) + c₃*(𝒪[4,0]+5𝒪[4,4]) +
            c₄*(-11𝒪[6,-6]-8𝒪[6,-3]+𝒪[6,-2]-8𝒪[6,-1]-8𝒪[6,1]+8𝒪[6,3]) + c₅*(𝒪[6,0]-21𝒪[6,4]) + c₆*((9/5)𝒪[6,-6]+(24/5)𝒪[6,-5]+𝒪[6,-2]+(8/5)𝒪[6,-1]+(8/5)𝒪[6,1]+(24/5)𝒪[6,5])
        
        Bond(1, 2, [0, 0, 0])
        Distance 0.35355339059327, coordination 6
        Connects [0, 0, 0] to [1/4, 1/4, 0]
        Allowed exchange matrix: [A C -D
                                  C A -D
                                  D D  B]
        Allowed DM vector: [-D D 0]
        
        Bond(3, 5, [0, 0, 0])
        Distance 0.61237243569579, coordination 12
        Connects [1/2, 1/2, 0] to [1/4, 0, 1/4]
        Allowed exchange matrix: [  A  C-E  D-F
                                  C+E    B -C+E
                                  D+F -C-E    A]
        Allowed DM vector: [E F -E]
        
        Bond(1, 3, [-1, 0, 0])
        Distance 0.70710678118655, coordination 6
        Connects [0, 0, 0] to [-1/2, 1/2, 0]
        Allowed exchange matrix: [A D C
                                  D A C
                                  C C B]
        
        Bond(1, 3, [0, 0, 0])
        Distance 0.70710678118655, coordination 6
        Connects [0, 0, 0] to [1/2, 1/2, 0]
        Allowed exchange matrix: [A D C
                                  D A C
                                  C C B]
        
        Bond(1, 2, [-1, 0, 0])
        Distance 0.79056941504209, coordination 12
        Connects [0, 0, 0] to [-3/4, 1/4, 0]
        Allowed exchange matrix: [A  D -F
                                  D  B  E
                                  F -E  C]
        Allowed DM vector: [E F 0]

        """

    capt = IOCapture.capture() do
        print_suggested_frame(cryst, 2)
    end
    @test capt.output == """
        R = [1/√2      0  1/√2
             1/√6 -√2/√3 -1/√6
             1/√3   1/√3 -1/√3]
        """
    
    R = [1/√2 0 1/√2; 1/√6 -√2/√3 -1/√6; 1/√3 1/√3 -1/√3]
    capt = IOCapture.capture() do
        print_site(cryst, 2; R)
    end
    @test capt.output == """
        Atom 2
        Position [1/4, 1/4, 0], multiplicity 16
        Allowed g-tensor: [A-B   0    0
                             0 A-B    0
                             0   0 A+2B]
        Allowed anisotropy in Stevens operators:
            c₁*𝒪[2,0] +
            c₂*𝒪[4,-3] + c₃*𝒪[4,0] +
            c₄*𝒪[6,-3] + c₅*𝒪[6,0] + c₆*𝒪[6,6]
        
        Modified reference frame! Transform using `rotate_operator(op, R)` where
        R = [1/√2      0  1/√2
             1/√6 -√2/√3 -1/√6
             1/√3   1/√3 -1/√3]
        """

    capt = IOCapture.capture() do
        cryst = Sunny.hyperkagome_crystal()
        print_suggested_frame(cryst, 2)
    end
    @test capt.output == """
        [ Info: Could not find a symmetry axis orthogonal to [1/√2, 1/√2, 0].
        R = [1/√2 -1/√2  0
                0     0 -1
             1/√2  1/√2  0]
        """
end


@testitem "Renormalization" begin
    latvecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)
    warnstr = "Found a nonconventional tetragonal unit cell. Consider using `lattice_vectors(a, a, c, 90, 90, 90)`"
    cryst = @test_warn warnstr Crystal(latvecs, [[0, 0, 0]])    
    
    # Dipole system with renormalized anisotropy
    sys0 = System(cryst, (1,1,1), [SpinInfo(1, S=3, g=2)], :dipole)
    randomize_spins!(sys0)

    i = 1
    O = stevens_matrices(spin_label(sys0, i))
    Λ = randn()*(O[2,0]+3O[2,2]) +
        randn()*(O[4,0]-5O[4,2]) + randn()*(O[4,0]+5O[4,4]) +
        randn()*(O[6,0]-21O[6,4]) + randn()*(O[6,0]+(105/16)O[6,2]+(231/16)O[6,6])
    set_onsite_coupling!(sys0, Λ, i)
    E0 = energy(sys0)
    
    # Corresponding SU(N) system
    sys = System(cryst, (1,1,1), [SpinInfo(1, S=3, g=2)], :SUN)
    for site in eachsite(sys)
        set_dipole!(sys, sys0.dipoles[site], site)
    end
    set_onsite_coupling!(sys, Λ, i)
    E = energy(sys)
    
    @test E ≈ E0    
end
