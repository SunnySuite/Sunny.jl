@testitem "Wyckoff table" begin
    # Tests depending on Crystalline are disabled as this is a heavy dependency.
    #=
    import Crystalline

    for sgnum in 1:230
        wyckoffs = Crystalline.wyckoffs(sgnum, 3)
        wyckoffs = sort(wyckoffs; by = wp -> (!isuppercase(wp.letter), -Int(wp.letter)))
        for (i, wp) in enumerate(wyckoffs)
            # From Crystalline
            (; letter, mult, v) = wp
            (; cnst, free) = v

            # From Sunny table, sourced from Spglib
            (mult2, letter2, symb, pos) = Sunny.wyckoff_table[sgnum][i]
            (cnst2, free2) = Sunny.WyckoffExpr(pos)
            if sgnum == 98 && i == 3
                free2 *= -1
            end

            @test mult == mult2
            @test letter == letter2
            @test free == free2
            @test cnst == cnst2
        end
    end
    =#

    # Test that the orbit of a Wyckoff matches the multiplicity data.
    for sgnum in 1:230
        sg = Sunny.Spacegroup(Sunny.standard_setting[sgnum])
        for (mult, letter, sitesym, pos) in Sunny.wyckoff_table[sgnum]
            orbit = Sunny.crystallographic_orbit(Sunny.WyckoffExpr(pos); sg.symops)
            @test length(orbit) == mult
        end
    end

    # Test that Wyckoffs can be correctly inferred from a position
    niters = 20
    for sgnum in rand(1:230, niters)
        for (_, letter, _, pos) in Sunny.wyckoff_table[sgnum]
            w = Sunny.WyckoffExpr(pos)
            θ = 10 * randn(3)
            r = w.F * θ + w.c
            @test letter == Sunny.find_wyckoff(sgnum, r; tol=1e-12).letter
        end
    end
end

@testitem "Crystal Construction" begin
    using IOCapture

    ### Test construction of diamond lattice

    # Spglib inferred symmetry
    latvecs = [1 1 0; 0 1 1; 1 0 1]' / 2
    positions = [[1, 1, 1], [-1, -1, -1]] / 8
    cryst = Crystal(latvecs, positions)
    ref_bonds = reference_bonds(cryst, 2.)
    dist1 = [Sunny.global_distance(cryst, b) for b in ref_bonds]

    # From spacegroup symmetry
    latvecs = Sunny.Mat3(latvecs)
    positions = [Sunny.Vec3(1, 1, 1) / 8]
    types = [""]
    cryst = Sunny.crystal_from_spacegroup(latvecs, positions, types, cryst.sg; symprec=1e-8)
    ref_bonds = reference_bonds(cryst, 2.)
    dist2 = [Sunny.global_distance(cryst, b) for b in ref_bonds]

    # Using international symbol
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90) # must switch to standard cubic unit cell
    positions = [[1, 1, 1] / 4]
    @test_throws "Symbol \"F d -3 m\" is ambiguous; consider choice=\"2\" (standard) or choice=\"1\"" Crystal(latvecs, positions, "F d -3 m")
    cryst = Crystal(latvecs, positions, "F d -3 m"; choice="1")
    ref_bonds = reference_bonds(cryst, 2.)
    dist3 = [Sunny.global_distance(cryst, b) for b in ref_bonds]

    @test dist1 ≈ dist2 ≈ dist3

    ### FCC lattice, primitive vs. standard unit cell

    latvecs = [1 1 0; 0 1 1; 1 0 1]' / 2
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)

    latvecs = [1 0 0; 0 1 0; 0 0 1]'
    positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    cryst′ = Crystal(latvecs, positions)

    @test Sunny.get_wyckoff(cryst, 1) == Sunny.get_wyckoff(cryst′, 1)

    # Calculate interaction table
    ref_bonds = reference_bonds(cryst, 2.)
    bond = ref_bonds[2]
    basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, bond)
    bs = Sunny.all_symmetry_related_bonds_for_atom(cryst, bond.i, bond)
    @test length(bs) == Sunny.coordination_number(cryst, bond.i, bond)

    ### Triangular lattice, primitive unit cell

    c = 10
    latvecs = [1 0 0; -1/2 √3/2 0; 0 0 c]'
    @test collect(lattice_params(latvecs)) ≈ [1., 1., c, 90., 90., 120.]
    @test Sunny.cell_type(latvecs) == Sunny.hexagonal
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    @test Sunny.natoms(cryst) == 1
    @test Sunny.cell_volume(cryst) ≈ c * √3 / 2

    ### Kagome lattice

    latvecs = [1 0 0; -1/2 √3/2 0; 0 0 c]'
    positions = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
    cryst = Crystal(latvecs, positions)
    @test Sunny.natoms(cryst) == 3
    @test Sunny.cell_volume(cryst) ≈ c * √3 / 2

    ### Arbitrary monoclinic

    latvecs = lattice_vectors(6, 7, 8, 90, 90, 40)
    @test Sunny.cell_type(latvecs) == Sunny.monoclinic
    positions = [[0,0,0]]
    msg = "Symbol \"C2\" is ambiguous; consider \"A 1 1 2\" or \"B 1 1 2\" or \"I 1 1 2\""
    @test_throws msg Crystal(latvecs, positions, "C2")
    cryst = Crystal(latvecs, positions, "C 2/c"; choice="c1")
    @test Sunny.natoms(cryst) == 4
    @test_throws "Cell is nonstandard for spacegroup 5; consider \"A 1 1 2\" or \"B 1 1 2\" or \"I 1 1 2\"" Crystal(latvecs, positions, 5)
    @test_throws "Monoclinic axis choice (c) is incompatible with \"C 1 2 1\"" Crystal(latvecs, positions, "C 1 2 1")
    Crystal(latvecs, positions, "A 1 1 2") # No error
    Crystal(lattice_vectors(6, 7, 8, 90, 40, 90), positions, 5) # No error

    ### Arbitrary trigonal

    latvecs = lattice_vectors(5, 5, 6, 90, 90, 120)
    @test Sunny.cell_type(latvecs) == Sunny.hexagonal
    positions = [[0,0,0]]
    cryst1 = Crystal(latvecs, positions, "P -3")
    @test Sunny.natoms(cryst1) == 1
    cryst2 = Crystal(latvecs, positions, "R -3")
    @test Sunny.natoms(cryst2) == 3
    cryst3 = Crystal(latvecs, positions, 147) # same as "P -3"
    @test Sunny.natoms(cryst3) == 1
    msg = "Expected hexagonal or rhombohedral lattice system but got triclinic"
    @test_throws msg Crystal(lattice_vectors(5, 5, 6, 90, 95, 120), positions, 148)

    ### Arbitrary triclinic

    latvecs = lattice_vectors(6, 7, 8, 70, 80, 90)
    @test Sunny.cell_type(latvecs) == Sunny.triclinic
    positions = [[0,0,0]]
    cryst1 = Crystal(latvecs, positions, "P 1")
    @test Sunny.natoms(cryst1) == 1
    cryst2 = Crystal(latvecs, positions) # Infers 'P -1'
    @test Sunny.natoms(cryst1) == Sunny.natoms(cryst2) == 1

    ### Orthorhombic test, found by Ovi Garlea

    latvecs = lattice_vectors(13.261, 7.718, 6.278, 90.0, 90.0, 90.0)
    types = ["Yb1", "Yb2"]
    positions = [[0,0,0], [0.266,0.25,0.02]] # Locations of atoms as multiples of lattice vectors
    cryst = Crystal(latvecs, positions, 62; types, symprec=1e-4)
    @test count(==(1), cryst.classes) == 4
    @test count(==(2), cryst.classes) == 4
end


@testitem "Crystal errors" begin
    latvecs = lattice_vectors(1, 1, 1.5, 90, 90, 120)
    x = 0.15

    positions = [[x, 2x, 1/4], [-x, -2x, 3/4 + 1e-3]]
    msg = "Equivalent positions [0.1500, 0.3000, 1/4] and [-0.1500, -0.3000, 0.7510] in Wyckoff 6h at symprec=0.001"
    @test_throws msg Crystal(latvecs, positions, 194; symprec=1e-3)

    positions = [[x, 2x, 1/4], [-x, -2x, 3/4 + 2e-3]]
    msg = "Near-equivalent positions [0.1500, 0.3000, 1/4] and [-0.1500, -0.3000, 0.7520] in Wyckoff 6h at symprec=0.001"
    @test_throws msg Crystal(latvecs, positions, 194; symprec=1e-3)

    positions = [[x, 2x, 1/4], [-x, -2x, 3/4 + 5e-3]]
    Crystal(latvecs, positions, 194; symprec=1e-3) # No error

    positions = [[-x, -2x, 3/4], [-x, -2x, 3/4 + 1e-3]]
    msg = "Overlapping positions [-0.1500, -0.3000, 3/4] and [-0.1500, -0.3000, 0.7510] at symprec=0.001"
    @test_throws msg Crystal(latvecs, positions; symprec=1e-3)

    positions = [[-x, -2x, 3/4], [-x, -2x, 3/4 + 2e-3]]
    msg = "Near-overlapping positions [-0.1500, -0.3000, 3/4] and [-0.1500, -0.3000, 0.7520] at symprec=0.001"
    @test_throws msg Crystal(latvecs, positions; symprec=1e-3)

    positions = [[-x, -2x, 3/4], [-x, -2x, 3/4 + 5e-3]]
    Crystal(latvecs, positions; symprec=1e-3)
end


@testitem "Spacegroup settings" begin
    using LinearAlgebra
    import Spglib

    # Check conversions between settings for different Hall numbers
    for hall_c in 1:530
        hall_s = Sunny.standard_setting_for_hall_number(hall_c)
        setting = Sunny.mapping_to_standard_setting(hall_c)
        g_c = Sunny.SymOp.(Spglib.get_symmetry_from_database(hall_c)...)
        g_s = Sunny.SymOp.(Spglib.get_symmetry_from_database(hall_s)...)
        @test Sunny.map_symop_to_setting.(g_s; setting) ≈ g_c
    end

    ### Check settings for trigonal spacegroup

    # Trigonal spacegroup in standard hexagonal setting
    latvecs = lattice_vectors(1, 1, 1.2, 90, 90, 120)
    cryst = Crystal(latvecs, [[0, 0, 0]], 160)
    @test Sunny.get_wyckoff(cryst, 1).sitesym == "3m"

    # Same spacegroup in rhombohedral setting, which is the primitive cell
    prim_latvecs = cryst.latvecs * primitive_cell(cryst)
    cryst2 = Crystal(prim_latvecs, [[0, 0, 0]], 160; choice="R")
    @test primitive_cell(cryst2) ≈ I

    # Meaningful error if using spacegroup with non-default setting
    @test_throws "Cell is nonstandard for spacegroup 160; consider \"R3m\"" Crystal(prim_latvecs, [[0, 0, 0]], 160)
    @test_throws "Expected rhombohedral lattice system but got hexagonal" Crystal(latvecs, [[0, 0, 0]], 160; choice="R")

    # Check equivalence of positions
    @test norm(cryst.latvecs * cryst.positions[1]) < 1e-12
    @test cryst.latvecs * cryst.positions[2] ≈ cryst2.latvecs[:, 1]
    @test cryst.latvecs * cryst.positions[3] ≈ cryst2.latvecs[:, 1] + cryst2.latvecs[:, 2]

    # Inference of Wyckoff symbols
    lat_vecs = lattice_vectors(1, 1, 1.2, 90, 90, 120)
    cryst = Crystal(lat_vecs, [[0.2, 0.2, 1/2]], 164)
    @test Sunny.get_wyckoff(cryst, 1).letter == 'h'

    ### Check settings for monoclinic spacegroup

    # Standard setting for monoclinic spacegroup 5
    latvecs = lattice_vectors(1, 1.1, 1.2, 90, 100, 90)
    cryst = Crystal(latvecs, [[0, 0.2, 1/2]], "C 1 2 1")
    @test cryst.sg.label == "'C 2 = C 1 2 1' (5)"
    @test Sunny.get_wyckoff(cryst, 1).letter == 'b'

    # Alternative setting
    latvecs2 = reduce(hcat, eachcol(latvecs)[[3, 1, 2]])
    cryst2 = Crystal(latvecs2, [[1/2, 0, 0.2]], "A 1 1 2")
    @test cryst2.sg.label == "'C 2 = A 1 1 2' (5)"
    @test Sunny.get_wyckoff(cryst, 1).letter == 'b'

    # Verify `cryst` is already in standard setting
    @test cryst.sg.setting.R ≈ I

    # Equivalence of standardized lattice vectors
    @test cryst.latvecs ≈ cryst2.latvecs * inv(cryst2.sg.setting.R)

    # Primitive lattice vectors are alway standardized
    prim_latvecs1 = cryst.latvecs * primitive_cell(cryst)
    prim_latvecs2 = cryst2.latvecs * primitive_cell(cryst2)
    @test prim_latvecs1 ≈ prim_latvecs2
end


@testitem "Idealize setting and positions" begin
    latvecs = lattice_vectors(1, 1, 2, 90, 90, 120)
    positions = [[0.3333, 0.6667, 0.0]]
    cryst = Crystal(latvecs, positions, 191; symprec=1e-3)
    @test cryst.positions ≈ [[2/3, 1/3, 0], [1/3, 2/3, 0]]

    latvecs = lattice_vectors(1, 1, 2, 90, 90, 120)
    positions = [[0.6667, 0.3333, 0.0], [0.3333, 0.6667, 0.0]]
    cryst = Crystal(latvecs, positions; symprec=1e-3)
    @test cryst.positions ≈ [[2/3, 1/3, 0], [1/3, 2/3, 0]]

    latvecs = lattice_vectors(9.091, 9.091, 13.285, 90.0, 90.0, 120.0)
    z = 0.3333
    # Idealizes to [1/3, 2/3, 2/3 + z] where z is held at low precision
    positions = [[0.3333, 0.6667, 2/3 + z]] 
    cryst = Crystal(latvecs, positions, 166; symprec=1e-3)
    ref = [[2/3, 1/3, 1/3 - z], [0, 0, 0 + z], [1/3, 2/3, 2/3 - z], [2/3, 1/3, 1/3 + z], [0, 0, 0 - z + 1], [1/3, 2/3, 2/3 + z]]
    @test cryst.positions ≈ ref
end

@testitem "Conventional settings" begin
    misses = Int[]

    for hall in 1:530
        latvecs = if Sunny.cell_type(hall) == Sunny.hexagonal
            lattice_vectors(1, 1, 1, 90, 90, 120)
        else
            lattice_vectors(1, 1, 1, 90, 90, 90)
        end

        sgnum = Int(Sunny.all_spacegroup_types[hall].number)
        setting = Sunny.mapping_to_standard_setting(hall)

        sg = Sunny.Spacegroup(hall)
        setting′ = Sunny.conventionalize_setting(latvecs, setting, sgnum)

        # The `setting` defines the map from the given crystal cell to an ITA
        # "standard" one. Crystallographically, the choice is not unique:
        # left-applying any symop to `setting` would give another one that is
        # also valid. The setting convention _does_, however, matter for a
        # `cryst` loaded from an mCIF. Specifically, `cryst.sg.setting` will
        # define the indexing convention for a system that has been reshaped to
        # match the magnetic cell of that mCIF. For this reason, we need a fixed
        # and unambiguous convention. The function `conventionalize_setting`
        # includes heuristics that aim for consistency with the table
        # `mapping_to_standard_setting`, which was sourced from PyXTal.
        if !(sg.setting ≈ setting′)
            push!(misses, hall)
        end
    end

    # Any ITA standard setting must conventionalize to itself. Conversely, all
    # misses must be non-standard settings.
    @test all(misses) do hall
        Sunny.standard_setting_for_hall_number(hall) != hall
    end

    # The simple heuristics in `conventionalize_setting` miss for the 54 Hall
    # numbers listed below. These misses are to some extent arbitrary (depending
    # on the choices made in the PyXTal-sourced table). Attempts to reduce this
    # miss-count were only successful with significant complication to the
    # heuristics, which does not seem worthwhile. The purpose of this test is to
    # ensures consistency of the `conventionalize_setting` behavior.
    @test misses == [5, 11, 14, 15, 59, 65, 68, 69, 111, 114, 118, 121, 127, 136, 141, 153, 163, 172, 175, 184, 189, 195, 201, 207, 211, 214, 217, 220, 225, 232, 237, 240, 242, 243, 252, 255, 265, 277, 282, 299, 302, 305, 308, 312, 315, 317, 320, 324, 328, 330, 340, 344, 347, 527]
end

@testitem "Standardize Crystal" begin
    cryst = Crystal([1 0 1; 1 1 0; 0 1 1], [[0.1, 0.2, 0.3]])
    cryst2 = standardize(cryst)
    @test cryst2.positions ≈ [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]

    msg = "Found a nonconventional hexagonal unit cell. Consider using `lattice_vectors(a, a, c, 90, 90, 120)`."
    cryst = @test_logs (:warn, msg) Crystal(lattice_vectors(1, 1, 1, 90, 90, 60), [[0.1, 0.2, 0.3]])
    cryst2 = standardize(cryst)
    @test cryst2.positions ≈ [[0.0, 0.0, 0.0]]
end


@testitem "Allowed anisotropy" begin
    using LinearAlgebra

    # Test some inferred anisotropy matrices
    let
        # This test should also work for S = Inf, but there is a false negative
        # within second call to `is_anisotropy_valid`. TODO: Create minimized
        # test for `isapprox` bug and report to DynamicPolynomials repo. Cf.
        # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/created_by/kbarros

        s = 3
        k = 6
        i = 1
        cryst = Sunny.diamond_crystal()
        O = stevens_matrices(s)

        # print_site(cryst, i)
        Λ = O[6,0]-21O[6,4]
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)

        R = [normalize([1 1 -2]); normalize([-1 1 0]); normalize([1 1 1])]
        # print_site(cryst, i; R)
        Λ = O[6,0]-(35/√8)*O[6,3]+(77/8)*O[6,6]
        Λ′ = rotate_operator(Λ, R)
        @test Sunny.is_anisotropy_valid(cryst, i, Λ′)

        latvecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)

        msg = "Found a nonconventional tetragonal unit cell. Consider using `lattice_vectors(a, a, c, 90, 90, 90)`."
        cryst = @test_logs (:warn, msg) Crystal(latvecs, [[0, 0, 0]])

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
            sys = System(cryst, [1 => Moment(s=2, g=2)], mode)
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
            sys.dipoles .= [R*S for S in sys.dipoles]

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
    cryst = Crystal(Sunny.Mat3(I), [[0, 0, 0]], 227)
    @test cryst.positions ≈ [
        [0, 0, 0], [1/4, 1/4, 0], [1/2, 1/2, 0], [3/4, 3/4, 0], [1/4, 0, 1/4], [0, 1/4, 1/4], [3/4, 1/2, 1/4], [1/2, 3/4, 1/4], [1/2, 0, 1/2], [3/4, 1/4, 1/2], [0, 1/2, 1/2], [1/4, 3/4, 1/2], [3/4, 0, 3/4], [1/2, 1/4, 3/4], [1/4, 1/2, 3/4], [0, 3/4, 3/4],
    ]
    capt = IOCapture.capture() do
        print_symmetry_table(cryst, 0.8)
    end
    @test capt.output == """
        Atom 1
        Position [0, 0, 0], Wyckoff 16c
        Allowed g-tensor: [A B B
                           B A B
                           B B A]
        Allowed anisotropy in Stevens operators:
            c₁*(𝒪[2,-2]+2𝒪[2,-1]+2𝒪[2,1]) +
            c₂*(7𝒪[4,-3]+2𝒪[4,-2]-𝒪[4,-1]-𝒪[4,1]-7𝒪[4,3]) + c₃*(𝒪[4,0]+5𝒪[4,4]) +
            c₄*(11𝒪[6,-6]+8𝒪[6,-3]-𝒪[6,-2]+8𝒪[6,-1]+8𝒪[6,1]-8𝒪[6,3]) + c₅*(-𝒪[6,0]+21𝒪[6,4]) + c₆*(9𝒪[6,-6]+24𝒪[6,-5]+5𝒪[6,-2]+8𝒪[6,-1]+8𝒪[6,1]+24𝒪[6,5])

        Bond(1, 2, [0, 0, 0])
        Distance 0.3535533906, coordination 6
        Connects [0, 0, 0] to [1/4, 1/4, 0]
        Allowed exchange matrix: [A C -D
                                  C A -D
                                  D D  B]
        Allowed DM vector: [-D D 0]

        Bond(3, 5, [0, 0, 0])
        Distance 0.6123724357, coordination 12
        Connects [1/2, 1/2, 0] to [1/4, 0, 1/4]
        Allowed exchange matrix: [  A  C-E  D-F
                                  C+E    B -C+E
                                  D+F -C-E    A]
        Allowed DM vector: [E F -E]

        Bond(1, 3, [0, 0, 0])
        Distance 0.7071067812, coordination 6
        Connects [0, 0, 0] to [1/2, 1/2, 0]
        Allowed exchange matrix: [A D C
                                  D A C
                                  C C B]

        Bond(1, 3, [-1, 0, 0])
        Distance 0.7071067812, coordination 6
        Connects [0, 0, 0] to [-1/2, 1/2, 0]
        Allowed exchange matrix: [A D C
                                  D A C
                                  C C B]

        Bond(1, 2, [-1, 0, 0])
        Distance 0.7905694150, coordination 12
        Connects [0, 0, 0] to [-3/4, 1/4, 0]
        Allowed exchange matrix: [A  D -F
                                  D  B  E
                                  F -E  C]
        Allowed DM vector: [E F 0]

        """

    capt = IOCapture.capture() do
        print_site(cryst, 5; i_ref=2)
    end
    @test capt.output == """
        Atom 5
        Position [1/4, 0, 1/4], Wyckoff 16c
        Allowed g-tensor: [ A -B  B
                           -B  A -B
                            B -B  A]
        Allowed anisotropy in Stevens operators:
            c₁*(𝒪[2,-2]+2𝒪[2,-1]-2𝒪[2,1]) +
            c₂*(-7𝒪[4,-3]-2𝒪[4,-2]+𝒪[4,-1]-𝒪[4,1]-7𝒪[4,3]) + c₃*(𝒪[4,0]+5𝒪[4,4]) +
            c₄*(-11𝒪[6,-6]-8𝒪[6,-3]+𝒪[6,-2]-8𝒪[6,-1]+8𝒪[6,1]-8𝒪[6,3]) + c₅*(-𝒪[6,0]+21𝒪[6,4]) + c₆*(9𝒪[6,-6]+24𝒪[6,-5]+5𝒪[6,-2]+8𝒪[6,-1]-8𝒪[6,1]-24𝒪[6,5])
        """

    𝒪 = stevens_matrices(4)
    @test Sunny.is_anisotropy_valid(cryst, 5, 7𝒪[4,-3]+2𝒪[4,-2]-𝒪[4,-1]+𝒪[4,1]+7𝒪[4,3])

    capt = IOCapture.capture() do
        print_suggested_frame(cryst, 2)
    end
    @test capt.output == """
        R = [ 1/√2 -1/√2      0
             -1/√6 -1/√6 -√2/√3
              1/√3  1/√3  -1/√3]
        """

    R =  [1/√2 -1/√2 0; -1/√6 -1/√6 -√2/√3; 1/√3  1/√3  -1/√3]
    capt = IOCapture.capture() do
        print_site(cryst, 2; R)
    end
    @test capt.output == """
        Atom 2
        Position [1/4, 1/4, 0], Wyckoff 16c
        Allowed g-tensor: [A 0 0
                           0 A 0
                           0 0 B]
        Allowed anisotropy in Stevens operators:
            c₁*𝒪[2,0] +
            c₂*𝒪[4,-3] + c₃*𝒪[4,0] +
            c₄*𝒪[6,-3] + c₅*𝒪[6,0] + c₆*𝒪[6,6]
        Modified reference frame! Use R*g*R' or rotate_operator(op, R).
        """

    cryst = Sunny.hyperkagome_crystal()
    @assert Sunny.get_wyckoff(cryst, 1).letter == 'd'
    capt = IOCapture.capture() do
        print_suggested_frame(cryst, 2)
    end
    @test capt.output == """
        [ Info: Could not find a symmetry axis orthogonal to [1/√2, 1/√2, 0].
        R = [1/√2 -1/√2  0
                0     0 -1
             1/√2  1/√2  0]
        """

    # Test for https://github.com/SunnySuite/Sunny.jl/issues/260

    a = 6.22
    distortion = 0.15
    latvecs = lattice_vectors(a, a, a, 90+distortion, 90+distortion, 90+distortion)
    positions = Sunny.fcc_crystal().positions
    cryst = Crystal(latvecs, positions; types = ["A", "B", "B", "B"])

    capt = IOCapture.capture() do
        print_suggested_frame(cryst, 1)
    end
    @test capt.output == """
        R = [0.70803177573023 -0.70618057503467                 0
             0.40878233631266  0.40985392929053 -0.81542428107328
             0.57583678770556  0.57734630170186  0.57886375066688]
        """

    R = [0.70803177573023 -0.70618057503467                 0
         0.40878233631266  0.40985392929053 -0.81542428107328
         0.57583678770556  0.57734630170186  0.57886375066688]
    capt = IOCapture.capture() do
        print_site(cryst, 1; R)
        print_site(cryst, 2; R)
    end
    @test capt.output == """
        Atom 1
        Type 'A', position [0, 0, 0], Wyckoff 3a
        Allowed g-tensor: [A 0 0
                           0 A 0
                           0 0 B]
        Allowed anisotropy in Stevens operators:
            c₁*𝒪[2,0] +
            c₂*𝒪[4,-3] + c₃*𝒪[4,0] +
            c₄*𝒪[6,-3] + c₅*𝒪[6,0] + c₆*𝒪[6,6]
        Modified reference frame! Use R*g*R' or rotate_operator(op, R).
        Atom 2
        Type 'B', position [1/2, 1/2, 0], Wyckoff 9e
        Allowed g-tensor: [A   0   0
                           0   B D+E
                           0 D-E   C]
        Allowed anisotropy in Stevens operators:
            c₁*𝒪[2,-1] + c₂*𝒪[2,0] + c₃*𝒪[2,2] +
            c₄*𝒪[4,-3] + c₅*𝒪[4,-1] + c₆*𝒪[4,0] + c₇*𝒪[4,2] + c₈*𝒪[4,4] +
            c₉*𝒪[6,-5] + c₁₀*𝒪[6,-3] + c₁₁*𝒪[6,-1] + c₁₂*𝒪[6,0] + c₁₃*𝒪[6,2] + c₁₄*𝒪[6,4] + c₁₅*𝒪[6,6]
        Modified reference frame! Use R*g*R' or rotate_operator(op, R).
        """

    # These operators should be symmetry allowed
    @test Sunny.is_anisotropy_valid(cryst, 2, 𝒪[6,-1]+0.997385420𝒪[6,1])
    @test Sunny.is_anisotropy_valid(cryst, 2, rotate_operator(𝒪[6,2], R))
end


@testitem "Renormalization" begin
    latvecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)
    msg = "Found a nonconventional tetragonal unit cell. Consider using `lattice_vectors(a, a, c, 90, 90, 90)`."
    cryst = @test_logs (:warn, msg) Crystal(latvecs, [[0, 0, 0]])

    # Dipole system with renormalized anisotropy
    sys0 = System(cryst, [1 => Moment(s=3, g=2)], :dipole)
    randomize_spins!(sys0)

    i = 1
    O = stevens_matrices(spin_label(sys0, i))
    Λ = randn()*(O[2,0]+3O[2,2]) +
        randn()*(O[4,0]-5O[4,2]) + randn()*(O[4,0]+5O[4,4]) +
        randn()*(O[6,0]-21O[6,4]) + randn()*(O[6,0]+(105/16)O[6,2]+(231/16)O[6,6])
    set_onsite_coupling!(sys0, Λ, i)
    E0 = energy(sys0)

    # Corresponding SU(N) system
    sys = System(cryst, [1 => Moment(s=3, g=2)], :SUN)
    for site in eachsite(sys)
        set_dipole!(sys, sys0.dipoles[site], site)
    end
    set_onsite_coupling!(sys, Λ, i)
    E = energy(sys)

    @test E ≈ E0
end
