@testitem "Crystal Construction" begin
    include("test_shared.jl")


    ### Test construction of diamond lattice

    # Spglib inferred symmetry
    lat_vecs = [1 1 0; 1 0 1; 0 1 1]' / 2
    positions = [[1, 1, 1], [-1, -1, -1]] / 8
    cryst = Crystal(lat_vecs, positions)
    ref_bonds = reference_bonds(cryst, 2.)
    dist1 = [distance(cryst, b) for b in ref_bonds]

    # Using explicit symops
    lat_vecs = Sunny.Mat3(lat_vecs)
    positions = [Sunny.Vec3(1, 1, 1) / 8]
    types = [""]
    cryst = Sunny.crystal_from_symops(lat_vecs, positions, types, cryst.symops, cryst.spacegroup)
    ref_bonds = reference_bonds(cryst, 2.)
    dist2 = [distance(cryst, b) for b in ref_bonds]

    # Using Hall number
    lat_vecs = lattice_vectors(1, 1, 1, 90, 90, 90) # must switch to standard cubic unit cell
    positions = [Sunny.Vec3(1, 1, 1) / 4]
    cryst = Sunny.crystal_from_hall_number(lat_vecs, positions, types, 525)
    ref_bonds = reference_bonds(cryst, 2.)
    dist3 = [distance(cryst, b) for b in ref_bonds]

    # Using international symbol
    positions = [[1, 1, 1] / 4]
    # cryst = Crystal(lat_vecs, positions, "F d -3 m") # Ambiguous!
    cryst = Crystal(lat_vecs, positions, "F d -3 m"; setting="1")
    ref_bonds = reference_bonds(cryst, 2.)
    dist4 = [distance(cryst, b) for b in ref_bonds]

    @test dist1 ≈ dist2 ≈ dist3 ≈ dist4



    ### FCC lattice, primitive vs. standard unit cell

    lat_vecs = [1 1 0; 1 0 1; 0 1 1]' / 2
    positions = [[0, 0, 0]]
    cryst = Crystal(lat_vecs, positions)

    lat_vecs = [1 0 0; 0 1 0; 0 0 1]'
    positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    cryst′ = Crystal(lat_vecs, positions)

    @test cryst.sitesyms[1] == cryst′.sitesyms[1]

    # Calculate interaction table
    ref_bonds = reference_bonds(cryst, 2.)
    b = ref_bonds[2]
    basis = basis_for_symmetry_allowed_couplings(cryst, b)
    J = basis' * randn(length(basis))
    (bs, Js) = all_symmetry_related_couplings_for_atom(cryst, b.i, b, J)
    @test length(Js) == coordination_number(cryst, b.i, b)


    ### Triangular lattice, primitive unit cell

    c = 10
    lat_vecs = [1 0 0;  -1/2 √3/2 0;  0 0 c]'
    positions = [[0, 0, 0]]
    cryst = Crystal(lat_vecs, positions)
    @test cell_type(cryst) == Sunny.hexagonal
    @test nbasis(cryst) == 1
    @test cell_volume(cryst) ≈ c * √3 / 2 
    @test all(lattice_params(cryst) .≈ (1., 1., c, 90., 90., 120.))

    ### Kagome lattice

    lat_vecs = [1 0 0;  -1/2 √3/2 0;  0 0 c]'
    positions = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
    cryst = Crystal(lat_vecs, positions)
    @test cell_type(cryst) == Sunny.hexagonal
    @test nbasis(cryst) == 3
    @test cell_volume(cryst) ≈ c * √3 / 2 
    @test all(lattice_params(cryst) .≈ (1., 1., c, 90., 90., 120.))


    ### Arbitrary monoclinic

    mono_lat_params = (6, 7, 8, 90, 90, 40)
    lat_vecs = lattice_vectors(mono_lat_params...)
    positions = [[0,0,0]]
    # cryst = Crystal(lat_vecs, positions, "C 2/c")
    cryst = Crystal(lat_vecs, positions, "C 2/c", setting="c1")
    @test cell_type(cryst) == Sunny.monoclinic
    @test nbasis(cryst) == 4
    @test all(lattice_params(cryst) .≈ mono_lat_params)


    ### Arbitrary trigonal

    lat_vecs = lattice_vectors(5, 5, 6, 90, 90, 120)
    positions = [[0,0,0]]
    cryst1 = Crystal(lat_vecs, positions, "P -3")
    @test nbasis(cryst1) == 1
    @test cell_type(cryst1) == Sunny.hexagonal
    cryst2 = Crystal(lat_vecs, positions, "R -3")
    @test nbasis(cryst2) == 3
    cryst3 = Crystal(lat_vecs, positions, 147) # spacegroup number
    @test cell_type(cryst1) == cell_type(cryst2) == cell_type(cryst3) == Sunny.hexagonal


    ### Arbitrary triclinic

    lat_vecs = lattice_vectors(6, 7, 8, 70, 80, 90)
    positions = [[0,0,0]]
    cryst1 = Crystal(lat_vecs, positions, "P 1")
    @test nbasis(cryst1) == 1
    cryst2 = Crystal(lat_vecs, positions) # Infers 'P -1'
    @test nbasis(cryst1) == nbasis(cryst2) == 1
    @test cell_type(cryst1) == cell_type(cryst2) == Sunny.triclinic

    ### Orthorhombic test, found by Ovi Garlea

    lat_vecs = lattice_vectors(13.261, 7.718, 6.278, 90.0, 90.0, 90.0);
    types = ["Yb1","Yb2"];
    basis_vecs = [[0,0,0], [0.266,0.25,0.02]]; # Locations of atoms as multiples of lattice vectors
    crysts = Crystal(lat_vecs, basis_vecs, 62; types, symprec=1e-4)
    @test length(crysts) == 6
    cryst = Crystal(lat_vecs, basis_vecs,62; types, symprec=1e-4, setting="-cba")
    @test count(==(1), cryst.classes) == 4
    @test count(==(2), cryst.classes) == 4    
end


### TODO: Merge this with David's code
@testitem "Spin operators" begin
    include("test_shared.jl")

    function infer_ket_from_dipole(S, n::Sunny.Vec3)
        # Find a ket (up to an irrelevant phase) that corresponds to a pure dipole.
        # TODO, we can do this much faster by using the exponential map of spin
        # operators, expressed as a polynomial expansion,
        # http://www.emis.de/journals/SIGMA/2014/084/
        (evals, evecs) = eigen(n'*S)
        return normalize(evecs[:, argmax(evals)])
    end
    function spin_bilinear(S, Z)
        return Sunny.Vec3(real(Z'*S[1]*Z), real(Z'*S[2]*Z), real(Z'*S[3]*Z))
    end
    
    # Levi-Civita symbol
    ϵ = [(i-j)*(j-k)*(k-i)/2 for i=1:3, j=1:3, k=1:3]

    # Kronecker delta
    δ(i,j) = (i==j) ? 1 : 0

    ### Verify 𝔰𝔲(2) irreps
    for N = 2:5
        S₀ = (N-1)/2
        S = Sunny.gen_spin_ops(N)

        for i=1:3, j=1:3
            # Test commutation relations
            @test S[i]*S[j] - S[j]*S[i] ≈ im * sum(ϵ[i,j,k]*S[k] for k=1:3)

            # Test orthonormality
            @test tr(S[i]*S[j]) ≈ (2/3)*S₀*(S₀+1/2)*(S₀+1)*δ(i,j)
        end

        # Test magnitude
        @test sum(S[i]^2 for i=1:3) ≈ S₀*(S₀+1)*I

        # Test dipole -> ket -> dipole round trip
        n = S₀ * normalize(randn(Sunny.Vec3))
        ψ = infer_ket_from_dipole(S, n)
        @test spin_bilinear(S, ψ) ≈ n
    end    
end


@testitem "Spherical tensors" begin
    include("test_shared.jl")

    # Lie bracket, aka matrix commutator
    bracket(A, B) = A*B - B*A

    for N=2:7
        S = Sunny.gen_spin_ops(N)
        Sp = S[1] + im*S[2]
        Sm = S[1] - im*S[2]
        
        for k = 0:N-1
            # Spherical tensors acting on N-dimensional Hilbert space
            T = Sunny.spherical_tensors(N, k)

            # Generators of rotations in the spin-k representation
            K = Sunny.gen_spin_ops(2k+1)

            # The selected basis is q ∈ [|k⟩, |k-1⟩, ... |-k⟩]. This function
            # converts from a q value to a 1-based index.
            idx(q) = k-q+1

            # A random axis-angle
            θ = randn(3)
            # Different representations of the same physical rotation
            D = exp(-im * θ' * K)
            U = exp(-im * θ' * S)

            for q = -k:k
                # Racah's commutation relations
                @test bracket(S[3], T[idx(q)]) ≈ q * T[idx(q)]
                q < +k && @test bracket(Sp, T[idx(q)]) ≈ sqrt((k-q)*(k+q+1)) * T[idx(q+1)]
                q > -k && @test bracket(Sm, T[idx(q)]) ≈ sqrt((k+q)*(k-q+1)) * T[idx(q-1)]

                # Wigner D matrix encodes rotation
                @test U' * T[idx(q)] * U ≈ (conj(D) * T)[idx(q)]
            end
        end
    end
end

@testitem "Stevens operators" begin
    include("test_shared.jl")

    for N=2:7
        for k = 0:N-1
            𝒪 = Sunny.stevens_ops(N, k)
            T = Sunny.spherical_tensors(N, k)

            # Check that two ways of calculating Stevens operators agree
            @test 𝒪 ≈ Sunny.stevens_ops_alt(N, k)

            # Check conversion of coefficients
            c = randn(2k+1)
            b = Sunny.transform_spherical_to_stevens_coefficients(k, c)
            A1 = sum(c[i]*T[i] for i in eachindex(c))
            A2 = sum(b[q]*𝒪[q] for q in eachindex(b))
            @test A1 ≈ A2
        end
    end

    # Test some inferred anisotropy matrices
    begin
        N = 7
        k = 6
        i = 1
        cryst = Sunny.diamond_crystal()

        # print_allowed_anisotropy(cryst, i; k)
        𝒪 = Sunny.stevens_ops(N, k)
        Λ = 𝒪[0]-21𝒪[4]
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)

        R = hcat(normalize([1, 1, -2]), normalize([-1, 1, 0]), normalize([1, 1, 1]))
        R = Sunny.Mat3(R)
        # print_allowed_anisotropy(cryst, i; R)
        𝒪 = Sunny.stevens_ops(N, k; R)
        Λ = 𝒪[0]-(35/√8)*𝒪[3]+(77/8)*𝒪[6]
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)

        lat_vecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)
        cryst = Crystal(lat_vecs, [[0., 0., 0.]])
        # print_allowed_anisotropy(cryst, i)
        𝒪 = Sunny.stevens_ops(N, k)
        Λ = randn()*(𝒪[0]-21𝒪[4]) + randn()*(𝒪[2]+(16/5)*𝒪[4]+(11/5)*𝒪[6])
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)
    end
end

@testitem "Rotations" begin
    include("test_shared.jl")

    A = randn(3,3)
    R = exp(A - A')
    (n, θ) = Sunny.axis_angle(Sunny.Mat3(R))
    @test 1 + 2cos(θ) ≈ tr(R)
    @test norm(n) ≈ 1
    @test R*n ≈ n
end
