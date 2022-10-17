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
    basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, b)
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


@testitem "Spin matrices" begin
    include("test_shared.jl")
    
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
        ψ = Sunny._get_coherent_from_dipole(n, Val(N))
        @test Sunny.expected_spin(ψ) ≈ n
    end    
end


@testitem "Spherical tensors" begin
    include("test_shared.jl")

    # Spherical tensors that satisfy `norm(T) =  √ tr T† T = 1`.
    function spherical_tensors_normalized(N, k)
        S = (N-1)/2
        ret = Matrix{Float64}[]
        for q = k:-1:-k
            T = zeros(Float64, N, N)
            for i = 1:N, i′ = 1:N
                m  = S - i + 1
                m′ = S - i′+ 1
                T[i, i′] = clebschgordan(S, m′, k, q, S, m) * sqrt((2k+1)/N)
            end
            push!(ret, T)
        end
        return ret
    end

    # Spherical tensors T(k,q) as NxN matrices. The result is ambiguous up to an
    # overall (k,N)-dependent scaling factor. Here we're using the normalization
    # convention of KS/BCS.
    function spherical_tensors(N, k)
        j = (N-1)/2
        ret = Matrix{Float64}[]
        for q = k:-1:-k
            Tq = zeros(Float64, N, N)
            for i′ = 1:N, i = 1:N
                m′ = j - i′+ 1
                m  = j - i + 1

                # By the Wigner-Eckhardt theorem, the spherical tensor T must have
                # this m and m′ dependence. An overall (j, k)-dependent rescaling
                # factor is arbitrary, however.
                Tq[i′, i] = (-1)^(j-m′) * wigner3j(j, k, j, -m′, q, m)
            end

            # Below we will apply two rescaling factors obtained from Rudowicz and
            # Chung, J. Phys.: Condens. Matter 16 (2004) 5825–5847.

            # With this rescaling factor, we get the Buckmaster and Smith & Thornley
            # (BST) operator
            Tq .*= 2.0^(-k) * sqrt(factorial((N-1)+k+1) / factorial((N-1)-k))

            # With this additional rescaling factor, we get the Koster and Statz
            # (1959) and Buckmaster et al (1972) operator (KS/BCS)
            Tq ./= sqrt(factorial(2k) / (2^k * factorial(k)^2))

            push!(ret, Tq)
        end
        return ret
    end

    # Lie bracket, aka matrix commutator
    bracket(A, B) = A*B - B*A

    for N=2:7
        S = Sunny.gen_spin_ops(N)
        Sp = S[1] + im*S[2]
        Sm = S[1] - im*S[2]
        
        for k = 0:N-1
            # Spherical tensors acting on N-dimensional Hilbert space
            T = spherical_tensors(N, k)

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

    # Stevens operators
    for N=2:7
        for k = 1:N-1
            𝒪 = Sunny.stevens_matrices(N, k)
            T = spherical_tensors(N, k)

            # Check that Stevens operators are proper linear combination of
            # spherical tensors
            @test 𝒪 ≈ Sunny.stevens_α[k] * T
    
            # Check conversion of coefficients
            c = randn(2k+1)
            b = Sunny.transform_spherical_to_stevens_coefficients(k, c)
            @test transpose(c)*T ≈ transpose(b)*𝒪
        end
    end
end

@testitem "Local operator symbols" begin
    include("test_shared.jl")

    A = randn(3,3)
    R = Sunny.Mat3(exp(A - A'))
    N = 5

    # Test axis-angle decomposiiton
    let
        (n, θ) = Sunny.axis_angle(R)
        @test 1 + 2cos(θ) ≈ tr(R)
        @test norm(n) ≈ 1
        @test R*n ≈ n
    end

    # Test that Stevens operator symbols transform properly
    let
        p = randn(3)' * Sunny.stevens_operator_symbols[1] +
            randn(5)' * Sunny.stevens_operator_symbols[2] +
            randn(7)' * Sunny.stevens_operator_symbols[3]
        @test Sunny.operator_to_matrix(rotate_operator(p, R); N) ≈ rotate_operator(Sunny.operator_to_matrix(p; N), R)
    end

    # Test that spin operator symbols transform properly
    let
        J = randn(3, 3)
        J = (J+J')/2
        p = randn(3)'*𝒮 + 𝒮'*J*𝒮
        @test Sunny.operator_to_matrix(rotate_operator(p, R); N) ≈ rotate_operator(Sunny.operator_to_matrix(p; N), R)
    end

    # Test that a linear combination transforms properly
    let
        p = randn(3)'*𝒮 + randn(5)'*Sunny.stevens_operator_symbols[2]
        @test Sunny.operator_to_matrix(rotate_operator(p, R); N) ≈ rotate_operator(Sunny.operator_to_matrix(p; N), R)
    end

    # Internal conversion between spin and Stevens operators
    let
        J = randn(3,3)
        J = J+J'
        p = randn(3)'*𝒮 + 𝒮'*J*𝒮 +
            randn(11)' * Sunny.stevens_operator_symbols[5] +
            randn(13)' * Sunny.stevens_operator_symbols[6]
        cp = Sunny.operator_to_classical_polynomial(p)
        @test cp |> Sunny.classical_polynomial_to_classical_stevens |> Sunny.operator_to_classical_polynomial ≈ cp
    end

    # Test some inferred anisotropy matrices
    let
        N = 7
        k = 6
        i = 1
        cryst = Sunny.diamond_crystal()

        # print_site(cryst, i)
        Λ = 𝒪[6,0]-21𝒪[6,4]
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)

        R = hcat(normalize([1, 1, -2]), normalize([-1, 1, 0]), normalize([1, 1, 1]))
        R = Sunny.Mat3(R)
        # print_site(cryst, i; R)
        Λ = 𝒪[6,0]-(35/√8)*𝒪[6,3]+(77/8)*𝒪[6,6]
        Λ′ = rotate_operator(Λ, R)
        @test Sunny.is_anisotropy_valid(cryst, i, Λ′)

        lat_vecs = lattice_vectors(1.0, 1.1, 1.0, 90, 90, 90)
        cryst = Crystal(lat_vecs, [[0., 0., 0.]])
        # print_site(cryst, i)
        Λ = randn()*(𝒪[6,0]-21𝒪[6,4]) + randn()*(𝒪[6,2]+(16/5)*𝒪[6,4]+(11/5)*𝒪[6,6])
        @test Sunny.is_anisotropy_valid(cryst, i, Λ)
    end

    # Test fast evaluation of Stevens operators
    let
        import DynamicPolynomials

        s = randn(Sunny.Vec3)
        p = randn(5)' * Sunny.stevens_operator_symbols[2] + 
            randn(9)' * Sunny.stevens_operator_symbols[4] +
            randn(13)' * Sunny.stevens_operator_symbols[6]
        (_, c2, _, c4, _, c6) = Sunny.operator_to_classical_stevens_coefficients(p)

        p_classical = Sunny.operator_to_classical_polynomial(p)
        grad_p_classical = DynamicPolynomials.differentiate(p_classical, Sunny.spin_classical_symbols)

        E_ref = p_classical(Sunny.spin_classical_symbols => s)

        gradE_ref = [g(Sunny.spin_classical_symbols => s) for g = grad_p_classical]
        gradE_ref -= (gradE_ref⋅s)*s / (s⋅s) # Orthogonalize to s

        E, gradE = Sunny.energy_and_gradient_for_classical_anisotropy(s, c2, c4, c6)
        gradE -= (gradE⋅s)*s / (s⋅s)         # Orthogonalize to s

        @test E ≈ E_ref
        @test gradE_ref ≈ gradE
    end
end
