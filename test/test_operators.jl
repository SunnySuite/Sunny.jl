@testitem "Spin operators" begin
    using LinearAlgebra

    # Levi-Civita symbol
    ϵ = [(i-j)*(j-k)*(k-i)/2 for i=1:3, j=1:3, k=1:3]

    ### Verify 𝔰𝔲(2) irreps
    for N = 2:5
        s = (N-1)/2
        S = spin_matrices(s)

        for i in 1:3, j in 1:3
            # Test commutation relations
            @test S[i]*S[j] - S[j]*S[i] ≈ im * sum(ϵ[i,j,k]*S[k] for k=1:3)

            # Test orthonormality
            @test tr(S[i]*S[j]) ≈ (2/3)*s*(s+1/2)*(s+1) * (i==j)
        end

        # Test magnitude
        @test sum(S[i]^2 for i=1:3) ≈ s*(s+1)*I

        # Test dipole -> ket -> dipole round trip
        n = s * normalize(randn(Sunny.Vec3))
        Z = Sunny.ket_from_dipole(n, Val(N))
        @test Sunny.expected_spin(Z) ≈ n

        # Test time reversal operator
        Z = randn(Sunny.CVec{N})
        @test Sunny.flip_ket(Z) ≈ exp(-im*π*S[2]) * conj(Z)
    end

    # Test action of mul_spin_matrices(B, Z)
    for N = 4:6
        Λ = randn(ComplexF64, N, N)
        B = randn(Sunny.Vec3)
        Z = randn(Sunny.CVec{N})
        @test Sunny.mul_spin_matrices(Λ, B, Z) ≈ (Λ + B'*spin_matrices((N-1)/2)) * Z
    end    
end

@testitem "Quadrupole operators" begin
    # Test expected_quadrupole generated function
    for N = 4:6
        Z = randn(Sunny.CVec{N})
        Qs = Sunny.stevens_matrices_of_dim(2; N)
        @test Sunny.expected_quadrupole(Z) ≈ Sunny.Vec5(real(Z'*Q*Z) for Q in Qs)
    end

    # Test action of mul_quadrupole_matrices(dE_dQ, Z)
    for N = 4:6
        dE_dQ = randn(Sunny.Vec5)
        Z = randn(Sunny.CVec{N})
        @test Sunny.mul_quadrupole_matrices(dE_dQ, Z) ≈ (dE_dQ' * Sunny.stevens_matrices_of_dim(2; N)) * Z
    end
end


@testitem "Stevens operators" begin
    include("shared.jl")
    import WignerSymbols: clebschgordan, wigner3j

    # Spherical tensors satisfying `norm(T) = √tr T† T = 1` (currently unused).
    function spherical_tensors_normalized(k; N)
        s = (N-1)/2
        ret = Matrix{Float64}[]
        for q in k:-1:-k
            T = zeros(Float64, N, N)
            for i = 1:N, i′ = 1:N
                m  = s - i + 1
                m′ = s - i′+ 1
                T[i, i′] = clebschgordan(s, m′, k, q, s, m) * sqrt((2k+1)/N)
            end
            push!(ret, T)
        end
        return ret
    end

    # KS/BCS spherical tensors T(k,q) as N×N matrices
    function spherical_tensors(k; N)
        j = (N-1)/2
        ret = Matrix{Float64}[]
        for q in k:-1:-k
            Tq = zeros(Float64, N, N)
            for i′ in 1:N, i in 1:N
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

    # Check transformation properties of spherical tensors
    for N in 2:7
        S = spin_matrices((N-1)/2)
        Sp = S[1] + im*S[2]
        Sm = S[1] - im*S[2]
        
        for k in 0:N-1
            # Spherical tensors acting on N-dimensional Hilbert space
            T = spherical_tensors(k; N)

            # Generators of rotations in the spin-k representation
            K = spin_matrices(k)

            # The selected basis is q ∈ [|k⟩, |k-1⟩, ... |-k⟩]. This function
            # converts from a q value to a 1-based index.
            idx(q) = k-q+1

            # A random axis-angle
            θ = randn(3)
            # Different representations of the same physical rotation
            D = exp(-im * θ' * K)
            U = exp(-im * θ' * S)

            for q in -k:k
                # Racah's commutation relations
                @test bracket(S[3], T[idx(q)]) ≈ q * T[idx(q)]
                q < +k && @test bracket(Sp, T[idx(q)]) ≈ sqrt((k-q)*(k+q+1)) * T[idx(q+1)]
                q > -k && @test bracket(Sm, T[idx(q)]) ≈ sqrt((k+q)*(k-q+1)) * T[idx(q-1)]

                # Wigner D matrix encodes rotation
                @test U' * T[idx(q)] * U ≈ (conj(D) * T)[idx(q)]
            end
        end
    end

    # Check mapping between spherical tensors and Stevens operators
    for N in 2:7
        for k in 1:N-1
            O = Sunny.stevens_matrices_of_dim(k; N)
            T = spherical_tensors(k; N)

            # Check that Stevens operators are proper linear combination of
            # spherical tensors
            @test O ≈ Sunny.stevens_α[k] * T
    
            # Check conversion of coefficients
            c = randn(2k+1)
            b = Sunny.transform_spherical_to_stevens_coefficients(k, c)
            @test transpose(c)*T ≈ transpose(b)*O
        end
    end

    # Test decomposition of a random Hermitian matrix into Stevens coefficients
    let
        N = 7 # big enough to yield contributions at k=6
        A = Hermitian(randn(ComplexF64, N, N))
        c = Sunny.matrix_to_stevens_coefficients(A)

        acc = zeros(ComplexF64, N, N)
        acc += (tr(A)/N) * I
        for k in 1:6
            acc += c[k]' * Sunny.stevens_matrices_of_dim(k; N)
        end
        @test acc ≈ A
    end
end


@testitem "Rotations" begin
    include("shared.jl")

    rng = Random.Xoshiro(0)
    R = Sunny.Mat3(Sunny.random_orthogonal(rng, 3; special=true))

    # Test axis-angle decomposition
    let
        (n, θ) = Sunny.matrix_to_axis_angle(R)
        @test 1 + 2cos(θ) ≈ tr(R)
        @test norm(n) ≈ 1
        @test R*n ≈ n
        @test R ≈ Sunny.axis_angle_to_matrix(n, θ)
    end

    # Test that spin matrices rotate as vectors
    let
        for s in (3, Inf)
            S = spin_matrices(s)
            for α in 1:3
                @test (R * S)[α] ≈ rotate_operator(S[α], R)
            end
        end
    end

    # Test that Stevens quadrupoles rotate correctly
    let
        for s in (3, Inf)
            O = stevens_matrices(s)
            
            # Cannot use [O[2, q] for q in 2:-1:-2] because:
            # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/149
            Q = [O[2, 2], O[2, 1], O[2, 0], O[2, -1], O[2, -2]]

            V = Sunny.operator_for_stevens_rotation(2, R)
            A = V * Q
            A′ = rotate_operator.(Q, Ref(R))

            # Cannot use A ≈ A′ because:
            # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/148
            for (Aq, Aq′) in zip(A, A′)
                @test Aq ≈ Aq′
            end
        end
    end

    # Test that Stevens coefficients rotate properly
    let
        N = 7
        A = Hermitian(randn(ComplexF64, N, N))
        c = Sunny.matrix_to_stevens_coefficients(A)

        # Rotate coefficients directly
        c′1 = Sunny.rotate_stevens_coefficients.(c, Ref(R))

        # Rotate matrix and recalculate coefficients
        A′ = rotate_operator(A, R)
        c′2 = Sunny.matrix_to_stevens_coefficients(A′)

        @test c′1 ≈ c′2
    end

    # Test rotation of 5×5 matrix of quadrupole-quadrupole interactions
    let
        R = Sunny.Mat3(Sunny.random_orthogonal(rng, 3; special=true))
        A = randn(5, 5)
        
        N = 3
        O = Sunny.stevens_matrices_of_dim(2; N)
        Oi, Oj = to_product_space(O, O)
        
        O′ = Sunny.rotate_operator.(O, Ref(R))
        Oi′, Oj′ = to_product_space(O′, O′)
        
        A′ = Sunny.transform_coupling_by_symmetry(Sunny.Mat5(A), R, true)
        @test Oi′' * A′ * Oj′ ≈ Oi' * A * Oj
        
        A″ = Sunny.transform_coupling_by_symmetry(Sunny.Mat5(A), R, false)
        @test Oj′' * A″ * Oi′ ≈ Oi' * A * Oj
    end

    # Test evaluation of the classical Stevens functions (i.e. spherical
    # harmonics) and their gradients
    let 
        using LinearAlgebra, FiniteDifferences, OffsetArrays

        # Random dipole and Stevens coefficients
        s = normalize(randn(Sunny.Vec3))
        c = map(OffsetArrays.OffsetArray(0:6, 0:6)) do k
            iseven(k) ? randn(2k+1) : zero(2k+1)
        end
        stvexp = Sunny.StevensExpansion(c)

        # Rotate dipole and Stevens coefficients
        s′ = R*s
        stvexp′ = Sunny.rotate_operator(stvexp, R)

        # Verify that the energy is the same regardless of which is rotated
        E1, _ = Sunny.energy_and_gradient_for_classical_anisotropy(s′, stvexp)
        E2, _ = Sunny.energy_and_gradient_for_classical_anisotropy(s, stvexp′)
        @test E1 ≈ E2

        # Verify that gradient agrees with finite differences
        _, gradE1 = Sunny.energy_and_gradient_for_classical_anisotropy(s, stvexp)
        f(s) = Sunny.energy_and_gradient_for_classical_anisotropy(s, stvexp)[1]
        gradE2 = grad(central_fdm(5, 1), f, s)[1]

        # When calculating gradE2, the value X = |S|^2 is treated as varying
        # with S, such that dX/dS = 2S. Conversely, when calculating gradE1, the
        # value X is treated as a constant, such that dX/dS = 0. In practice,
        # gradE will be used to drive spin dynamics, for which |S| is constant,
        # and the component of gradE parallel to S will be projected out anyway.
        # Therefore we only need agreement in the components perpendicular to S.
        gradE1 -= (gradE1⋅s)*s
        gradE2 -= (gradE2⋅s)*s
        @test gradE1 ≈ gradE2
    end
end

@testitem "Symbolics" begin
    using LinearAlgebra
    import IOCapture, OffsetArrays

    @test repr(stevens_matrices(Inf)[3,1]) == "-𝒮ˣ³ - 𝒮ʸ²𝒮ˣ + 4𝒮ᶻ²𝒮ˣ"

    capt = IOCapture.capture() do
        𝒪 = stevens_matrices(Inf)
        𝒮 = spin_matrices(Inf)
        Sunny.pretty_print_operator((1/4)𝒪[4,4] + (1/20)𝒪[4,0] + (3/5)*(𝒮'*𝒮)^2)
    end
    @test capt.output == "𝒮ˣ⁴ + 𝒮ʸ⁴ + 𝒮ᶻ⁴\n"

    capt = IOCapture.capture() do
        𝒮 = spin_matrices(Inf)
        print_stevens_expansion(𝒮[1]^4 + 𝒮[2]^4 + 𝒮[3]^4)
        print_stevens_expansion(0*𝒮[1])
        print_stevens_expansion(0*𝒮[1]+3)
    end
    @test capt.output == "(1/20)𝒪₄₀ + (1/4)𝒪₄₄ + (3/5)𝒮⁴\n0\n3\n"

    capt = IOCapture.capture() do
        S = spin_matrices(2)
        print_stevens_expansion(S[1]^4 + S[2]^4 + S[3]^4)
        print_stevens_expansion(0*S[1])
        print_stevens_expansion(0*S[1]+3*I)
    end
    @test capt.output == "(1/20)𝒪₄₀ + (1/4)𝒪₄₄ + 102/5\n0\n3\n"

    # Test Stevens coefficients extraction
    S = spin_matrices(Inf)
    O = stevens_matrices(Inf)
    S² = π
    p = S'*S * O[4, 2]
    c = Sunny.operator_to_stevens_coefficients(p, S²)
    @test iszero(c[1]) && iszero(c[2]) && iszero(c[3]) && iszero(c[5]) && iszero(c[6])
    @test c[4] ≈ [0.0, 0.0, S², 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Test round trip Stevens -> spin -> Stevens
    c_ref = map(OffsetArrays.OffsetArray(0:6, 0:6)) do k
        randn(2k+1)
    end
    p = sum(c_ref[k]'*Sunny.stevens_symbols[k] for k in 0:6)
    p = Sunny.expand_as_spin_polynomial(p)
    c = Sunny.operator_to_stevens_coefficients(p, 1.0)
    @test c ≈ c_ref
end
