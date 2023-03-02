@testitem "Spin operators" begin
    include("shared.jl")
    
    ### Verify 𝔰𝔲(2) irreps
    for N = 2:5
        S₀ = (N-1)/2
        S = Sunny.spin_matrices(N)

        for i in 1:3, j in 1:3
            # Test commutation relations
            @test S[i]*S[j] - S[j]*S[i] ≈ im * sum(ϵ[i,j,k]*S[k] for k=1:3)

            # Test orthonormality
            @test tr(S[i]*S[j]) ≈ (2/3)*S₀*(S₀+1/2)*(S₀+1)*δ(i,j)
        end

        # Test magnitude
        @test sum(S[i]^2 for i=1:3) ≈ S₀*(S₀+1)*I

        # Test dipole -> ket -> dipole round trip
        n = S₀ * normalize(randn(Sunny.Vec3))
        Z = Sunny.ket_from_dipole(n, Val(N))
        @test Sunny.expected_spin(Z) ≈ n

        # Test time reversal operator
        Z = randn(Sunny.CVec{N})
        @test Sunny.flip_ket(Z) ≈ exp(-im*π*S[2]) * conj(Z)
    end

    # Test action of apply_spin_matrices!(B, Z)
    for N = 4:6
        Λ = randn(ComplexF64, N, N)
        B = randn(Sunny.Vec3)
        Z = randn(Sunny.CVec{N})
        @test Sunny.mul_spin_matrices(Λ, B, Z) ≈ (Λ + B'*Sunny.spin_matrices(N)) * Z
    end    
end


@testitem "Stevens operators" begin
    include("shared.jl")
    import WignerSymbols: clebschgordan, wigner3j

    # Spherical tensors satisfying `norm(T) = √tr T† T = 1` (currently unused).
    function spherical_tensors_normalized(k; N)
        S = (N-1)/2
        ret = Matrix{Float64}[]
        for q in k:-1:-k
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
        S = Sunny.spin_matrices(N)
        Sp = S[1] + im*S[2]
        Sm = S[1] - im*S[2]
        
        for k in 0:N-1
            # Spherical tensors acting on N-dimensional Hilbert space
            T = spherical_tensors(k; N)

            # Generators of rotations in the spin-k representation
            K = Sunny.spin_matrices(2k+1)

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
            𝒪 = Sunny.stevens_matrices(k; N)
            T = spherical_tensors(k; N)

            # Check that Stevens operators are proper linear combination of
            # spherical tensors
            @test 𝒪 ≈ Sunny.stevens_α[k] * T
    
            # Check conversion of coefficients
            c = randn(2k+1)
            b = Sunny.transform_spherical_to_stevens_coefficients(k, c)
            @test transpose(c)*T ≈ transpose(b)*𝒪
        end
    end

    # Test decomposition of a random Hermitian matrix into Stevens coefficients
    let
        N = 7 # big enough to yield contributions at k=6
        A = randn(ComplexF64, N, N)
        A = A + A'

        c = Sunny.matrix_to_stevens_coefficients(A)

        acc = zeros(ComplexF64, N, N)
        acc += (tr(A)/N) * I
        for k in 1:6
            acc += c[k]' * Sunny.stevens_matrices(k; N)
        end
        @test acc ≈ A
    end
end


@testitem "Rotations" begin
    include("shared.jl")

    rng = Random.Xoshiro(0)
    R = Sunny.Mat3(Sunny.random_orthogonal(rng, 3; special=true))
    N = 5

    # Test axis-angle decomposition
    let
        (n, θ) = Sunny.axis_angle(R)
        @test 1 + 2cos(θ) ≈ tr(R)
        @test norm(n) ≈ 1
        @test R*n ≈ n

        # Rodrigues formula
        R2 = zeros(3,3)
        for i=1:3, j=1:3
            R2[i,j] = δ(i,j)*cos(θ) + (1-cos(θ))*n[i]*n[j] - sin(θ)*sum(ϵ[i,j,k]*n[k] for k=1:3)
        end
        @test R2 ≈ R
    end

    # Test that spin matrices rotate as vectors
    let
        S = Sunny.spin_matrices(N)
        @test R * S ≈ rotate_operator.(S, Ref(R))
    end

    # Test that Stevens coefficients rotate properly
    let 
        A = randn(ComplexF64, N, N)
        A = A + A'
        c = Sunny.matrix_to_stevens_coefficients(A)

        # Rotate coefficients directly
        c′1 = Sunny.rotate_stevens_coefficients.(c, Ref(R))

        # Rotate matrix and recalculate coefficients
        A′ = rotate_operator(A, R)
        c′2 = Sunny.matrix_to_stevens_coefficients(A′)

        @test c′1 ≈ c′2
    end

    # Test that a symbolic operators rotate properly
    let
        p = randn(3)'*𝒮 + randn(5)'*Sunny.stevens_operator_symbols[2]
        @test Sunny.operator_to_matrix(rotate_operator(p, R); N) ≈ rotate_operator(Sunny.operator_to_matrix(p; N), R)
    end        
end


@testitem "Symbolic operators" begin
    using LinearAlgebra

    # Internal conversion between spin and Stevens operators
    let
        J = randn(3,3)
        J = J+J'
        p = randn(3)'*𝒮 + 𝒮'*J*𝒮 +
            randn(11)' * Sunny.stevens_operator_symbols[5] +
            randn(13)' * Sunny.stevens_operator_symbols[6]
        cp1 = p |> Sunny.operator_to_classical_polynomial
        cp2 = p |> Sunny.operator_to_classical_stevens |> Sunny.operator_to_classical_polynomial
        @test cp1 ≈ cp2
    end


    # Test fast evaluation of Stevens functions
    let
        import DynamicPolynomials as DP

        s = randn(Sunny.Vec3)
        p = randn(5)' * Sunny.stevens_operator_symbols[2] + 
            randn(9)' * Sunny.stevens_operator_symbols[4] +
            randn(13)' * Sunny.stevens_operator_symbols[6]
        (_, c2, _, c4, _, c6) = Sunny.operator_to_classical_stevens_coefficients(p, 1.0)

        p_classical = Sunny.operator_to_classical_polynomial(p)
        grad_p_classical = DP.differentiate(p_classical, Sunny.spin_classical_symbols)

        E_ref = p_classical(Sunny.spin_classical_symbols => s)
        gradE_ref = [g(Sunny.spin_classical_symbols => s) for g = grad_p_classical]

        stvexp = Sunny.StevensExpansion(c2, c4, c6)
        E, gradE = Sunny.energy_and_gradient_for_classical_anisotropy(s, stvexp)

        @test E ≈ E_ref

        # Above, when calculating gradE_ref, the value X = |S|^2 is treated
        # as varying with S, such that dX/dS = 2S. Conversely, when calculating
        # gradE, the value X is treated as a constant, such that dX/dS = 0. In
        # practice, gradE will be used to drive spin dynamics, for which |S| is
        # constant, and the component of gradE parallel to S will be projected
        # out anyway. Therefore we only need agreement between the parts of
        # gradE and gradE_ref that are perpendicular to S.
        gradE_ref -= (gradE_ref⋅s)*s / (s⋅s) # Orthogonalize to s
        gradE -= (gradE⋅s)*s / (s⋅s)         # Orthogonalize to s
        @test gradE_ref ≈ gradE
    end
end


