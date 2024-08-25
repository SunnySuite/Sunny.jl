
@testitem "Tensors basic" begin
    Ni, Nj = (3, 4)
    Si0 = spin_matrices((Ni-1)/2)
    Sj0 = spin_matrices((Nj-1)/2)
    Si, Sj = Sunny.to_product_space(Si0, Sj0)

    # Basic properties of Kronecker product
    A1 = randn(ComplexF64, Ni, Ni)
    A2 = randn(ComplexF64, Ni, Ni)
    B1 = randn(ComplexF64, Nj, Nj)
    B2 = randn(ComplexF64, Nj, Nj)
    @test kron(A1, B1) * kron(A2, B2) ≈ kron(A1*A2, B1*B2)
    @test kron(A1, A2, B1, B2) ≈ kron(kron(A1, A2), kron(B1, B2))

    # Check transpose
    @test Sunny.reverse_kron(kron(A1, B1), Ni, Nj) ≈ kron(B1, A1)

    # Check factorization: S₁ˣ⊗S₂ˣ + S₁ˣ⊗S₂ʸ == S₁ˣ⊗(S₂ˣ+S₂ʸ)
    B = Si[1] * Sj[1] + Si[1] * Sj[2]
    D = Sunny.svd_tensor_expansion(B, Ni, Nj)
    @test length(D) == 1
    @test sum(kron(d...) for d in D) ≈ B

    # Check complicated SVD decomposition
    B = Si' * randn(3, 3) * Sj + (Si' * randn(3, 3) * Sj)^2
    D = Sunny.svd_tensor_expansion(B, Ni, Nj)
    @test length(D) == 9 # a nice factorization is always possible
    @test sum(kron(d...) for d in D) ≈ B
end


@testitem "General interactions" begin
    cryst = Sunny.diamond_crystal()
    sys = System(cryst, [1 => Moment(s=2, g=2)], :SUN; dims=(2, 2, 2))
    randomize_spins!(sys)
    
    J = 0.5
    K = 1.0
    Γ = 0.2
    D = 0.4
    J_exch = [J   Γ   -D;
              Γ   J   -D;
              D   D  J+K]
    bond = Bond(1, 2, [0, 0, 0])
    
    set_exchange!(sys, J_exch, bond)
    E = energy(sys)
    dE_dZ = Sunny.energy_grad_coherents(sys)
    
    set_pair_coupling!(sys, (Si, Sj) -> Si'*J_exch*Sj, bond; extract_parts=false)
    E′ = energy(sys)
    dE_dZ′ = Sunny.energy_grad_coherents(sys)
    
    @test E ≈ E′
    @test dE_dZ ≈ dE_dZ′

    set_pair_coupling!(sys, (Si, Sj) -> (Si'*Sj)^2, bond)
    E = energy_per_site(sys)
    dE_dZ = Sunny.energy_grad_coherents(sys)

    set_pair_coupling!(sys, (Si, Sj) -> (Si'*Sj)^2, bond; extract_parts=false)
    E′ = energy_per_site(sys)
    dE_dZ′ = Sunny.energy_grad_coherents(sys)

    @test E ≈ E′
    @test dE_dZ ≈ dE_dZ′
end
