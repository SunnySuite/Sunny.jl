@testitem "Stereographic projection" begin
    using LinearAlgebra, FiniteDifferences

    # Test gradients (limited to real-valued vectors due to FiniteDifferences)
    let
        N = 10
        α = randn(Float64, N)
        n = normalize(randn(Float64, N))

        u(α) = Sunny.stereographic_projection(α, n)

        x = randn(Float64, N)
        x̄J  = Sunny.vjp_stereographic_projection(x, α, n)
        x̄J′ = x'*jacobian(central_fdm(5, 1), u, α)[1]
        @test x̄J ≈ x̄J′
    end

    # Test inverse projection
    for T in (Sunny.Vec3, Sunny.CVec{10})
        n = normalize(randn(T))
        v = (I - n*n') * randn(T)
        u = Sunny.stereographic_projection(v, n)
        @test Sunny.inverse_stereographic_projection(u, n) ≈ v
    end
end


@testitem "Optimization" begin
    # H = -∑Sᵢ⋅Sⱼ - ∑(Sᵢᶻ)² on 2D square lattice (z-polarized ground state)
    function simple_sys(; dims=(4,4,1), mode, seed, S)
        cryst = Crystal(lattice_vectors(1,1,2,90,90,90), [[0,0,0]])
        sys = System(cryst, dims, [SpinInfo(1; S, g=2)], mode; seed) 
        set_exchange!(sys, -1, Bond(1,1,[1,0,0]))
        set_onsite_coupling!(sys, S -> -S[3]^2, 1)
        sys
    end

    S = 3/2
    is_z_polarized(sys) = all(s -> abs(s[3]) ≈ S, sys.dipoles)

    seed = 101
    sys_dip = simple_sys(; mode=:dipole, seed, S)
    sys_sun = simple_sys(; mode=:SUN, seed, S)

    # Thermalize near ground state
    dt = 0.05
    integrator = Langevin(dt; damping=0.1, kT=0.1)
    for _ in 1:1000
        step!(sys_dip, integrator)
        step!(sys_sun, integrator)
    end
    @test minimize_energy!(sys_dip) > 0
    @test minimize_energy!(sys_sun) > 0

    @test is_z_polarized(sys_dip)
    @test is_z_polarized(sys_sun)

    # Randomize initial condition and use `minimize_energy!` 
    randomize_spins!(sys_dip)
    randomize_spins!(sys_sun)

    @test minimize_energy!(sys_dip) > 0
    @test minimize_energy!(sys_sun) > 0

    @test is_z_polarized(sys_dip)
    @test is_z_polarized(sys_sun)
end

@testitem "Optimization Coverage" begin
    # Make sure optimization works on system with dipole-dipole
    function simple_sys(; dims=(4,4,1), mode, seed, S)
        cryst = Crystal(lattice_vectors(1,1,2,90,90,90), [[0,0,0]])
        sys = System(cryst, dims, [SpinInfo(1; S, g=2)], mode; seed) 
        set_exchange!(sys, -1, Bond(1,1,[1,0,0]))
        set_onsite_coupling!(sys, S -> -S[3]^2, 1)
        enable_dipole_dipole!(sys)
        sys
    end

    S = 3/2

    seed = 101
    sys_dip = simple_sys(; mode=:dipole, seed, S)
    sys_sun = simple_sys(; mode=:SUN, seed, S)
    randomize_spins!(sys_dip)
    randomize_spins!(sys_sun)

    @test minimize_energy!(sys_dip) > 0
    @test minimize_energy!(sys_sun) > 0
end