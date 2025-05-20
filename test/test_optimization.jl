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
    function simple_sys(; dims=(4, 4, 1), mode, seed, s)
        cryst = Crystal(lattice_vectors(1, 1, 2, 90, 90, 90), [[0, 0, 0]])
        sys = System(cryst, [1 => Moment(; s, g=2)], mode; dims, seed) 
        set_exchange!(sys, -1, Bond(1, 1, [1, 0, 0]))
        set_onsite_coupling!(sys, S -> -S[3]^2, 1)
        sys
    end

    s = 3/2
    is_z_polarized(sys) = all(d -> abs(d[3]) ≈ s, sys.dipoles)

    seed = 101
    sys_dip = simple_sys(; mode=:dipole, seed, s)
    sys_sun = simple_sys(; mode=:SUN, seed, s)

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
    function simple_sys(; dims=(4,4,1), mode, seed, s)
        cryst = Crystal(lattice_vectors(1,1,2,90,90,90), [[0,0,0]])
        sys = System(cryst, [1 => Moment(; s, g=2)], mode; dims, seed)
        set_exchange!(sys, -1, Bond(1,1,[1,0,0]))
        set_onsite_coupling!(sys, S -> -S[3]^2, 1)
        enable_dipole_dipole!(sys, 1.0)
        sys
    end

    s = 3/2

    seed = 101
    sys_dip = simple_sys(; mode=:dipole, seed, s)
    sys_sun = simple_sys(; mode=:SUN, seed, s)
    randomize_spins!(sys_dip)
    randomize_spins!(sys_sun)

    @test minimize_energy!(sys_dip) > 0
    @test minimize_energy!(sys_sun) > 0
end


@testitem "Vacancies" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
    cryst = Crystal(latvecs, [[0,0,0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(3, 3, 1), seed=1)
    set_exchange!(sys, 1, Bond(1, 1, [1, 0, 0]))
    sys2 = to_inhomogeneous(sys)
    sys2.κs[1, 1, 1, 1] = 0.01
    randomize_spins!(sys2)
    minimize_energy!(sys2)
    energy_per_site(sys2)
    @test energy_per_site(sys2) ≈ -0.8889
    
    set_vacancy_at!(sys2, (1, 1, 1, 1))
    randomize_spins!(sys2)
    minimize_energy!(sys2)
    @test energy_per_site(sys2) ≈ -8/9
end
