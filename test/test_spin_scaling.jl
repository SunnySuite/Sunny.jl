@testset "Spin Scaling" begin


function make_exchange_interactions()
    J  = 1.0   # Anti-ferro nearest neighbor
    J′ = -1.0  # Ferro next-nearest neighbor
    K  = 1.0   # Scale of Kitaev term
    Γ  = 0.0   # Off-diagonal exchange, not used
    J_exch = [J     Γ   0.0;
              Γ     J   0.0;
              0.0  0.0  J+K]
    return [exchange(J_exch, Bond(1, 2, [0,0,0])),
            heisenberg(J′, Bond(1, 1, [1,0,0]))]
end


function make_test_system_lld(;κ=1.0)
    cryst = Sunny.fcc_crystal()

    # Exchange interactions
    exchange_interactions = make_exchange_interactions()

    # Quartic anisotropy
    D = 1.0 
    Jquar = zeros(3,3,3,3)
    Jquar[1,1,1,1] = Jquar[2,2,2,2] = Jquar[3,3,3,3] = D
    quartic_interactions = [quartic_anisotropy(Jquar, i, "quartic") for i ∈ 1:4]

    interactions_all = vcat(exchange_interactions..., quartic_interactions...) 
    dims = (4,4,4)

    return SpinSystem(cryst,
                      interactions_all,
                      dims,
                      [SiteInfo(1, 0, 2*I(3), κ)]
    )
end


function make_test_system_gsd(; κ=1.0, N=2)
    cryst = Sunny.fcc_crystal()

    # Exchange interactions
    exchange_interactions = make_exchange_interactions()

    # Quartic anisotropy
    S = Sunny.gen_spin_ops(N)
    quartic_sun = [SUN_anisotropy(-S[3]^4, i, "quartic") for i ∈ 1:4]

    dims = (4,4,4)
    interactions_all = vcat(exchange_interactions..., quartic_sun...) 

    return SpinSystem(cryst,
                      interactions_all,
                      dims,
                      [SiteInfo(1, N, 2*I(3), κ)]
    )
end

function spin_magnitude_stability_tester(sys_maker, integrators, num_kappas)
    Δt = 0.01
    κs = 3.0 * rand(num_kappas) 
    for integrator in integrators
        for κ in κs
            sys = sys_maker(; κ)
            int = integrator(sys)
            rand!(sys)
            mags = norm.(sys._dipoles)
            for i ∈ 1:100
                evolve!(int, Δt)
            end
            @test mags ≈ norm.(sys._dipoles)
        end
    end
end

function test_spin_magnitude_stability()
    kT = 0.1
    α  = 0.1
    num_kappas = 3

    integrators_lld = [sys -> Sunny.LangevinHeunP(sys, kT, α),
                       sys -> Sunny.SphericalMidpoint(sys)]
    integrators_gsd = [sys -> Sunny.LangevinHeunPSUN(sys, kT, α),
                       sys -> Sunny.SchrodingerMidpoint(sys)]

    spin_magnitude_stability_tester(make_test_system_lld, integrators_lld, num_kappas)
    spin_magnitude_stability_tester(make_test_system_gsd, integrators_gsd, num_kappas)
end

test_spin_magnitude_stability()


function test_energy_scaling_lld()
    N = 0
    num_scalings = 3    # number of κs to try

    cryst = Sunny.fcc_crystal()
    dims = (4,4,4)
    J_quad = I(3) 
    J_quar = zeros(3,3,3,3)
    J_quar[3,3,3,3] = 1.0

    interactions_lld = [heisenberg(1.0, Bond(1,2,[0,0,0])),
                        quadratic_anisotropy(J_quad, 1, ""),
                        quartic_anisotropy(J_quar, 1, "")]
    powers_lld = [2, 2, 4]

    for (interaction, power) in zip(interactions_lld, powers_lld)
        κs = 5.0 * rand(num_scalings)
        for κ in κs

            # Get energy for system when κ=1.0
            sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1, N, 2*I(3), 1.0)])
            rand!(sys)
            E₀ = energy(sys)

            # Get energy for same configuration but with κ scaling
            S₀ = copy(sys._dipoles)
            sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1, N, 2*I(3), κ)])
            sys._dipoles .= S₀
            Sunny.normalize_dipoles!(sys)
            E₁ = energy(sys)

            @test (E₁/E₀) ≈ κ^power
        end
    end
end

test_energy_scaling_lld()


function test_energy_scaling_gsd()
    N = 3
    num_scalings = 3    # number of κs to try

    cryst = Sunny.fcc_crystal()
    dims = (4,4,4)

    S = Sunny.gen_spin_ops(N)
    interactions_gsd = [heisenberg(1.0, Bond(1,2,[0,0,0])),
                        SUN_anisotropy(S[3]^4, 1, "")]
    powers_gsd = [2, 1]

    for (interaction, power) in zip(interactions_gsd, powers_gsd)
        κs = 5.0 * rand(num_scalings)
        for κ ∈ κs
            sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1, N, 2*I(3), 1.0)])
            rand!(sys)
            E₀ = energy(sys)

            Z₀ = copy(sys._coherents)
            sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1, N, 2*I(3), κ)])
            sys._coherents .= Z₀
            Sunny.set_expected_spins!(sys)
            E₁ = energy(sys)

            @test (E₁/E₀) ≈ κ^power
        end
    end
end

test_energy_scaling_gsd()



end