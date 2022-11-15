@testitem "Spin Scaling" begin
include("test_shared.jl")


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


function make_test_system_lld(; spin_rescaling=1.0)
    cryst = Sunny.fcc_crystal()

    # Exchange interactions
    exchange_interactions = make_exchange_interactions()

    # Quartic anisotropy
    D = 1.0 
    quartic_interactions = [anisotropy(D*(𝒮[1]^4+𝒮[2]^4+𝒮[3]^4), 1, "quartic")]

    interactions_all = [exchange_interactions..., quartic_interactions...]
    dims = (3,3,3)

    return SpinSystem(cryst,
                      interactions_all,
                      dims,
                      [SiteInfo(1; spin_rescaling)]
    )
end


function make_test_system_gsd(; spin_rescaling=1.0, N=2)
    cryst = Sunny.fcc_crystal()

    # Exchange interactions
    exchange_interactions = make_exchange_interactions()

    # Quartic anisotropy
    quartic_sun = anisotropy(-𝒮[3]^4, 1, "quartic")

    dims = (3,3,3)
    interactions_all = [exchange_interactions..., quartic_sun]

    return SpinSystem(cryst,
                      interactions_all,
                      dims,
                      [SiteInfo(1; N, spin_rescaling)]
    )
end

function spin_magnitude_stability_tester(sys_maker, integrators, num_rescalings)
    Δt = 0.01
    spin_rescalings = 3.0 * rand(num_rescalings) 
    for integrator in integrators
        for spin_rescaling in spin_rescalings
            sys = sys_maker(; spin_rescaling)
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
    num_rescalings = 2  

    cryst = Sunny.fcc_crystal()
    dims = (2,2,2)

    interactions_lld = [heisenberg(1.0, Bond(1,2,[0,0,0])),
                        anisotropy(+𝒮[1]^4+𝒮[2]^4+𝒮[3]^4, 1, "quartic")]
    powers_lld = [2, 4]

    for (interaction, power) in zip(interactions_lld, powers_lld)
        spin_rescalings = 5.0 * rand(num_rescalings)
        for spin_rescaling in spin_rescalings

            # Get energy for system when spin_rescaling=1.0
            sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1; N)])
            rand!(sys)
            E₀ = energy(sys)

            # Get energy for same configuration but with a spin rescaling 
            S₀ = copy(sys._dipoles)
            sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1; N, spin_rescaling)])
            sys._dipoles .= S₀
            Sunny.normalize_dipoles!(sys)
            E₁ = energy(sys)

            @test (E₁/E₀) ≈ spin_rescaling^power
        end
    end
end

test_energy_scaling_lld()


function test_energy_scaling_gsd()
    Ns = [5, 6]
    num_rescalings = 2    # number of rescalings to try

    cryst = Sunny.fcc_crystal()
    dims = (2,2,2)

    Λ = 𝒪[4,0]+5𝒪[4,4]

    interactions_gsd = [heisenberg(1.0, Bond(1,2,[0,0,0])),
                        anisotropy(Λ, 1)]
    powers_gsd = [2, 1]

    for N in Ns
        for (interaction, power) in zip(interactions_gsd, powers_gsd)
            spin_rescalings = 5.0 * rand(num_rescalings)
            for spin_rescaling ∈ spin_rescalings
                sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1; N)])
                rand!(sys)
                E₀ = energy(sys)

                Z₀ = copy(sys._coherents)
                sys = SpinSystem(cryst, [interaction], dims, [SiteInfo(1; N, spin_rescaling)])
                sys._coherents .= Z₀
                Sunny.set_expected_spins!(sys)
                E₁ = energy(sys)

                @test (E₁/E₀) ≈ spin_rescaling^power
            end
        end
    end
end

test_energy_scaling_gsd()

"""Generates a trajectory for a single spin in the presence of an 
external magnetic field. Rescales resulting spin magnitude so trajectories
with different scalings can be directly compared.
"""
function generate_scaled_zeeman_trajectory(spin_rescaling, θ, Δt; N=0, dur=10.0)
    cryst = Sunny.cubic_crystal()
    dims = (1,1,1)
    interactions = [external_field([0.0, 0.0, 10.0])]

    sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N, spin_rescaling)])

    spin = [0.0, sin(θ), cos(θ)] .* spin_rescaling 
    dpv = Sunny.DipoleView(sys)
    dpv[1] = Sunny.Vec3(spin)

    Integrator = N == 0 ? SphericalMidpoint : SchrodingerMidpoint
    integrator = Integrator(sys)

    numsteps = round(Int, dur/Δt) 
    ts = (0:numsteps) .* Δt
    S = zeros(Sunny.Vec3, numsteps+1)
    S[1] = sys._dipoles[1]

    for i in 1:numsteps
        evolve!(integrator, Δt)
        S[i+1] = sys._dipoles[1]
    end

    return (;
        xs = [S[1]/spin_rescaling for S ∈ S],
        ys = [S[1]/spin_rescaling for S ∈ S],
        zs = [S[1]/spin_rescaling for S ∈ S],
        ts
    ) 
end

"""Tests invariance of spin dynamics under spin rescaling 
in the presence of a Zeeman term. Tests both LLD and GSD. 
"""
function test_scaling_zeeman()
    Δt = 0.001
    θ = (π/4 - π/32)*rand() + π/32  # amount to tilt spin in zy-plane
    spin_rescaling = 3.0*rand()
    Ns = [0, 2]

    for N ∈ Ns
        (; xs) = generate_scaled_zeeman_trajectory(1.0, θ, Δt; N)
        xs_1 = xs
        (; xs) = generate_scaled_zeeman_trajectory(spin_rescaling, θ, Δt; N)
        xs_2 = xs

        rms = √sum( (xs_2 .- xs_1) .^2 )

        @test rms < 1e-10 
    end
end

test_scaling_zeeman()

"""Generate a trajectory for a system with only quadratic interactions. Results are rescaled 
so results with different spin magnitudes can be compared directly.
"""
function generate_scaled_quadratic_trajectory(spin_rescaling, Δt; N=0, dur=10.0)
    rng = Random.MersenneTwister(111)
    cryst = Sunny.cubic_crystal()
    dims = (4,4,3)
    interactions = [
        heisenberg(1.0, Bond(1,1,[1,0,0])),
        dipole_dipole()
    ]
    if N == 0   # "Quadratic anisotropy" only scales quadratically for old dynamics
        push!(interactions, quadratic_anisotropy(1.0*I(3), 1))
    end

    sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N, spin_rescaling)]; rng)
    rand!(sys)

    Integrator = N == 0 ? SphericalMidpoint : SchrodingerMidpoint
    integrator = Integrator(sys)

    numsteps = round(Int, dur/Δt) 
    ts = (0:numsteps) .* Δt
    S = zeros(Sunny.Vec3, numsteps+1)
    S[1] = sys._dipoles[1]

    for i in 1:numsteps
        evolve!(integrator, Δt/spin_rescaling)
        S[i+1] = sys._dipoles[1]
    end

    return (;
        xs = [S[1]/spin_rescaling for S ∈ S],
        ys = [S[1]/spin_rescaling for S ∈ S],
        zs = [S[1]/spin_rescaling for S ∈ S],
        ts 
    ) 
end

"""Test invariance of dynamics (with Hamiltonian that is quadratic in spins) under 
the rescaling of spin magnitudes.
"""
function test_scaling_quadratic()
    Δt = 0.01
    spin_rescaling = 3.0*rand()
    Ns = [0, 2]

    for N ∈ Ns
        (; xs) = generate_scaled_quadratic_trajectory(1.0, Δt; N)
        xs_1 = xs
        (; xs) = generate_scaled_quadratic_trajectory(spin_rescaling, Δt; N)
        xs_2 = xs

        rms = √sum( (xs_2 .- xs_1) .^2 )

        @test rms < 1e-8
    end
end

test_scaling_quadratic()


end
