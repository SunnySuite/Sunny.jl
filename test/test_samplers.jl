## TODO: Add LocalSampler to the tests below

@testitem "Anisotropy" begin
    
    # Analytical mean energy for SU(3) model with Λ = D*(Sᶻ)^2
    function su3_mean_energy(kT, D)
        a = D/kT
        return D * (2 - (2 + 2a + a^2)*exp(-a)) / (a * (1 - (1+a)*exp(-a)))
    end 

    # Analytical mean energy for SU(5) model with Λ = D*((Sᶻ)^2-(1/5)*(Sᶻ)^4)
    function su5_mean_energy(kT, D)
        a = 4D/(5kT)
        return 4D*(exp(-a)*(-a*(a*(a*(a+4)+12)+24)-24)+24) / (5a*(exp(-a)*(-a*(a*(a+3)+6)-6)+6))
    end

    # Eliminate all spacegroup symmetries
    function asymmetric_crystal()
        latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
        positions = [[0,0,0]]
        Crystal(latvecs, positions, 1)
    end


    function su3_anisotropy_model(; L=20, D=1.0, seed)
        cryst = asymmetric_crystal()
        sys = System(cryst, [1 => Moment(s=1, g=2)], :SUN; dims=(L, 1, 1), seed)
        set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
        randomize_spins!(sys)

        return sys
    end

    function su5_anisotropy_model(; L=20, D=1.0, seed)
        cryst = asymmetric_crystal()
        sys = System(cryst, [1 => Moment(s=2, g=2)], :SUN; dims=(L, 1, 1), seed)
        randomize_spins!(sys)

        S = spin_matrices(spin_label(sys, 1))
        R = Sunny.random_orthogonal(sys.rng, 3; special=true)
        Λ = Sunny.rotate_operator(D*(S[3]^2-(1/5)*S[3]^4), R)
        set_onsite_coupling!(sys, Λ, 1)

        return sys
    end

    function thermalize!(sys, integrator, dur)
        numsteps = round(Int, dur/integrator.dt)
        for _ in 1:numsteps
            step!(sys, integrator)
        end
    end

    function calc_mean_energy(sys, integrator, dur)
        numsteps = round(Int, dur/integrator.dt)
        Es = zeros(numsteps)
        for i in 1:numsteps
            step!(sys, integrator)
            Es[i] = energy_per_site(sys)
        end
        sum(Es)/length(Es) 
    end

    function test_su3_anisotropy_energy()
        D = 1.0
        L = 20   # number of (non-interacting) sites
        damping = 1.0
        dt = 0.01
        kTs = [0.125, 0.5]
        thermalize_dur = 10.0
        collect_dur = 100.0

        sys = su3_anisotropy_model(; D, L, seed=0)
        heun = Langevin(dt; damping, kT=0)
        midpoint = Sunny.ImplicitMidpoint(dt; damping, kT=0)

        for integrator in [heun, midpoint]
            for kT in kTs
                integrator.kT = kT
                thermalize!(sys, integrator, thermalize_dur)
                E = calc_mean_energy(sys, integrator, collect_dur)
                E_ref = su3_mean_energy(kT, D)
                @test isapprox(E, E_ref; rtol=0.1)
            end
        end
    end

    test_su3_anisotropy_energy()
        

    function test_su5_anisotropy_energy()
        D = 1.0
        L = 20   # number of (non-interacting) sites
        damping = 0.1
        dt = 0.01
        kTs = [0.125, 0.5]
        thermalize_dur = 10.0
        collect_dur = 200.0

        sys = su5_anisotropy_model(; D, L, seed=0)
        heun = Langevin(dt; damping, kT=0)
        midpoint = Sunny.ImplicitMidpoint(dt; damping, kT=0)

        for integrator in [heun, midpoint]
            for kT ∈ kTs
                integrator.kT = kT
                thermalize!(sys, integrator, thermalize_dur)
                E = calc_mean_energy(sys, integrator, collect_dur)
                E_ref = su5_mean_energy(kT, D)
                @test isapprox(E, E_ref; rtol=0.1)
            end
        end
    end

    test_su5_anisotropy_energy()
end


# Test energy statistics of a two-site spin chain (LLD and GSD).
@testitem "Spin chain" begin

    # Consider a hypercube [0, 1]ᵏ (coordinates satisfying 0 ≤ xᵢ ≤ 1) and a hyperplane
    # defined by (x₁ + x₂ + ... xₖ) = α. The volume of the hypercube "beyond" this hyperplane
    # (∑ᵢ xᵢ > α) is given by
    #   V = 1 - ∑_{i=0..floor(α)} (-1)ⁱ binomial(k, i) (α - i)^k / k!
    # https://math.stackexchange.com/a/455711/660903
    # Taking the derivative of V with respect to α gives the (k-1)-dimensional "area" of
    # intersection between the hyperplane and the hypercube.
    function cubic_slice_area(α, k)
        sum([(-1)^i * binomial(k, i) * (α - i)^(k-1) / factorial(k-1) for i=0:floor(Int, α)])
    end

    # Energy distribution for an open-ended spin chain
    function P(E, kT; n=2, J=1.0)
        E_min = -J * max(1., n - 1.)
        return (2J)^(n-2) * cubic_slice_area((E - E_min)/2J, n-1) * exp(-E/kT) / (2kT * sinh(J/kT))^(n-1)
    end

    # Generates an empirical probability distribution from `data`.
    function empirical_distribution(data, numbins)
        N = length(data)
        lo, hi = minimum(data), maximum(data)
        Δ = (hi-lo)/numbins
        boundaries = collect(0:numbins) .* Δ .+ lo

        counts = zeros(Float64, numbins)
        for x in data
            idx = min(round(Int, floor((x - lo)/Δ) + 1), numbins)
            counts[idx] += 1.0
        end

        Ps = counts / N
        (; Ps, boundaries)
    end

    # Produces a discrete probability distribution from the continous one for
    # comparison with the empirical distribution
    function discretize_P(boundaries, kT; n=2, J=1.0, Δ = 0.001)
        numbins = length(boundaries) - 1
        Ps = zeros(numbins)
        for i in 1:numbins 
            Es = boundaries[i]:Δ:boundaries[i+1]
            Ps[i] = sum([P(E, kT; n, J)*Δ for E in Es])
        end
        Ps
    end

    # Generates a two-site spin chain spin system
    function two_site_spin_chain(; mode, seed)
        latvecs = lattice_vectors(1,1,1,90,90,90)
        cryst = Crystal(latvecs, [[0,0,0]])
        
        s = mode==:SUN ? 1/2 : 1
        κ = mode==:SUN ? 2 : 1
        sys = to_inhomogeneous(System(cryst, [1 => Moment(; s, g=2)], mode; dims=(2, 1, 1), seed))
        sys.κs .= κ
        set_exchange_at!(sys, 1.0, (1,1,1,1), (2,1,1,1); offset=(-1,0,0))
        randomize_spins!(sys)

        return sys
    end

    # Checks that the Langevin sampler produces the appropriate energy
    # distribution for a two-site spin chain.
    function test_spin_chain_energy()
        for mode in (:SUN, :dipole)
            sys = two_site_spin_chain(; mode, seed=0)

            damping = 0.1
            kT = 0.1
            dt = 0.02
            heun = Langevin(dt; damping, kT)
            # midpoint = Sunny.ImplicitMidpoint(dt; damping, kT)

            n_equilib = 1000
            n_samples = 2000
            n_decorr = 500

            for integrator in (heun,)

                # Initialize the Langevin sampler and thermalize the system
                for _ in 1:n_equilib
                    step!(sys, integrator)
                end

                # Collect samples of energy
                Es = Float64[]
                for _ in 1:n_samples
                    for _ in 1:n_decorr
                        step!(sys, integrator)
                    end
                    push!(Es, energy(sys))
                end

                # Generate empirical distribution and discretize analytical distribution
                n_bins = 10
                (; Ps, boundaries) = empirical_distribution(Es, n_bins)
                Ps_analytical = discretize_P(boundaries, kT) 

                # RMS error between empirical distribution and discretized analytical distribution
                @test isapprox(Ps, Ps_analytical, atol=0.05)
            end
        end
    end

    test_spin_chain_energy()
end
