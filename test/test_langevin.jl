@testitem "Langevin Dynamics" begin
    include("shared.jl")

    #= Test energy statistics for an SU(3) single ion problem with anisotropy. (GSD only.) =#

    "Analytical mean energy for SU(3) model with Λ = D*(Sᶻ)^2"
    function su3_mean_energy(kT, D)
        a = D/kT
        return D * (2 - (2 + 2a + a^2)*exp(-a)) / (a * (1 - (1+a)*exp(-a))) # - Λ₀
    end 

    "Analytical mean energy for SU(5) model with Λ = D*((Sᶻ)^2-(1/5)*(Sᶻ)^4)"
    function su5_mean_energy(kT, D)
        a = 4D/(5kT)
        return 4D*(exp(-a)*(-a*(a*(a*(a+4)+12)+24)-24)+24) / (5a*(exp(-a)*(-a*(a*(a+3)+6)-6)+6)) # - Λ₀
    end

    #= Generates an FeI2 cyrstal (Fe+ ions only). This crystal supports
    the anisotropies in the tests below =#
    function FeI2_crystal()
        a = b = 4.05012
        c = 6.75214
        lat_vecs = lattice_vectors(a, b, c, 90, 90, 120)
        basis_vecs = [[0,0,0]]
        Crystal(lat_vecs, basis_vecs, 164; setting="1")
    end


    function su3_anisotropy_model(; L=20, D=1.0, seed)
        N = 3
        Λ = D*𝒮[3]^2
        cryst = FeI2_crystal()
        interactions = [anisotropy(Λ, 1)]
        dims = (L,1,1)

        sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N)]; seed)
        rand!(sys)

        return sys
    end

    function su5_anisotropy_model(; L=20, D=1.0, seed)
        N = 5
        Λ = D*(𝒮[3]^2-(1/5)*𝒮[3]^4)
        
        cryst = FeI2_crystal()
        interactions = [anisotropy(Λ, 1)]
        dims = (L,1,1)

        sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N)]; seed)
        rand!(sys)

        return sys
    end

    function thermalize!(sys, integrator, dur)
        Δt = integrator.Δt
        numsteps = round(Int, dur/Δt)
        for _ in 1:numsteps
            step!(sys, integrator)
        end
    end

    function calc_mean_energy(sys, integrator, dur)
        L = size(sys.dipoles)[1]
        numsteps = round(Int, dur/integrator.Δt)
        Es = zeros(numsteps)
        for i in 1:numsteps
            step!(sys, integrator)
            Es[i] = energy(sys) / L
        end
        sum(Es)/length(Es) 
    end

    function test_su3_anisotropy_energy()
        D = 1.0
        L = 20   # number of (non-interacting) sites
        λ = 1.0
        Δt = 0.01
        kTs = [0.125, 0.5]
        thermalize_dur = 10.0
        collect_dur = 100.0
        seed = 111

        sys = su3_anisotropy_model(; D, L, seed)
        integrator = LangevinHeunP(0.0, λ, Δt)

        for kT ∈ kTs
            integrator.kT = kT
            thermalize!(sys, integrator, thermalize_dur)
            E = calc_mean_energy(sys, integrator, collect_dur)
            E_ref = su3_mean_energy(kT, D)

            #= No more than 5% error with respect to reference. =#
            @test abs(E - E_ref) < 0.05*E_ref    
        end
    end

    test_su3_anisotropy_energy()
        

    function test_su5_anisotropy_energy()
        D = 1.0
        L = 20   # number of (non-interacting) sites
        λ = 0.1
        Δt = 0.01
        kTs = [0.125, 0.5]
        thermalize_dur = 10.0
        collect_dur = 100.0
        seed = 111

        sys = su5_anisotropy_model(; D, L, seed)
        integrator = LangevinHeunP(0.0, λ, Δt)

        for kT ∈ kTs
            integrator.kT = kT
            thermalize!(sys, integrator, thermalize_dur)
            E = calc_mean_energy(sys, integrator, collect_dur)
            E_ref = su5_mean_energy(kT, D)

            #= No more than 5% error with respect to reference. =#
            @test abs(E - E_ref) < 0.05*E_ref    
        end
    end

    test_su5_anisotropy_energy()




    #= Test energy statistics of a two-site spin chain (LLD and GSD). =#

    function cubic_slice_area(α, k)
        V = sum([(-1)^i * binomial(k, i) * (α - i)^(k-1) / factorial(k-1) for i=0:floor(Int, α)])
    end

    "Energy distribution for an open-ended spin chain"
    function P(E, kT; n=2, J=1.0)
        E_min = -J * max(1., n - 1.)
        return (2J)^(n-2) * cubic_slice_area((E - E_min)/2J, n-1) * exp(-E/kT) / (2kT * sinh(J/kT))^(n-1)
    end

    "Generates an empirical probability distribution from `data`."
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

    "Produces a discrete probability distribution from the continous one for comparison
    with the empirical distribution"
    function discretize_P(boundaries, kT; n=2, J=1.0, Δ = 0.001)
        numbins = length(boundaries) - 1
        Ps = zeros(numbins)
        for i in 1:numbins 
            Es = boundaries[i]:Δ:boundaries[i+1]
            Ps[i] = sum([P(E, kT; n, J)*Δ for E ∈ Es])
        end
        Ps
    end

    "Generates a two-site spin chain spin system."
    function two_site_spin_chain(; N=0, J=1.0, spin_rescaling=1.0, seed)
        a = 1.0
        b = 1.1
        c = 1.2

        lat_vecs = lattice_vectors(a,b,c,90,90,90)
        basis_vecs = [[0,0,0], [0.45, 0.0, 0.0]]
        cryst = Crystal(lat_vecs, basis_vecs)
        interactions = [heisenberg(J, Bond(1,2,[0,0,0]))]
        dims = (1,1,1)
        sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N, spin_rescaling)]; seed)
        rand!(sys)

        return sys
    end

    "Checks that the Langevin sampler produces the appropriate energy distribution 
    for a two-site spin chain."
    function test_spin_chain_energy()
        Ns = [0, 2]
        spin_rescalings = [1.0, 2.0]
        seed = 111
        for (N, spin_rescaling) in zip(Ns, spin_rescalings)
            sys = two_site_spin_chain(; N, spin_rescaling, seed)

            λ = 0.1
            kT = 0.1
            Δt = 0.01

            n_decorr = 500  # Decorrelation steps between samples
            n_samples = 1000
            n_bins = 10  # Number of bins in empirical distribution

            # Initialize the Langevin sampler and thermalize the system
            integrator = LangevinHeunP(kT, λ, Δt)
            sampler = LangevinSampler(integrator, 1000)
            sample!(sys, sampler)
            sampler.nsteps = n_decorr

            # Collect samples of energy
            Es = zeros(n_samples)
            for i ∈ 1:n_samples
                sample!(sys, sampler)
                Es[i] = energy(sys)
            end

            # Generate empirical distribution and discretize analytical distribution
            (; Ps, boundaries) = empirical_distribution(Es, n_bins)
            Ps_analytical = discretize_P(boundaries, kT) 

            # RMS error between empirical distribution and discretized analytical distribution
            rms = sqrt(sum( (Ps .- Ps_analytical) .^ 2))

            @test rms < 0.05
        end
    end

    test_spin_chain_energy()

end