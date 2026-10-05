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

@testitem "Planck noise generator" begin
    import FFTW: fft!, fftfreq
    using Random, Statistics

    # Generate a trajectory with the Planck noise generator
    function generate_pn_traj(N, dt, damping, kT; numburnin=1000)
        rng = Random.Xoshiro(0)
        cng = Sunny.PlanckNoiseGenerator(dt; kT, damping, dims=(1,1,1,1))
        buf = zeros(N)

        # Burnin
        for _ in 1:numburnin
            Sunny.step_pn!(rng, cng)
        end

        # Trajectory
        for i in 1:N
            Sunny.step_pn!(rng, cng)
            buf[i] = cng.ζ[1]
        end

        buf
    end

    # Hann window for Welch averaging.
    hann(L) = 0.5 .* (1 .- cos.(2π .* (0:L-1) ./ (L-1)))

    # Simple routine for calculating the two-sided Welch-averaged estimate of the
    # power spectrum.
    function welch_psd_twosided(x, dt; nperseg=1024, overlap=0.5)
        N = nperseg
        step = round(Int64, N * (1 - overlap))
        fs = 1/dt

        w = hann(N)
        specnorm = fs * sum(abs2, w)
        X = zeros(ComplexF64, N)
        pspec = zeros(N)

        count = 0
        for n in 1:step:(length(x)-N+1)
            @. X = x[n:n+N-1] * w
            fft!(X)
            @. pspec += abs2(X) / specnorm
            count += 1
        end
        pspec ./= count

        fs = fftfreq(N, fs)   # frequencies in cycles / unit time

        return fs, pspec
    end

    function test_planck_noise_power_spectrum()
        kTs = [0.1, 1.0, 10.0, 100.0]
        damping = 1.0

        for kT in kTs
            dt = 1/(10kT)   # Roughly, largest step side that maintains stability of filter dynamics
            nperseg = 1024  # Window size for Welch averaging estimate of power spectrum

            sig_pn = generate_pn_traj(1_000_000, dt, damping, kT)
            fs, pspec = welch_psd_twosided(sig_pn, dt; nperseg)

            ωs = 2π*fs[1:(nperseg÷2)]
            estspec = pspec[1:(nperseg÷2)]
            refspec = map(ω -> Sunny.planck_spectrum(ω, kT), ωs)

            # Examine the relative error at frequencies not too low (requires very
            # large trajectory) and not too high (relative error hard to measure
            # because of small denominator).
            ilo = findfirst(ω -> ω > kT/5, ωs)
            ihi = findfirst(ω -> ω > 5kT, ωs) - 1
            relerr = @. abs(estspec[ilo:ihi] - refspec[ilo:ihi]) / refspec[ilo:ihi]

            # The analytical filter spectrum agrees with the Planck spectrum to
            # within about 4% over this range. The remaining error is
            # statistical noise in the Welch estimate.
            @test maximum(relerr) < 0.15 && sum(relerr)/length(relerr) < 0.05

            # The variance is a convention-independent check of the spectral
            # normalization. A one-sided/two-sided mix-up would change it by a
            # factor of 2. For the Planck spectrum, Var ζ = ∫ |ω| n(|ω|) dω/2π =
            # (π/6) kT². Each noise-driven filter contributes (c/Ω)² exactly.
            (; c₁, c₂, Ω₁, Ω₂) = Sunny.planck_noise_params(kT)
            @test isapprox(var(sig_pn), (c₁/Ω₁)^2 + (c₂/Ω₂)^2; rtol=0.02)
            @test isapprox(var(sig_pn), (π/6)*kT^2; rtol=0.05)
        end
    end

    test_planck_noise_power_spectrum()

    # With kT = 0 the noise vanishes identically.
    cng = Sunny.PlanckNoiseGenerator(0.01; kT=0.0, damping=1.0, dims=(1,1,1,1))
    for _ in 1:100
        Sunny.step_pn!(Random.Xoshiro(0), cng)
    end
    @test all(iszero, cng.ζ)
end


@testitem "LangevinPlanck single-spin statistics" begin
    using LinearAlgebra

    # A classical spin in a field precesses at a single frequency ω₀, for any
    # tilt angle. With weak damping, the tilt is driven by the noise spectrum at
    # ω₀, so the tilt distribution is Boltzmann with the effective temperature ε
    # = S(ω₀), where S is the spectrum of the noise filters (≈ ω₀ n(ω₀)). This
    # holds exactly for all spin magnitudes s, not only in the harmonic limit.
    # Many independent spins improve the statistics.
    function mean_projection(s, kT; seed)
        cryst = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0, 0, 0]])
        sys = System(cryst, [1 => Moment(; s, g=2)], :dipole; dims=(8, 8, 8), seed)
        set_field!(sys, [0, 0, 0.5])
        polarize_spins!(sys, [0, 0, -1])
        ∇E, = Sunny.get_dipole_buffers(sys, 1)
        Sunny.set_energy_grad_dipoles!(∇E, sys.dipoles, sys)
        ω₀ = norm(∇E[1])

        damping = 0.05 / s          # Effective relaxation rate s⋅λ⋅ω₀ is weak
        dt = min(0.02/ω₀, 0.1/kT)
        integrator = LangevinPlanck(dt; damping, kT)
        relax = 1 / (s * damping * ω₀)
        for _ in 1:round(Int, 10relax/dt)
            step!(sys, integrator)
        end
        acc = 0.0
        nsteps = round(Int, 50relax/dt)
        for _ in 1:nsteps
            step!(sys, integrator)
            acc += -sum(S -> S[3], sys.dipoles) / (s * length(sys.dipoles))
        end
        return ω₀, acc / nsteps
    end

    langevin_function(x) = coth(x) - 1/x

    # Nonlinear quantum regime and near-classical regime. For comparison, the
    # classical (white noise) predictions differ by -35% and -5%, respectively.
    # Residual errors of order 1% are statistical, or come from finite damping.
    for (s, kT, seed, rtol) in ((1.0, 1.0, 0, 0.04), (10.0, 5.0, 1, 0.025))
        ω₀, proj = mean_projection(s, kT; seed)
        ε = Sunny.filter_spectrum(ω₀, collect(Sunny.planck_noise_params(kT)))
        @test isapprox(proj, langevin_function(s*ω₀/ε); rtol)
    end
end


@testitem "LangevinPlanck interface" begin
    cryst = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2, 2, 2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_field!(sys, [0, 0, 0.5])
    polarize_spins!(sys, [0, 0, -1])
    E₀ = energy(sys)

    # Changing the temperature updates the noise source
    integrator = LangevinPlanck(0.01; damping=0.1, kT=1.0)
    integrator.kT = 2.0
    @test integrator.noisesource.kT == 2.0
    @test integrator.noisesource.Ω₂ ≈ 2 * Sunny.planck_noise_params_dimensionless.Ω₂

    # Changing dt updates the noise source
    integrator.dt = 0.02
    @test integrator.noisesource.dt == 0.02

    # Noise buffers are sized on first use, and copy preserves parameters
    step!(sys, integrator)
    @test size(integrator.noisesource.ζ) == (3, size(sys.dipoles)...)
    integrator2 = copy(integrator)
    @test (integrator2.dt, integrator2.damping, integrator2.kT) == (0.02, 0.1, 2.0)

    # A timestep too large for the noise filters triggers a warning
    @test_logs (:warn, r"exceeds 0.1/kT") LangevinPlanck(0.1; damping=0.1, kT=10.0)

    # Keyword-only constructor requires dt to be set before stepping
    integrator3 = LangevinPlanck(; damping=0.1, kT=1.0)
    @test_throws "Set integration timestep" step!(sys, integrator3)
    integrator3.dt = 0.01
    step!(sys, integrator3)

    # At kT = 0 the ground state is stationary
    polarize_spins!(sys, [0, 0, -1])
    integrator4 = LangevinPlanck(0.01; damping=0.1, kT=0.0)
    for _ in 1:100
        step!(sys, integrator4)
    end
    @test energy(sys) ≈ E₀

    # A timestep suggestion accounts for the noise filters
    @test Sunny.suggest_timestep_aux(sys, LangevinPlanck(; damping=0.1, kT=100.0); tol=1e-2) <= 0.1/100

end


@testitem "LangevinPlanck SU(N) statistics" begin
    # Spin-1 with ℋ = D Sz² has levels (0, D, D). With weak damping, every
    # energy-changing transition is thermalized at the same effective
    # temperature ε = S(D), so the stationary state is exactly the classical
    # SU(3) Boltzmann distribution at ε. Its energy is known in closed form
    # [Dahlbom et al., PRB 106, 235154 (2022), Eq. 76]. At kT = D/2 the
    # classical (white noise) prediction is 80% larger, and an implementation
    # using colored complex-vector noise would be about 35% too hot.
    cryst = Crystal(lattice_vectors(1, 1, 1.5, 90, 90, 90), [[0, 0, 0]])
    D = 1.0
    kT = 0.5
    damping = 0.1
    sys = System(cryst, [1 => Moment(s=1, g=2)], :SUN; dims=(6, 6, 6), seed=0)
    set_onsite_coupling!(sys, S -> D * S[3]^2, 1)
    randomize_spins!(sys)

    function mean_energy(sys, integrator, relax)
        (; dt) = integrator
        for _ in 1:round(Int, 10relax/dt)
            step!(sys, integrator)
        end
        acc = 0.0
        nsteps = round(Int, 60relax/dt)
        for _ in 1:nsteps
            step!(sys, integrator)
            acc += energy_per_site(sys)
        end
        return acc / nsteps
    end
    E = mean_energy(sys, LangevinPlanck(0.02; damping, kT), 1/(damping*kT))

    ε = Sunny.filter_spectrum(D, collect(Sunny.planck_noise_params(kT)))
    E_exact = 2ε + (D^2/ε) / (1 - exp(D/ε) + D/ε)
    # At this damping there is a systematic deficit of about 3%, which vanishes
    # linearly as damping → 0.
    @test isapprox(E, E_exact; rtol=0.08)
end
