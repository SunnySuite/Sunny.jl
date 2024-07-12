using Sunny, GLMakie

# Square lattice
latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
crystal = Crystal(latvecs, [[0,0,0]])

# 20×20 sites, ferromagnetic exchange
L = 20
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=-1)], :dipole, seed=0)
polarize_spins!(sys, (0,0,1))
set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

# Temperature schedule for parallel tempering
kT_min = 0.5
kT_max = 10.0
n_replicas = 40 # can be larger than number of cores on machine
kT_sched = collect(range(kT_min, kT_max, length=n_replicas))

# Metropolis sampler with local spin flips
sampler = LocalSampler(; kT=0, propose=propose_flip)

# Parallel tempering sampler
PT = Sunny.ParallelTempering(sys, sampler, kT_sched)
n_therm = 1000
n_measure = 200
measure_interval = 10
exch_interval = 5

# Energy histograms for each PT replica
E_hists = [Sunny.Histogram() for _ in 1:PT.n_replicas]

# Initial equilibration
Sunny.step_ensemble!(PT, n_therm, exch_interval)

# Start PT simulation
for _ in 1:n_measure
    # Run some sweeps and replica exchanges
    Sunny.step_ensemble!(PT, measure_interval, exch_interval)

    # Measurements - assuming LocalSampler used
    for (j, sampler) in enumerate(PT.samplers)
        E_hists[j][sampler.ΔE] += 1
    end
end

# Acceptance rates for (rank, rank+1)
A = PT.n_accept ./ PT.n_exch

# Apply multiple histogram reweighting to PT data to get ln_g(E)
E, ln_g = Sunny.WHAM(E_hists, PT.kT_sched; n_iters=1000)

# Derived thermodynamic properties
kT = collect(range(kT_min, kT_max, length=1000))
U = Sunny.ensemble_average(E, ln_g, E, kT)
U² = Sunny.ensemble_average(E, ln_g, E.^2, kT)
C = @. abs( (U² - U^2) / kT^2 )

# Plots
fig = Figure()
lines(fig[1,1], E/L^2, ln_g .- minimum(ln_g); axis=(xlabel="E/N", ylabel="ln[g(E)]"))
lines(fig[2,1], kT, U/L^2; axis=(xlabel="kT/J", ylabel="U/N"))
lines(fig[3,1], kT, C/L^2; axis=(xlabel="kT/J", ylabel="C/N"))
