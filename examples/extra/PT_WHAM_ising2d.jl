using Sunny, Plots

# create 2D Ising model
a = 1
latvecs = lattice_vectors(a,a,10a,90,90,90)
crystal = Crystal(latvecs, [[0,0,0]])

L = 20
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=1)], :dipole, units=Units.theory, seed=0)

set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

# start with randomized state
for site in all_sites(sys)
    polarize_spin!(sys, (0, 0, rand([-1,1])), site)
end

# temperature schedule for thermodynamics
kT_min = 0.5
kT_max = 10.0
n_replicas = 40 # can be larger than number of cores on machine
kT_sched = collect(range(kT_min, kT_max, length=n_replicas))

# use Metropolis sampler for Ising system
sampler = LocalSampler(; kT=kT_min, propose=propose_flip)

# initialize parallel tempering 
PT = ParallelTempering(sys, sampler, kT_sched)

# sampling parameters
n_therm = 1000
n_measure = 200
measure_interval = 10
exch_interval = 5

# energy histograms for each PT replica
E_hists = [Sunny.BinnedArray{Float64, Int64}() for _ in 1:PT.n_replicas]

# equilibration
sample!(PT, n_therm, exch_interval)

# start PT simulation
for i in 1:n_measure
    # run some sweeps and replica exchanges
    sample!(PT, measure_interval, exch_interval)

    # measurements - assuming LocalSampler used
    for (j, sampler) in enumerate(PT.samplers)
        E_hists[j][sampler.ΔE] += 1
    end
end

# acceptance rates for (rank, rank+1)
A = PT.n_accept ./ PT.n_exch

# apply multiple histogram reweighting to PT data to get ln_g(E)
E, ln_g = Sunny.WHAM(E_hists, PT.kT_sched; n_iters=20_000)

# calculate thermodynamics
U = Sunny.ensemble_average(E, ln_g, E, PT.kT_sched)
U² = Sunny.ensemble_average(E, ln_g, E .^2, PT.kT_sched)
C = abs.( (U² .- U .^2) ./ (PT.kT_sched .^2) )

# plot density of states and thermodynamics
display( plot(E/L^2, ln_g .- minimum(ln_g), xlabel="E/N", ylabel="ln[g(E)]", legend=false) )
display( plot(PT.kT_sched, U/L^2, xlabel="kT/J", ylabel="U/N", legend=false) )
display( plot(PT.kT_sched, C/L^2, xlabel="kT/J", ylabel="C/N", legend=false) )
