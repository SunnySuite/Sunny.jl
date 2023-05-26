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

# create type for Wang-Landau simulation
WL = Sunny.WangLandau(; bin_size=1/L^2, bounds=(NaN, NaN), propose=propose_flip, ln_f=1.0)

# REWL parameters - won't see much speedup vs. serial WL for small system
n_wins = 4
win_overlap = 0.8
windows = Sunny.get_windows((-2.0, 2.0), n_wins, win_overlap)

# create type for REWL simulation
REWL = Sunny.ParallelWangLandau(sys, WL, windows)

# initialize systems to their respective bounds - use pad of 50 bins
max_mcs_init = 10_000
@Threads.threads for rank in 1:REWL.n_replicas
    Sunny.init_system!(REWL.systems[rank], REWL.samplers[rank], max_mcs_init; limit_pad=50/L^2)
end

# sampling parameters
n_iters = 10
max_hchecks_per_iter = 100
hcheck_interval = 10_000
exch_interval = 100
flat = zeros(Bool, REWL.n_replicas)

# start REWL sampling - 20 iterations (final ln_f <= 1e-6)
for i in 1:n_iters
    for mcs in 1:max_hchecks_per_iter
        Sunny.sample!(REWL, hcheck_interval, exch_interval)

        # check flatness and start next iteration if satisfied
        flat .= false
        @Threads.threads for rank in 1:REWL.n_replicas
            flat[rank] = Sunny.check_flat(REWL.samplers[rank].hist)
        end
        all(flat) && break
    end
    println("iteration $i complete.")
    @Threads.threads for rank in 1:REWL.n_replicas
        Sunny.reset!(REWL.samplers[rank].hist)
        REWL.samplers[rank].ln_f /= 2.0
    end
end

# REWL results - unmerged ln_g(E) from windows
E_wins = [Sunny.get_keys(wl.ln_g) for wl in REWL.samplers]
ln_g_wins = [Sunny.get_vals(wl.ln_g) for wl in REWL.samplers]

# merge into one ln_g(E)
E, ln_g = Sunny.merge(E_wins, ln_g_wins)
E .*= L^2

# acceptance rates for replica exchanges
println("A = ", REWL.n_accept ./ REWL.n_exch)

# calculate thermodynamics by reweighting WL results to temperatures
kT = collect(range(0.1, 20, length=1000))
U = Sunny.ensemble_average(E, ln_g, E, kT)
U² = Sunny.ensemble_average(E, ln_g, E .^2, kT)
C = (U² .- U .^2) ./ (kT .^2)

# plot density of states and thermodynamics
display( plot(E/L^2, ln_g .- minimum(ln_g), xlabel="E/N", ylabel="ln[g(E)]", legend=false) )
display( plot(kT, U/L^2, xlabel="kT/J", ylabel="U/N", legend=false) )
display( plot(kT, C/L^2, xlabel="kT/J", ylabel="C/N", legend=false) )