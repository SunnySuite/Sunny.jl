using Sunny, Plots

# create 2D Ising model
a = 1
latvecs = lattice_vectors(a,a,10a,90,90,90)
crystal = Crystal(latvecs, [[0,0,0]])

L = 20
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=1)], :dipole, units=Units.theory, seed=0)
polarize_spins!(sys, (0,0,1))

set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

# REWL parameters - won't see much speedup vs. serial WL for small system
n_wins = 4
win_overlap = 0.8
windows = Sunny.get_windows((-2.0, 2.0), n_wins, win_overlap)

# create type for REWL simulation
REWL = Sunny.ParallelWangLandau(; sys, bin_size=1/L^2, propose=propose_flip, windows)

# sampling parameters
n_iters = 12
max_hchecks_per_iter = 100
hcheck_interval = 10_000
exch_interval = 100

# start REWL sampling
for i in 1:n_iters
    for mcs in 1:max_hchecks_per_iter
        Sunny.sample!(REWL, hcheck_interval, exch_interval)

        # check flatness and start next iteration if satisfied
        flat = fill(false, length(REWL.samplers))
        @Threads.threads for i in eachindex(REWL.samplers)
            flat[i] = Sunny.check_flat(REWL.samplers[i].hist)
        end
        all(flat) && break
    end
    println("iteration $i complete.")
    @Threads.threads for sampler in REWL.samplers
        Sunny.reset!(sampler.hist)
        sampler.ln_f /= 2
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