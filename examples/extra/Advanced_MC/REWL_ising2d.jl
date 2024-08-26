using Sunny, GLMakie

# Square lattice
latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
crystal = Crystal(latvecs, [[0,0,0]])

# 20×20 sites, ferromagnetic exchange
L = 20
sys = System(crystal, [1 => Moment(s=1, g=-1)], :dipole; dims=(L, L, 1), seed=0)
polarize_spins!(sys, (0,0,1))
set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

# REWL parameters - won't see much speedup vs. serial WL for small system
n_wins = 4
win_overlap = 0.8
windows = Sunny.get_windows((-2.0, 2.0), n_wins, win_overlap)
REWL = Sunny.ParallelWangLandau(; sys, bin_size=1/L^2, propose=propose_flip, windows)

# Sampling parameters
n_iters = 18
max_hchecks_per_iter = 100
hcheck_interval = 1_000
exch_interval = 100

# REWL sampling loop
for i in 1:n_iters
    for mcs in 1:max_hchecks_per_iter
        Sunny.step_ensemble!(REWL, hcheck_interval, exch_interval)
        # If flat, go to next iteration
        flat = fill(false, length(REWL.samplers))
        @Threads.threads for i in eachindex(REWL.samplers)
            flat[i] = Sunny.check_flat(REWL.samplers[i].hist; p=0.5)
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

# Merge into one ln_g(E)
E, ln_g = Sunny.merge(E_wins, ln_g_wins)
E .*= L^2

# Acceptance rates for replica exchanges
println("A = ", REWL.n_accept ./ REWL.n_exch)

# Derived thermodynamic properties
kT = collect(range(0.1, 20, length=1000))
U = Sunny.ensemble_average(E, ln_g, E, kT)
U² = Sunny.ensemble_average(E, ln_g, E .^2, kT)
C = @. (U² - U^2) / kT^2

# Plots
fig = Figure()
lines(fig[1,1], E/L^2, ln_g .- minimum(ln_g); axis=(xlabel="E/N", ylabel="ln[g(E)]"))
lines(fig[2,1], kT, U/L^2; axis=(xlabel="kT/J", ylabel="U/N"))
lines(fig[3,1], kT, C/L^2; axis=(xlabel="kT/J", ylabel="C/N"))
