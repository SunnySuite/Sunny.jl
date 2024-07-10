using Sunny, GLMakie

# Square lattice
latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
crystal = Crystal(latvecs, [[0,0,0]])

# 20×20 sites, ferromagnetic exchange
L = 20
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=-1)], :dipole, seed=0)
polarize_spins!(sys, (0, 0, 1))
set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

# WL parameters
WL = Sunny.WangLandau(; sys, bin_size=1/L^2, bounds=(-2.0, 2.0), propose=propose_flip)
n_iters = 18
max_hchecks_per_iter = 100
nsweeps = 1_000

# WL sampling loop
for i in 1:n_iters
    for mcs in 1:max_hchecks_per_iter
        Sunny.step_ensemble!(WL, nsweeps)
        # If flat, go to next iteration
        if Sunny.check_flat(WL.hist; p=0.5)
            break
        end
    end
    println("iteration $i complete.")
    Sunny.reset!(WL.hist)
    # Can change to 1/t algorithm or other
    WL.ln_f /= 2
end

# WL density of states
E = Sunny.get_keys(WL.ln_g) .* L^2
ln_g = Sunny.get_vals(WL.ln_g)

# Derived thermodynamic properties
kT = collect(range(0.1, 8, length=1000))
U = Sunny.ensemble_average(E, ln_g, E, kT)
U² = Sunny.ensemble_average(E, ln_g, E.^2, kT)
C = @. (U² - U^2) / kT^2

# Plots
fig = Figure()
lines(fig[1,1], E/L^2, ln_g .- minimum(ln_g); axis=(xlabel="E/N", ylabel="ln[g(E)]"))
lines(fig[2,1], kT, U/L^2; axis=(xlabel="kT/J", ylabel="U/N"))
lines(fig[3,1], kT, C/L^2; axis=(xlabel="kT/J", ylabel="C/N"))
