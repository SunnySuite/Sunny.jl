using Sunny, Plots

# create 2D Ising model
a = 1
latvecs = lattice_vectors(a,a,10a,90,90,90)
crystal = Crystal(latvecs, [[0,0,0]])

L = 20
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=1)], :dipole, units=Units.theory, seed=0)

set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

# start with randomized state
for i in 1:length(sys.dipoles)
    sys.dipoles[i] = Sunny.Vec3(0,0,rand([-1,1]))
end

# create type for Wang-Landau simulation
WL = WangLandau(;bin_size=1/L^2, bounds=[-2.0, 2.0], propose=propose_flip, ln_f=1.0)

# get system in bounded energy range - set pad to 50 bins
max_mcs_init = 10_000
Sunny.init_system!(sys, WL, max_mcs_init; limit_pad=50/L^2)

# sampling parameters
n_iters = 20
max_hchecks_per_iter = 100
hcheck_interval = 10_000

# start WL sampling - 20 iterations (final ln_f <= 1e-6)
for i in 1:n_iters
    # try to sample the histogram criterion (flat histogram)
    for mcs in 1:max_hchecks_per_iter
        sample!(sys, WL, hcheck_interval)

        # criterion satisfied - start next iteration
        if Sunny.check_flat(WL.hist; p=0.8)
            break
        end
    end
    println("iteration $i complete.")
    Sunny.reset!(WL.hist)
    # standard WL update to ln_f - can change to 1/t algorithm or other
    WL.ln_f /= 2.0
end

# WL results
E = Sunny.get_keys(WL.ln_g) .* L^2
ln_g = Sunny.get_vals(WL.ln_g)

# calculate thermodynamics by reweighting WL results to temperatures
kT = collect(range(0.1, 20, length=1000))
U = Sunny.ensemble_average(E, ln_g, E, kT)
U² = Sunny.ensemble_average(E, ln_g, E .^2, kT)
C = (U² .- U .^2) ./ (kT .^2)

# plot density of states and thermodynamics
display( plot(E/L^2, ln_g .- minimum(ln_g), xlabel="E/N", ylabel="ln[g(E)]", legend=false) )
display( plot(kT, U/L^2, xlabel="kT/J", ylabel="U/N", legend=false) )
display( plot(kT, C/L^2, xlabel="kT/J", ylabel="C/N", legend=false) )