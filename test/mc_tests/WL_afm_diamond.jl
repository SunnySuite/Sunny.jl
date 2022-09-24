using Sunny

## Make lattice
crystal = Sunny.diamond_crystal(a=1.0)

## Interactions -- units of K           
J = 28.28 
interactions = [
    heisenberg(J, Bond{3}(3, 6, [0,0,0])),
]

## Spin system
extent = (4,4,4)
sys = SpinSystem(crystal, interactions, extent)
sys_size = length(sys)
rand!(sys)

## Wang-Landau sampler
N_bins = 1000
bs = 2*J / N_bins
E_min = -53.0
E_max = 0.0
ps = true

WL = WangLandau(
    system=sys, 
    bin_size = bs,
    hcheck_interval=10_000,
    bounds=[E_min, E_max], 
    mc_step_size=0.2,
    flatness=0.6,
    ln_f_final=1e-6,
    per_spin=ps,
    max_mcsweeps=10_000_000_000
)

## Run WL simulation
@time run!(WL)

