using Sunny
using MPI

# run using the command: "mpiexec -n [n_procs] julia --project [PT_afm_heisenberg.jl]"

#######################################################
#   System setup
#######################################################

# Make lattice
crystal = Sunny.diamond_crystal(a=1.0)

# Interactions -- units of K           
J = Sunny.BOLTZMANN * 7.5413
interactions = [
    heisenberg(J, Bond(1, 3, [0,0,0])),
]

# Spin system
extent = (4,4,4)
S = 3/2
sys = SpinSystem(crystal, interactions, extent, [SiteInfo(1, 3/2)])
rand!(sys)

#######################################################
#   PT setup
#######################################################

kT_min = Sunny.BOLTZMANN * 6
kT_max = Sunny.BOLTZMANN * 7.33

# Temperature schedule: should take index (i) and number of processes (N) and return temperature (kT)
kT_sched(i, N) = 10 ^ (range(log10(kT_min), log10(kT_max), length=N))[i]

# Make replica for PT
replica = Replica(MetropolisSampler(sys, 1.0, 1))

# Run PT
run_parallel_temp!(
    replica, kT_sched;
    therm_mcs=10000, measure_interval=100, rex_interval=10,
    max_mcs=100000, bin_size=J, print_hist=true
)

