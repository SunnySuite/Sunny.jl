using Sunny
using MPI

# run using the command: "mpiexec -n [n_procs] julia --project [PT_afm_heisenberg.jl]"

#######################################################
#   system setup

# make lattice
crystal = Sunny.diamond_crystal()

# interactions -- units of K  
J = 28.28
interactions = [
    heisenberg(J, Bond{3}(1, 3, [0,0,0])),
]

# spin system -- setup to use K
S = 1 
extent = (4,4,4)
sys = SpinSystem(crystal, interactions, extent, [SiteInfo(1, S, 1)]; μB=0.67171381563034223582, μ0=17.3497470317891588)
sys_size = length(sys)
rand!(sys)

#######################################################
#   PT setup

# make replica for PT
replica = Replica(MetropolisSampler(sys, 1.0, 1))

# temperature schedule: should take index (i) and number of processes (N) and return temperature (T)
T_min = 3.0
T_max = 50.0

# units of K
T_sched(i, N) = (10 .^(range(log10(T_min), log10(T_max), length=N)) )[i]

# run parallel tempering
run_PT!(
    replica,
    T_sched;
    therm_mcs=5000,
    measure_interval=100,
    rex_interval=10,
    max_mcs=1_000_000,
    bin_size=J,
    print_hist=true
)
