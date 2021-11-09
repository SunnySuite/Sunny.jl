using Sunny

# run using the command: "mpiexec -n [n_procs] julia --project [PT_afm_heisenberg.jl]"

#######################################################
#   system setup

# make lattice
crystal = Sunny.diamond_crystal(a=1.0)

# interactions -- units of K           
J = 28.28
interactions = [
    heisenberg(J, Bond{3}(3, 6, [0,0,0])),
]

# spin system
extent = (4,4,4)
sys = SpinSystem(crystal, interactions, extent)
sys_size = length(sys)
rand!(sys)

#######################################################
#   PT setup

T_min = 22.5
T_max = 27.5

# temperature schedule: should take index (i) and number of processes (N) and return temperature (T)
T_sched(i, N) = (10 .^(range(log10(T_min), log10(T_max), length=N)) )[i]

# make replica for PT
replica = Replica(MetropolisSampler(sys, 1.0, 1))

# run PT
run!(replica, T_sched; therm_mcs=10000, measure_interval=1, rex_interval=100, max_mcs=10000, bin_size=J)

