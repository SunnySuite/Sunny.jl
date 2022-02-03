using Sunny
using LinearAlgebra
using MPI

# MPI rank is 0-based index
rank, N_ranks = init_MPI()

# Use kT as control parameter (α)
kT_min = 0.5
kT_max = 5.0

# Geometric kT schedule
#kT = kT_min * (kT_max/kT_min) ^ ((rank)/(N_ranks-1))

# Linear kT schedule
kT = collect(range(kT_min, kT_max, length=N_ranks))[rank+1]

# Nearest-neighbor ferromagnetic Ising (z-dir) interactions
J = -1.0
FM = exchange(diagm([J, J, J]), Bond(1, 1, [1, 0, 0]))
interactions = [FM]

# SC lattice -- extend lattice vector in z-direction to emulate 2D
lvecs = Sunny.lattice_vectors(1.0, 1.0, 2.0, 90, 90 ,90)
bvecs = [[0.0, 0.0, 0.0]]
crystal = Crystal(lvecs, bvecs)

# Make system and randomize spins
extent = (20, 20, 1)
system = SpinSystem(crystal, interactions, extent, [SiteInfo(1,1,1)])
randflips!(system)

# Make replica for REMC
α = kT
replica = Replica(IsingSampler(system, kT, 1), α)

# Run Replica Exchange MC (as PT)
run_REMC!(
    replica;
    therm_mcs=5000, 
    measure_interval=50, 
    rex_interval=10, 
    max_mcs=100_000, 
    bin_size=1.0, 
    print_hist=true
)

