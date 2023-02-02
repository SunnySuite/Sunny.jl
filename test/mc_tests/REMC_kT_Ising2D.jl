using Sunny
using LinearAlgebra
using MPI

# MPI rank is 0-based index
rank, N_ranks = init_MPI()

kT_min = 0.1
kT_max = 10.0

# Geometric temperature schedule - make smallest β first entry
kT = kT_min * (kT_max/kT_min) ^ ((N_ranks-1-rank)/(N_ranks-1))

# Linear temperature schedule
#kT = collect(range(kT_min, kT_max, length=N_ranks))[N_ranks-rank]

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
g = 1 / Sunny.BOHR_MAGNETON
system = SpinSystem(crystal, interactions, extent, [SpinInfo(1, S=1; g)])
for idx = all_sites(system)
    rand(sys.rng, Bool) && sys.dipoles[idx] *= -1
end
    
# Make replica for REMC
α = 1/kT
replica = Replica(IsingSampler(system, kT, 1), α)

# Define how the setting of α's changes
#  the sampling distribution or the system
function set_α!(replica::Replica, α::Float64)
    replica.α = α
    set_temp!(replica.sampler, 1/α)
end

# Run feedback-optimized replica exchange
α_opt = run_FBO!(
    replica,
    set_α!;
    max_mcs = 50_000,
    rex_interval = 1,
    update_interval = 10_000,
    w = 0.0
)

# Run Replica Exchange MC (as PT) using optimized α schedule
run_REMC!(
    replica;
    therm_mcs=5000,
    measure_interval=50,
    rex_interval=10,
    max_mcs=500_000,
    bin_size=1.0,
    print_hist=true
)
