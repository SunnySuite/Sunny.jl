using Sunny
using LinearAlgebra
using MPI

# MPI rank is 0-based index
rank, N_ranks = init_MPI()

# Use kT as control parameter (α)
kT_min = 0.5
kT_max = 5.0

# Geometric kT schedule
kT = kT_min * (kT_max/kT_min) ^ ((rank)/(N_ranks-1))

# Linear kT schedule
#kT = collect(range(kT_min, kT_max, length=N_ranks))[rank+1]

#=
# Linear schedule for external field strength
H_min = -4.0
H_max =  4.0
H = collect(range(H_min, H_max, length=N_ranks))[rank+1] / Sunny.BOHR_MAGNETON

# Fixed temperature
kT = 2.5
=#

# Nearest-neighbor ferromagnetic Ising (z-dir) interactions
J = -1.0
FM = exchange(diagm([J, J, J]), Bond(1, 1, [1, 0, 0]))
interactions = [FM]
#H_ext = external_field([0, 0, H])
#interactions = [FM, H_ext]

# SC lattice -- extend lattice vector in z-direction to emulate 2D
lvecs = Sunny.lattice_vectors(1.0, 1.0, 2.0, 90, 90 ,90)
bvecs = [[0.0, 0.0, 0.0]]
crystal = Crystal(lvecs, bvecs)

# Make system and randomize spins
extent = (20, 20, 1)
system = SpinSystem(crystal, interactions, extent, [SiteInfo(1,1,1)])
randflips!(system)

# Make replica for REMC
α = 1/kT
replica = Replica(IsingSampler(system, kT, 1), α)

# Define how the setting of α's changes
#  the sampling distribution or the Hamiltonian
function set_α!(replica::Replica, α::Float64)
    replica.α = α
    set_temp!(replica.sampler, 1/α)
end

# Run feedback-optimized replica exchange
α_opt = run_FBO!(
    replica,
    set_α!;
    max_mcs = 100_000,
    update_interval = 20_000
)

#=
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
=#
