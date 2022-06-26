using Sunny
using LinearAlgebra
using MPI

# MPI rank is 0-based index
rank, N_ranks = init_MPI()

# Linear schedule for external field strength
H_min = -4.0
H_max =  4.0
H = collect(range(H_min, H_max, length=N_ranks))[rank+1]

# Fixed temperature
kT = 1.0

function create_system(α::Float64)
    # Nearest-neighbor ferromagnetic Ising (z-dir) interactions
    J = -1.0
    FM = exchange(diagm([J, J, J]), Bond(1, 1, [1, 0, 0]))
    H_ext = external_field([0, 0, α])
    interactions = [FM, H_ext]

    # SC lattice -- extend lattice vector in z-direction to emulate 2D
    lvecs = Sunny.lattice_vectors(1.0, 1.0, 2.0, 90, 90 ,90)
    bvecs = [[0.0, 0.0, 0.0]]
    crystal = Crystal(lvecs, bvecs)

    # Make system and randomize spins
    extent = (120, 120, 1)
    g = 1 / Sunny.BOHR_MAGNETON
    system = SpinSystem(crystal, interactions, extent, [SiteInfo(1; spin_rescaling=1.0, g)])
    
    return system
end

α = H
system = create_system(α)   
randflips!(system)
    
# Make replica for REMC
replica = Replica(IsingSampler(system, kT, 1), α)

# Define how the setting of α's changes
#  the sampling distribution or the system
function set_α!(replica::Replica, α::Float64)
    replica.α = α
    copy_sites = deepcopy(replica.sampler.system.sites)
    replica.sampler.system = create_system(α)
    replica.sampler.system.sites .= copy_sites
    reset_running_energy!(replica.sampler)
    reset_running_mag!(replica.sampler)
end

# Run feedback-optimized replica exchange
α_opt = run_FBO!(
    replica,
    set_α!;
    max_mcs = 200_000,
    rex_interval = 1,
    update_interval = 10_000,
    w = 0.0,
    #print_ranks_trajectory=[0],
    #trajectory_interval=10
)

# Run Replica Exchange MC (as PT) using optimized α schedule
run_REMC!(
    replica;
    therm_mcs=5000,
    measure_interval=10,
    rex_interval=1,
    max_mcs=500_000,
    bin_size=1.0,
    print_hist=true
)
