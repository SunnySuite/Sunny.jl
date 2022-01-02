""" 
    mutable struct Replica{S <: AbstractSampler}

Replica type for parallel tempering (PT) Monte Carlo. A Sunny sampler type must
be provided during construction.
"""
mutable struct Replica{S <: AbstractSampler}
    rank     :: Int64           # MPI rank of sampler ∈ [0, N_ranks-1]
    N_ranks  :: Int64           # Total number of MPI ranks (samplers, replicas)
    N_rex    :: Int64           # Count of total number of attempted replica exchanges
    rex_dir  :: Int64           # Replica's current exchange direction ("up"=2, "down"=1)
    nn_ranks :: Vector{Int64}   # MPI ranks of nearest-neighbor samples, -1 signals no neighbor
    nn_βs    :: Vector{Float64} # Inverse temperatures of nearest-neighbor ranks
    label    :: Int64           # Label for replica temperature tracking. Starts as N_ranks+rank+1,
                                #   and updates as swaps occur. Label is < 0 when higher T is last
                                #   visited and > 0 when lowest temperature last visited
    N_up     :: Int64           # Number of exchanges with replica which last visited lowest T
    N_down   :: Int64           # Number of exchanges with replica which last visited highest T
    sampler  :: S               # Sampler containing/updating this replica's system
end

"""
Constructor which initializes MPI communicator and sets sampler. -> Use on
single communicator with no groups for now
"""
function Replica(sampler::S) where {S <: AbstractSampler}
    # Initialize MPI communicator and variables
    MPI.Init()
    MPI_COMM_WORLD = MPI.COMM_WORLD

    rank = MPI.Comm_rank(MPI_COMM_WORLD)
    N_ranks = MPI.Comm_size(MPI_COMM_WORLD)

    Random.seed!(round(Int64, time()*1000))

    # Each rank exchanges down when rex_dir==1, up when rex_dir==2
    nn_ranks = [rank-1, rank+1]

    # First and last replicas only exchange in one dir
    for i in 1:2
        if (nn_ranks[i] < 0) || (nn_ranks[i] > N_ranks-1)
            nn_ranks[i] = -1
        end
    end

    # Start extremal ranks with up/down trip direction
    if rank == 0
        label = 1
    elseif rank == N_ranks-1
        label = -N_ranks
    else
        label = N_ranks + rank+1
    end

    # Start even ranks exchanging down
    rex_dir = (rank % 2 == 0) ? 1 : 2

    return Replica(rank, N_ranks, 0, rex_dir, nn_ranks, [-1.0,-1.0], label, 0, 0, sampler)
end

"""
Propose and accept/reject a replica exchange. -> Swap spin configurations
instead of the β's for now (no bookkeeping or master proc. needed). -> Use
deterministic even/odd scheme (irreversible).
"""
function replica_exchange!(replica::Replica)
    # First and last replicas only exchange in one dir
    rex_rank = replica.nn_ranks[replica.rex_dir]
    if rex_rank < 0
        return false
    end
    replica.N_rex += 1

    rex_rank = replica.nn_ranks[replica.rex_dir]
    E_curr = [running_energy(replica.sampler)]
    rex_accept = [false]

    # Propose replica exchange
    if replica.rank < rex_rank
        # Send energy to neighbor
        MPI.Send(E_curr, rex_rank, 1, MPI.COMM_WORLD)

        # Receive accept/reject info from neighbor
        MPI.Recv!(rex_accept, rex_rank, 2, MPI.COMM_WORLD)
    else
        # Receive energy from neighbor
        E_rex = [0.0]
        MPI.Recv!(E_rex, rex_rank, 1, MPI.COMM_WORLD)

        # Replica exch. acceptance probability
        ln_P_rex = (replica.nn_βs[replica.rex_dir] - 1 / get_temp(replica.sampler)) * (E_rex[1] - E_curr[1])

        # Acceptance criterion
        if (ln_P_rex >= 0.0) || (rand() <= exp(ln_P_rex))
            rex_accept[1] = true
        end

        # Send accept/reject info to neighbor
        MPI.Send(rex_accept, rex_rank, 2, MPI.COMM_WORLD)
    end

    # Replica exchange rejected
    if !rex_accept[1]
        return false
    end

    # Accept replica exch.: swap configurations
    MPI.Sendrecv!(        
        deepcopy(replica.sampler.system.sites), rex_rank, 3,
                 replica.sampler.system.sites , rex_rank, 3, 
        MPI.COMM_WORLD
    )
    # Recalculate energy and magnetization of new state
    reset_running_energy!(replica.sampler)
    reset_running_mag!(replica.sampler)

    return true
end

"""
Exchange labels between two replicas. 

Each replica has a label that starts as (N_ranks + rank+1) until the replica
reaches an extremal temperature, after which the label ∈ [1, N_ranks]. Labels
are > 0 after last visiting T_min and < 0 after last visiting T_max.

The label trip direction (±) is used to calculate replica flow for FBOPT and
used to track all replica trajectories.

TODO: use to calculate average round-trip time for replicas
"""
function label_exchange!(replica::Replica, rex_rank::Int64)
    # Exchange labels
    rex_label = [0]
    MPI.Sendrecv!(
        [replica.label], rex_rank, 4,
        rex_label,       rex_rank, 4, 
        MPI.COMM_WORLD
    )
    replica.label = rex_label[1]

    # Set label of replica as initialized (has trip direction)
    if (replica.rank == 0) || (replica.rank == replica.N_ranks-1)
        if replica.label > replica.N_ranks
            replica.label -= replica.N_ranks
        end
    end

    # Change the trip direction when extrema are reached
    if (replica.rank == 0) && (replica.label < 0)
        replica.label *= -1
    elseif (replica.rank == replica.N_ranks-1) && (replica.label > 0)
        replica.label *= -1
    end

    return nothing
end

"""
Gather data from each process and print to specified filename
-> only handle 1D data buffer right now
"""
function gather_print(replica::Replica, data, fname::String)
    # Have root (rank 0) gather data from all processes
    data_all = MPI.Gather(data, 0, MPI.COMM_WORLD)

    # Write gathered data to file w/ column formatting
    if replica.rank == 0
        data_all = reshape(data_all, :, replica.N_ranks)

        f = open(fname, "w")
        for data_c in eachcol(data_all)
            for val in data_c 
                print(f, val, "\t")
            end
            print(f, "\n")
        end
        close(f)
    end

    return nothing
end

"""
Print xyz formatted (Lx, Ly, Lz, Sx, Sy, Sz) configurations to file 
"""
function xyz_to_file(sys::SpinSystem, output::IOStream)
    sites = reinterpret(reshape, Float64, sys.lattice)
    spins = reinterpret(reshape, Float64, sys.sites)
    xyz = vcat(sites, spins)
    xyz = reshape(xyz, 6, :)

    for c in eachcol(xyz)
        @printf(output, "%f\t%f\t%f\t%f\t%f\t%f\n", c...)
    end
    println(output, "")

    return nothing
end

""" 
Update temperature schedule using feedback-optimized parallel tempering from
Katzgraber et. al. (https://arxiv.org/pdf/cond-mat/0602085v3.pdf).

Note: use f(T) that goes from T_max -> T_min instead of original paper
Note: use (k-1)/(M-1) as target for integral in Eq. 11 instead of k/M 
"""
function feedback_update!(replica::Replica, N_updates::Int64; w::Float64=0.0, dkT′::Float64=1e-4)

    kT_new = zeros(replica.N_ranks)

    N_labeled = (replica.N_up + replica.N_down)    
    f = ((N_labeled > 0) ? replica.N_down/N_labeled : 0.0)
    
    # Gather temperatures and replica flow 'f' from all replicas
    all_kT_and_f = MPI.Gather([get_temp(replica.sampler), f], 0, MPI.COMM_WORLD)

    if replica.rank == 0

        all_kT_and_f = reshape(all_kT_and_f, :, replica.N_ranks)
        kT = all_kT_and_f[1, :]
        f  = all_kT_and_f[2, :]

        M = replica.N_ranks
        L = collect(range(0, 1.0, length=M))

        # Use 'smoothed replica flow' from Hamze et. al. (https://arxiv.org/pdf/1004.2840.pdf)
        fs = (1-w).*f .+ (w .* L)

        # Use crude linear approximation for df/dT
        Δf = fs[2:end] .- fs[1:end-1]
        Δf[Δf .< 0] .= 0

        ΔkT = kT[2:end] .- kT[1:end-1]

        # Calculated normalized density for adjusted temperatures
        η = sqrt.( (Δf ./ ΔkT) ./ ΔkT )
        C = sum(η .* ΔkT)
        η ./= C

        # New temperatures have fixed end points
        kT_new[1] = kT[1]
        kT_new[M] = kT[M]

        kT′= kT[1]
        cηf = 0.0
        pos = 2

        # Find new temperatures using Eq. 11 from Katzgraber et. al.
        for i in 1:M-1
            for j in 1:floor(ΔkT[i]/dkT′)
                cηf += η[i]*dkT′
                kT′ += dkT′

                if cηf >= L[pos]
                    kT_new[pos] = kT′
                    pos += 1
                end
            end
        end
    end

    # Send new temperatures to all replicas
    MPI.Bcast!(kT_new, 0, MPI.COMM_WORLD)

    # Set replica sampler to new temperature
    set_temp!(replica.sampler, kT_new[replica.rank+1])

    # Set replica neighbor β's  
    for i in 1:2
        if replica.nn_ranks[i] != -1
            replica.nn_βs[i] = 1.0 / kT_new[replica.nn_ranks[i]+1]
        end
    end

    return nothing
end

"""
Run a feedback-optimized parallel tempering simulation and return an optimized
temperature set. No measurements are made, but this is used to diagnose replica
flow in PT simulations. 

# Arguments
- `replica::Replica`: Replica type that contains MPI info and system data

- `kT_sched::Function`: Function that takes index i::Int64 and number of
  temperatures N::Int64 and returns temperature kT::Float64

- `max_mcs_opt::Int64`: Maximum number of MC sweeps to allow during optimization

- `update_interval::Int64`: Number of MC sweeps between feedback updates

- `therm_mcs::Int64`: Number of initial thermalization MC sweeps

- `rex_interval::Int64`: Number of MC sweeps between replica exchange attempts

- `w::Float64`: Damping factor (w ∈ [0,1]) for feedback updates: 0 gives
  aggressive updates, 1 gives leaves set of kT unchanged

- `dkT′::Float64`: Step size for integration during feedback update

- `print_ranks::Vector{Float64}`: If an MPI rank is in this array, timeseries
  for the replica starting with this rank will printed to file

- `print_interval::Int64`: Number of replica exchange attempts between printed
  timeseries points

"""
function run_FBOPT!(
    replica::Replica, 
    kT_sched::Function;
    max_mcs_opt::Int64=1_000_000, 
    update_interval::Int64=200_000, 
    therm_mcs::Int64=1000, 
    rex_interval::Int64=1, 
    w::Float64=0.0, 
    dkT′::Float64=1e-4, 
    print_ranks::Vector{Int64}=Int64[], 
    print_interval::Int64=1
)
    # Set replica sampling β using temperature (β = 1.0/(kT)) 
    kT = kT_sched(replica.rank+1, replica.N_ranks)
    set_temp!(replica.sampler, kT)

    # Set replica neighbor β's  
    for i in 1:2
        if replica.nn_ranks[i] != -1
            replica.nn_βs[i] = 1.0 / kT_sched(replica.nn_ranks[i]+1, replica.N_ranks)
        end
    end

    # Initialize energy and equilibrate replica to it's distribution
    reset_running_energy!(replica.sampler)
    reset_running_mag!(replica.sampler)

    thermalize!(replica.sampler, therm_mcs)
    
    replica.sampler.nsweeps = 1

    # Start parallel tempering for optimization run
    rex_accepts = zeros(Int64, 2)
    N_updates = 0

    for mcs in 1:max_mcs_opt
        sample!(replica.sampler)

        if (replica.rank == 0) && (mcs % 1000 == 0)
            println("FBOPT ", mcs, " mcs")
        end

        if mcs % rex_interval == 0 
            # Attempt replica exchange
            if replica_exchange!(replica)
                rex_accepts[replica.rex_dir] += 1

                # Exchange labels
                label_exchange!(replica, replica.nn_ranks[replica.rex_dir])
            end

            # Record replica flow  
            if replica.nn_ranks[replica.rex_dir] != -1
                if replica.label < 0
                    replica.N_down += 1
                elseif replica.label < replica.N_ranks+1
                    replica.N_up += 1
                end
            end

            # Print out replica exchange timeseries -- inefficient IO
            if (replica.N_rex % print_interval == 0) && (abs(replica.label)-1 in print_ranks)
                if replica.label <= replica.N_ranks
                    fname = @sprintf("P%03d_trajectory_FBOPT_%03d.dat", abs(replica.label)-1, N_updates)
                    rex_output = open(fname, "a")
                    println(rex_output, replica.rank)
                    close(rex_output)
                end
            end

            # Alternate up/down pairs of replicas for exchanges
            replica.rex_dir = 3 - replica.rex_dir
        end

        if mcs % update_interval == 0
            # Replica flow 'f' and smoothed replica flow 'fs'
            N_labeled = (replica.N_up + replica.N_down)    
            f = ((N_labeled > 0) ? replica.N_down/N_labeled : 0.0)
            fs = (1-w)*f + w*replica.rank/(replica.N_ranks-1)

            # Acceptance rate between replicas i and i+1
            A = ((replica.rank == replica.N_ranks-1) ? 0.0 : rex_accepts[2]/replica.N_rex)

            # Print out temperatures, replica flow, and acceptance rates for each update
            fname = @sprintf("FBOPT_%03d.dat", N_updates)
            gather_print(
                replica,
                [get_temp(replica.sampler), fs, A],
                fname
            )

            # Perform feedback optimization for temperatures
            feedback_update!(replica, N_updates; w=w, dkT′=dkT′)

            N_updates += 1
            rex_accepts[1] = rex_accepts[2] = 0
            replica.N_up = replica.N_down = 0
            replica.N_rex  = 0

            # Reset label
            if replica.rank == 0
                replica.label = 1
            elseif replica.rank == replica.N_ranks-1
                replica.label = -replica.N_ranks
            else
                replica.label = replica.N_ranks + replica.rank+1
            end
        end
    end

    return MPI.Allgather([get_temp(replica.sampler)], MPI.COMM_WORLD)
end

"""
Start a parallel tempering (PT) simulation that saves measured averages,
histograms, minimal energies, and configurations.

Run w/ the command "mpiexec -n [n_procs] julia --project [julia_script.jl]".

# Arguments
- `replica::Replica`: Replica type that contains MPI info and system data

- `kT_sched::Function`: Function that takes index i::Int64 and number of
  temperatures N::Int64 and returns temperature kT::Float64

- `therm_mcs::Int64`: Number of initial thermalization MC sweeps

- `measure_interval::Int64`: Number of MC sweeps between measurements

- `rex_interval::Int64`: Number of MC sweeps between replica exchange attempts

- `max_mcs::Int64`: Maximum number of MC sweeps to allow during simulation

- `bin_size::Float64`: Width in energy space for each histogram bin

- `print_hist::Bool`: Print energy histograms to file for each rank if true

- `print_xyz_ranks::Vector{Float64}`: If an MPI rank is in this array, the xyz
  coordinates for new minimal energies are printed to file

"""
function run_PT!(
    replica::Replica, 
    kT_sched::Function; 
    therm_mcs::Int64=1000, 
    measure_interval::Int64=10, 
    rex_interval::Int64=1, 
    max_mcs::Int64=1_000_000, 
    bin_size::Float64=1.0, 
    print_hist::Bool=false, 
    print_xyz_ranks::Vector{Int64}=Int64[]
)
    # Set replica sampling β using temperature (β = 1.0/(kT)) 
    kT = kT_sched(replica.rank+1, replica.N_ranks)
    set_temp!(replica.sampler, kT)

    # Set replica neighbor β's  
    for i in 1:2
        if replica.nn_ranks[i] != -1
            replica.nn_βs[i] = 1.0 / kT_sched(replica.nn_ranks[i]+1, replica.N_ranks)
        end
    end

    # Initialize energy and equilibrate replica to it's distribution
    reset_running_energy!(replica.sampler)
    reset_running_mag!(replica.sampler)

    thermalize!(replica.sampler, therm_mcs)
    
    replica.sampler.nsweeps = 1

    # Manually include these basic measurements for now
    M  = 0.0
    M2 = 0.0
    X  = 0.0
    U  = 0.0
    U2 = 0.0
    C  = 0.0

    system_size = length(replica.sampler.system)

    # Record timeseries for minimum energies 
    Emin = typemax(Float64)
    mcs_Emin = Vector{Tuple{Int64, Float64}}()
    if replica.rank in print_xyz_ranks
        xyz_output = open(@sprintf("P%03d_Emin.xyz", replica.rank), "w")
    end

    hist = BinnedArray{Float64, Int64}(bin_size=bin_size)
    rex_accepts = zeros(Int64, 2)
    N_measure = 0

    # Start PT with finite length
    for mcs in 1:max_mcs
        sample!(replica.sampler)

        if (replica.rank == 0) && (mcs % 1000 == 0)
            println("PT ", mcs, " mcs")
        end

        # Attempt replica exchange
        if mcs % rex_interval == 0 
            if replica_exchange!(replica)
                rex_accepts[replica.rex_dir] += 1
            end

            # Alternate up/down pairs of replicas for exchanges
            replica.rex_dir = 3 - replica.rex_dir
        end

        # Measurements
        if mcs % measure_interval == 0
            E = running_energy(replica.sampler)
            m = norm(running_mag(replica.sampler))

            # Record minimum recorded energies
            if E < Emin
                push!(mcs_Emin, (mcs, E))
                Emin = E
                if replica.rank in print_xyz_ranks
                    xyz_to_file(replica.sampler.system, xyz_output)
                end
            end

            U  += E
            U2 += E^2
            M  += m
            M2 += m^2
            N_measure += 1

            # Histogram for energies
            hist[E] += 1
        end
    end

    # Normalize averages
    U  /= N_measure
    U2 /= N_measure
    M  /= N_measure
    M2 /= N_measure
    C = 1.0/kT^2 * (U2 - U^2)
    X = 1.0/kT   * (M2 - M^2)

    # Acceptance rate between replicas i and i+1
    A = ((replica.rank == replica.N_ranks-1) ? 0.0 : rex_accepts[2]/replica.N_rex)

    # Gather and print: rank, "down" exch. accepts, "up" exch. accepts, acceptance ratio   
    fname = @sprintf("replica_exchanges.dat")
    gather_print(
        replica, 
        [get_temp(replica.sampler), A], 
        fname
    )

    # Gather and print thermodynamic measurements
    gather_print(
        replica, 
        [get_temp(replica.sampler), M/system_size, X/system_size, U/system_size, C/system_size], 
        "measurements.dat"
    )

    # Write histogram to file for each process if specified
    if print_hist
        f = open(@sprintf("P%03d_hist.dat", replica.rank), "w")
        println(f, hist)
        close(f)
    end

    # Write Emin timeseries to file for each process
    f = open(@sprintf("P%03d_Emin.dat", replica.rank), "w")
    for r in mcs_Emin
        @printf(f, "%f\t%f\n", r...)
    end
    close(f)

    if replica.rank in print_xyz_ranks
        close(xyz_output)
    end

    return :SUCCESS
end

function run_PT!(replica::Replica, kT_sched::Vector{Float64}; kwargs...)
    kT_func = (i, N) -> kT_sched[i]
    return run_PT!(replica, kT_func; kwargs...)
end
