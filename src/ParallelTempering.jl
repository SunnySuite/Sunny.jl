
""" 
Container for sampler type as there is no deriving from concrete types in Julia
"""
mutable struct Replica{S<:AbstractSampler}
    rank     :: Int64
    N_ranks  :: Int64
    rex_dir  :: Int64
    nn_ranks :: Vector{Int64}
    nn_βs    :: Vector{Float64}
    sampler  :: S
end

"""
Constructor which initializes MPI communicator and sets sampler.
-> Use on single communicator with no groups for now
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

    # First and last replicas only exchange in one dir.
    for i in 1:2
        if (nn_ranks[i] < 0) || (nn_ranks[i] >= N_ranks)
            nn_ranks[i] = -1
        end
    end

    # Start even ranks exchanging down, odd exchanging up
    rex_dir = (rank % 2 == 0) ? 1 : 2

    return Replica{S}(rank, N_ranks, rex_dir, nn_ranks, [-1,-1], sampler)
end


"""
Propose and accept/reject a replica exchange.
-> Swap spin configurations instead of the β's for now 
(no bookkeeping or master proc. needed).
-> Use deterministic even/odd scheme (irreversible).
"""
function replica_exchange!(replica::Replica)
    # First and last replicas only exchange in one dir.
    rex_rank = replica.nn_ranks[replica.rex_dir]
    if rex_rank < 0
        return false
    end
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

        # Replica exchange acceptance probability
        P_rex = exp( (replica.nn_βs[replica.rex_dir] - replica.sampler.β)*(E_rex[1] - E_curr[1]) )

        # Acceptance criterion
        if (P_rex > 1.0) || (rand() <= P_rex)
            rex_accept[1] = true
        end

        # Send accept/reject info to neighbor
        MPI.Send(rex_accept, rex_rank, 2, MPI.COMM_WORLD)
    end

    # Replica exchange rejected
    if !rex_accept[1]
        return false
    end

    # Accept replica exchange: swap configurations
    MPI.Sendrecv!(replica.sampler.system.sites, rex_rank, 3,
                  replica.sampler.system.sites, rex_rank, 3, MPI.COMM_WORLD)

    # Recalculate energy and magnetization of new state
    replica.sampler.E = energy(replica.sampler.system)
    replica.sampler.M = sum(replica.sampler.system)

    return true
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
        data_all = reshape(data_all, replica.N_ranks, :)

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
    run_parallel_temp!(replica::Replica, kt_sched<:Function;
                       therm_mcs=1000, measure_interval=10, rex_interval=100,
                       max_mcs=500_000, bin_size=1.0, print_hist=false)

Start a parallel tempering (PT) simulation. 
Run scripts containing this function w/ the command
    "mpiexec -n [n_procs] julia --project [julia_script.jl]".
"""
function run_parallel_temp!(
    replica::Replica, kT_sched<:Function;
    therm_mcs=1000, measure_interval=10, rex_interval=100,
    max_mcs=500_000, bin_size=1.0, print_hist=false
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

    # Set sampler's initial energy and magnetization values
    reset_running_energy!(replica.sampler)
    reset_running_mag!(replica.sampler)

    # Thermalize replica to its equilibrium distribution
    thermalize!(replica.sampler, therm_mcs)

    system_size = length(replica.sampler.system)

    replica.sampler.nsweeps = 1

    # Manually include these basic measurements for now
    M = zero(Vec3)
    U = 0.0
    U2 = 0.0
    C = 0.0

    N_rex = 0
    N_measure = 0

    hist = BinnedArray{Float64, Int64}(bin_size=bin_size)
    rex_accepts = zeros(Int64, 2)

    # Start PT with finite length
    for mcs in 1:max_mcs
        sample!(replica.sampler)

        # Replica exchanges
        if mcs % rex_interval == 0 
            # Attempt replica exchange
            if replica_exchange!(replica)
                rex_accepts[replica.rex_dir] += 1
            end

            N_rex += 1

            # Alternate up/down pairs of replicas for exchanges
            replica.rex_dir = 3 - replica.rex_dir
        end

        # Measurements
        if mcs % measure_interval == 0
            E = running_energy(replica.sampler)
            m = running_mag(replica.sampler)

            U += E
            U2 += E*E
            M += m

            N_measure += 1

            # histogram for energies
            hist[E] += 1
        end
    end

    # Normalize averages
    U  /= N_measure
    U2 /= N_measure
    M  /= N_measure
    C = 1.0 / kT^2 * (U2 - U*U)
 
    # Gather and print: temperature, mag. per spin, avg. energy per spin, specific heat   
    gather_print(
        replica, 
        [T, norm(M)/system_size, U/system_size, C/system_size], 
        "measurements.dat"
    )

    # Gather and print: rank, "down" exch. accepts, "up" exch. accepts, acceptance ratio   
    gather_print(
        replica, 
        [replica.rank, rex_accepts[1], rex_accepts[2], sum(rex_accepts)/N_rex], 
        "replica_exchanges.dat"
    )

    # Write histogram to file for each process if specified
    if print_hist
        f = open(@sprintf("P%03d_hist.dat", replica.rank), "w")
        println(f, hist)
        close(f)
    end

    return :SUCCESS
end

