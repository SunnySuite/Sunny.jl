using MPI

""" 
Container for sampler type as there is no deriving from concrete types in Julia
"""
mutable struct Replica
    rank::Int64
    N_ranks::Int64
    rex_dir::Int64
    nn_ranks::Vector{Int64}
    nn_βs::Vector{Float64}

    sampler::AbstractSampler
end


"""
Constructor which initializes MPI communicator and sets sampler.
-> Use on single communicator with no groups for now
"""
function Replica(s::AbstractSampler)
    # initialize MPI communicator and variables
    MPI.Init()
    MPI_COMM_WORLD = MPI.COMM_WORLD

    rank = MPI.Comm_rank(MPI_COMM_WORLD)
    N_ranks = MPI.Comm_size(MPI_COMM_WORLD)

    Random.seed!(round(Int64, time()*1000))

    # even rank exch. down when rex_dir==0, up when rex_dir==1
    nn_ranks = [
        rank + 2*(rank%2) - 1,
        rank + 1 - 2*(rank%2)
    ]
    # first and last replicas only exch. in one dir.
    for i in 1:2
        if (nn_ranks[i] < 0) || (nn_ranks[i] >= N_ranks)
            nn_ranks[i] = -1
        end
    end

    return Replica(rank, N_ranks, 1, nn_ranks, [-1,-1], s)
end


"""
Propose and accept/reject a replica exchange.
-> Swap spin configurations instead of the β's for now 
(no bookkeeping or master proc. needed).
-> Use deterministic even/odd scheme (irreversible).
"""
function replica_exchange!(replica::Replica)
    rex_ID = replica.rex_dir+1

    # first and last replicas only exch. in one dir.
    if replica.nn_ranks[rex_ID] < 0
        return false
    end

    rex_rank = replica.nn_ranks[rex_ID]
    E_curr = [running_energy(replica.sampler)]
    rex_accept = [false]

    # propose replica exch.
    if replica.rank < rex_rank
        # send energy to neighbor
        MPI.Send(E_curr, rex_rank, 1, MPI.COMM_WORLD)

        # receive accept/reject info from neighbor
        MPI.Recv!(rex_accept, rex_rank, 2, MPI.COMM_WORLD)
    else
        # receive energy from neighbor
        E_rex = [0.0]
        MPI.Recv!(E_rex, rex_rank, 1, MPI.COMM_WORLD)

        # replica exch. acceptance probability
        P_rex = exp( (replica.nn_βs[rex_ID] - replica.sampler.β)*(E_rex[1] - E_curr[1]) )

        # acceptance criterion.
        if (P_rex > 1.0) || (rand() <= P_rex)
            rex_accept[1] = true
        end

        # send accept/reject info to neighbor
        MPI.Send(rex_accept, rex_rank, 2, MPI.COMM_WORLD)
    end

    # replica exch. rejected
    if !rex_accept[1]
        return false
    end

    # accept replica exch.: swap configurations
    MPI.Sendrecv!(replica.sampler.system.sites, rex_rank, 3,
                  replica.sampler.system.sites, rex_rank, 3, MPI.COMM_WORLD)

    # recalculate energy of new state
    replica.sampler.E = energy(replica.sampler.system)

    return true
end


"""
Start a parallel tempering (PT) simulation. 
Run w/ the command "mpiexec -n [n_procs] julia --project [julia_script.jl]".
"""
function run!(replica::Replica, T_sched::Function;
                therm_mcs=1000, measure_interval=10, rex_interval=100, max_mcs=500_000, bin_size=1.0)
    # set replica sampling β using temperature (β = 1.0/(kT)) 
    T = T_sched(replica.rank+1, replica.N_ranks)
    set_temp!(replica.sampler, T)

    # set replica neighbor β's  
    for i in 1:2
        if replica.nn_ranks[i] != -1
            replica.nn_βs[i] = 1.0 / T_sched(replica.nn_ranks[i]+1, replica.N_ranks)
        end
    end

    # equilibrate replica to it's distribution
    thermalize!(replica.sampler, therm_mcs)

    N_rex = max_mcs / rex_interval
    N_measure = cld(rex_interval, measure_interval)
    replica.sampler.nsweeps = measure_interval

    A = BinnedArray{Float64, Int64}(bin_size=bin_size)
    rex_accepts = zeros(Int64, 2)

    # start PT with finite length
    for i in 1:N_rex

        # measurements between exchanges
        for j in 1:N_measure
            # decorrelation steps
            sample!(replica.sampler)

            # add measurements here 
            #...

            # histogram for WHAM (TODO)
            A[running_energy(replica.sampler)] += 1
        end

        # attempt replica exchange
        if replica_exchange!(replica)
            rex_accepts[replica.rex_dir+1] += 1
        end

        # alternate up/down pairs of replicas for exchanges
        replica.rex_dir = 1 - replica.rex_dir
    end

    # write replica exchange rates to file for each process
    f = open(@sprintf("P%03d_rex.dat", replica.rank), "w")
    println(f, "REX accepts (down, up) = (", rex_accepts[1], " ", rex_accepts[2],"), total = ", sum(rex_accepts), " / ", N_rex)
    close(f)

    # write histogram to file for each process
    f = open(@sprintf("P%03d_hist.dat", replica.rank), "w")
    println(f, A)
    close(f)

    return :SUCCESS
end

