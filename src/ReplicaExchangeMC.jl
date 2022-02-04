""" 
    mutable struct Replica{AS<:AbstractSampler}

Replica type for replica exchange Monte Carlo. A Sunny sampler type must
be provided during construction.
"""
mutable struct Replica{AS <: AbstractSampler}
    rank::Int64             # MPI rank of sampler is ∈ [0, N_ranks-1]
    N_ranks::Int64          # Total number of MPI ranks (samplers, replicas) 
    N_rex::Int64            # Count of total number of attempted replica exchanges
    rex_dir::Int64          # Replica's current exchange direction ("up"=2, "down"=1)
    α::Float64              # The current value for control parameter
    nn_ranks::Vector{Int64} # MPI ranks of nearest-neighbor samplers, -1 signals no neighbor
    label::Int64            # Label for replica control parameter tracking. Starts as N_ranks+rank+1,
                            #   and updates as swaps occur. Label is < 0 when last α is most recently
                            #   visited and > 0 when first α most recently visited
    N_up::Int64             # Number of exchanges with replica most recently visiting first α
    N_down::Int64           # Number of exchanges with replica most recently visiting last α
    sampler::AS             # Sampler containing/updating this replica's system
end

"""
Initialize MPI and return process rank, number of ranks
"""
function init_MPI()
    if !MPI.Initialized()
        MPI.Init()
    end
    return (MPI.Comm_rank(MPI.COMM_WORLD), MPI.Comm_size(MPI.COMM_WORLD))
end

"""
Constructor which initializes MPI communicator and sets sampler. -> Use on
single communicator with no groups for now
"""
function Replica(sampler::AS, α::Float64) where {AS <: AbstractSampler}
    # MPI must be initialized by user after importing MPI
    if !MPI.Initialized()
        println("""MPI not initialized. Please add "MPI.Init()" at start of script""")
        exit(-1)
    end

    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    N_ranks = MPI.Comm_size(MPI.COMM_WORLD)

    Random.seed!(round(Int64, time()*1000))

    # Even rank exch. down when rex_dir==1, up when rex_dir==2
    nn_ranks = [rank-1, rank+1]

    # First and last replicas only exch. in one dir.
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

    return Replica(rank, N_ranks, 0, rex_dir, α, nn_ranks, label, 0, 0, sampler)
end

"""
Propose and accept/reject a replica exchange. -> Swap spin configurations
instead of the α's for now (no bookkeeping or master proc. needed). -> Use
deterministic even/odd scheme.
"""
function replica_exchange!(replica::Replica)
    # First and last replicas only exch. in one dir.
    if replica.nn_ranks[replica.rex_dir] < 0
        return false
    end
    replica.N_rex += 1

    # Rank of exchange partner
    rex_rank = replica.nn_ranks[replica.rex_dir]

    # Action of this rank's configuration
    S_curr = running_energy(replica.sampler) / get_temp(replica.sampler)

    # Backup current configuration
    backup_sites = deepcopy(replica.sampler.system.sites)

    # Swap trial configuration with partner
    MPI.Sendrecv!(
                        backup_sites, rex_rank, 1,
        replica.sampler.system.sites, rex_rank, 1,
        MPI.COMM_WORLD
    )

    # This rank's action with partner's configuration
    S_exch = energy(replica.sampler.system) / get_temp(replica.sampler)

    # change in action due to exchange
    ΔS = S_curr - S_exch

    # Calculate exchange probability
    rex_accept = false
    if replica.rank < rex_rank
        rex_ΔS = MPI.Recv(Float64, rex_rank, 2, MPI.COMM_WORLD)[1]
        ln_P = ΔS + rex_ΔS

        # Replica exchange acceptance criterion
        if (ln_P >= 0) || (rand() < exp(ln_P))
            rex_accept = true
        end

        MPI.Send(rex_accept, rex_rank, 3, MPI.COMM_WORLD)
    else
        MPI.Send(ΔS, rex_rank, 2, MPI.COMM_WORLD)

        rex_accept = MPI.Recv(Bool, rex_rank, 3, MPI.COMM_WORLD)[1]
    end

    # Reject exchange
    if !rex_accept
        replica.sampler.system.sites .= backup_sites
        return false
    end

    # Recalculate energy and magnetization of new state
    reset_running_energy!(replica.sampler)
    reset_running_mag!(replica.sampler)

    return true
end

"""
Exchange labels between two replicas. 

Each replica has a label that starts as (N_ranks + rank+1) until the replica
reaches an extremal temperature, after which the label ∈ [1, N_ranks]. Labels
are > 0 after most recently visiting first α and < 0 after most recently 
visiting lasst α

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

    # Change the trip direction when α extrema are reached
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
Update α schedule using feedback-optimization scheme from
Katzgraber et. al. (https://arxiv.org/pdf/cond-mat/0602085v3.pdf).

Note: use f(θ) that goes from θ_max -> θ_min instead of original paper
Note: use (k-1)/(M-1) as target for integral in Eq. 11 instead of k/M 

Integration uses dummy variable for θ = 1/α (e.g. θ ↔ kT)    
"""
function feedback_update!(replica::Replica, set_α!::F; w::Float64=0.0, dθ′::Float64=1e-5) where {F <: Function}
    # Storage for updated α
    α_new = [replica.α]

    # Replica flow
    N_labeled = (replica.N_up + replica.N_down)
    f = ((N_labeled > 0) ? replica.N_down/N_labeled : 0.0)

    # Gather θ and replica flow 'f' from all replicas
    all_θ_and_f = MPI.Gather([1/replica.α, f], 0, MPI.COMM_WORLD)

    if replica.rank == 0
        all_θ_and_f = reshape(all_θ_and_f, :, replica.N_ranks)
        θ = all_θ_and_f[1, :]
        f = all_θ_and_f[2, :]

        M = replica.N_ranks
        L = collect(range(0, 1.0, length=M))

        # Use 'smoothed replica flow' from Hamze et. al. (https://arxiv.org/pdf/1004.2840.pdf)
        fs = (1-w).*f .+ (w .* L)

        # Use crude linear approximation for df/dθ
        Δf = [ ((fs[i+1]-fs[i]) < 0) ? 0 : abs(fs[i+1]-fs[i]) for i in 1:M-1 ]
        Δθ = abs.(θ[2:end] .- θ[1:end-1])

        # Calculated normalized density for adjusted θ's
        η = sqrt.( (Δf ./ Δθ) ./ Δθ )
        C = sum(η .* Δθ)
        η ./= C

        # New θ's (and α's) have fixed end points
        θ_new = zeros(M)
        θ_new[1] = θ[1]
        θ_new[M] = θ[M]

        θ′= θ[1]
        cηf = 0.0
        pos = 2

        # find new temperatures using Eq. 11 from Katzgraber et. al.
        for i in 1:M-1
            for j in 1:Int(round(Δθ[i]/dθ′))
                cηf += η[i]*dθ′
                θ′ += dθ′

                if (pos < M) && (cηf >= L[pos])
                    θ_new[pos] = θ′
                    pos += 1
                end
            end
        end
        # Send out updated α to all replicas
        MPI.Scatter!(MPI.UBuffer(1 ./ θ_new, 1), α_new, 0, MPI.COMM_WORLD)
    else
        MPI.Scatter!(nothing, α_new, 0, MPI.COMM_WORLD)
    end

    # Assign new value of α and set in replica sampler
    # -> corresponds to update of kT or other Hamiltonian parameter
    set_α!(replica, α_new[1])

    return nothing
end

"""
"""
function run_FBO!(
    replica::Replica,
    set_α!::F;
    therm_mcs::Int64=1000,
    rex_interval::Int64=1,
    max_mcs::Int64=500_000,
    update_interval::Int64=100_000,
    w::Float64=0.0,
    dθ′::Float64=1e-5,
    print_ranks_trajectory::Vector{Int64}=Int64[],
    trajectory_interval::Int64=500_000
) where {F <: Function}

    # Initialize energy and equilibrate replica to it's distribution
    reset_running_energy!(replica.sampler)
    reset_running_mag!(replica.sampler)

    thermalize!(replica.sampler, therm_mcs)

    replica.sampler.nsweeps = 1

    rex_accepts = zeros(Int64, 2)
    N_updates = 0

    # Start feedback optimization simulation
    for mcs in 1:max_mcs
        sample!(replica.sampler)

        if (replica.rank == 0) && (mcs % 1000 == 0)
            println("FBOPT ", mcs, " mcs")
        end

        # Attempt replica exchange
        if mcs % rex_interval == 0
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
            if (replica.N_rex % trajectory_interval == 0) && (abs(replica.label)-1 in print_ranks_trajectory)
                if replica.label <= replica.N_ranks
                    fname = @sprintf("P%03d_FBO_trajectory_%03d.dat", abs(replica.label)-1, N_updates)
                    output = open(fname, "a")
                    println(output, replica.rank)
                    close(output)
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

            # Print out α, replica flow, and acceptance rates for each update
            fname = @sprintf("FBO_%03d.dat", N_updates)
            gather_print(
                replica,
                [replica.α, fs, A],
                fname
            )

            # Perform feedback optimization for control parameters
            feedback_update!(replica, set_α!; w=w, dθ′=dθ′)

            N_updates += 1
            rex_accepts[1] = rex_accepts[2] = 0
            reset_stats!(replica)
        end
    end

    rex_accepts[1] = rex_accepts[2] = 0
    reset_stats!(replica)

    # Return a vector of optimized α's
    return MPI.Allgather([replica.α], MPI.COMM_WORLD)
end

"""
"""
function reset_stats!(replica::Replica)
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
    return nothing
end

"""
Start a parallel tempering (PT) simulation that saves measured averages,
histograms, minimal energies, and configurations.

Run w/ the command "mpiexec -n [n_procs] julia --project [julia_script.jl]".

# Arguments
- `replica::Replica`: Replica type that contains MPI info and system data

- `therm_mcs::Int64`: Number of initial thermalization MC sweeps

- `measure_interval::Int64`: Number of MC sweeps between measurements

- `rex_interval::Int64`: Number of MC sweeps between replica exchange attempts

- `max_mcs::Int64`: Maximum number of MC sweeps to allow during simulation

- `bin_size::Float64`: Width in energy space for each histogram bin

- `print_hist::Bool`: Print energy histograms to file for each rank if true

- `print_xyz_ranks::Vector{Float64}`: If an MPI rank is in this array, the xyz
  coordinates for new minimal energies are printed to file

- `print_ranks_trajectory::Vector{Float64}`: If an MPI rank is in this array, the 
  timeseries for the replica starting with this rank will printed to file

- `trajectory_interval::Int64`: Number of replica exchange attempts between printed
  timeseries points

- `measure_replica_flow::Bool`: If true, then the flux of downward replica flow in 
  α space is measured and printed after the simulation

"""
function run_REMC!(
    replica::Replica; 
    therm_mcs::Int64=1000, 
    measure_interval::Int64=10, 
    rex_interval::Int64=1, 
    max_mcs::Int64=1_000_000, 
    bin_size::Float64=1.0, 
    print_hist::Bool=false, 
    print_xyz_ranks::Vector{Int64}=Int64[],
    print_ranks_trajectory::Vector{Int64}=Int64[],
    trajectory_interval::Int64=1_000_000,
    measure_replica_flow::Bool=false
)
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
    Emin = typemax(Float64)
    hist = BinnedArray{Float64, Int64}(bin_size=bin_size)
    rex_accepts = zeros(Int64, 2)
    N_measure = 0

    # Start PT with finite length
    for mcs in 1:max_mcs
        sample!(replica.sampler)

        # Crude progress printing to stdout
        if (replica.rank == 0) && (mcs % 1000 == 0)
            println("PT ", mcs, " mcs")
        end

        # Attempt replica exchange
        if mcs % rex_interval == 0 
            if replica_exchange!(replica)
                rex_accepts[replica.rex_dir] += 1

                # Exchange labels
                if measure_replica_flow
                    label_exchange!(replica, replica.nn_ranks[replica.rex_dir])
                end
            end

            # Record replica flow (for quality evaluation) 
            if measure_replica_flow && replica.nn_ranks[replica.rex_dir] != -1
                if replica.label < 0
                    replica.N_down += 1
                elseif replica.label < replica.N_ranks+1
                    replica.N_up += 1
                end
            end

            # Print out replica exchange timeseries -- inefficient IO
            if (replica.N_rex % trajectory_interval == 0) && (abs(replica.label)-1 in print_ranks_trajectory)
                if replica.label <= replica.N_ranks
                    fname = @sprintf("P%03d_trajectory.dat", abs(replica.label)-1)
                    rex_output = open(fname, "a")
                    println(rex_output, replica.rank)
                    close(rex_output)
                end
            end

            # Alternate up/down pairs of replicas for exchanges
            replica.rex_dir = 3 - replica.rex_dir
        end

        # Measurements
        if mcs % measure_interval == 0
            E = running_energy(replica.sampler)
            m = (running_mag(replica.sampler))[3]

            # New minimum energy measured
            if E < Emin
                Emin = E

                # Print out evolution of Emin w.r.t MC time
                f = open(@sprintf("P%03d_Emin.dat", replica.rank), "a")
                println(f, mcs," ",Emin)
                close(f)

                # Overwrite Emin coordinates in saved file
                if replica.rank in print_xyz_ranks
                    xyz_output = open(@sprintf("P%03d_Emin.xyz", replica.rank), "w")
                    xyz_to_file(replica.sampler.system, xyz_output)
                    close(xyz_output)
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
    kT = get_temp(replica.sampler)
    U  /= N_measure
    U2 /= N_measure
    M  /= N_measure
    M2 /= N_measure
    C = 1.0/kT^2 * (U2 - U^2)
    X = 1.0/kT   * (M2 - M^2)

    # Acceptance rate between replicas i and i+1
    A = ((replica.rank == replica.N_ranks-1) ? 0.0 : rex_accepts[2]/replica.N_rex)

    # Gather and print: rank, "down" exch. accepts, "up" exch. accepts, acceptance ratio   
    gather_print(
        replica, 
        [replica.α, A], 
        "replica_exchanges.dat"
    )

    # Gather and print thermodynamic measurements
    gather_print(
        replica, 
        [replica.α, M/system_size, X/system_size, U/system_size, C/system_size], 
        "measurements.dat"
    )

    # Write histogram to file for each process if specified
    if print_hist
        f = open(@sprintf("P%03d_hist.dat", replica.rank), "w")
        println(f, hist)
        close(f)
    end
    
    # Print replica flow 'f' to file
    if measure_replica_flow
        N_labeled = (replica.N_up + replica.N_down)    
        f = ((N_labeled > 0) ? replica.N_down/N_labeled : 0.0)

        gather_print(
            replica,
            [replica.α, f],
            "replica_flow.dat"
        )
    end

    return :SUCCESS
end

