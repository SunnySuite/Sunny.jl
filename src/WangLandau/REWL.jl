""" 
    mutable struct WLReplica

Contains data for a Wang-Ladau 'replica' used in 
a 1D replica-exchange simulation
"""
mutable struct WLReplica
    # Spin system
    sys::SpinSystem

    # Running energy and magnetization
    E::Float64

    # Per-spin normalization
    norm::Int64

    # Binned 2D histogram
    hist::BinnedArrayND{Float64, Int64}

    # Natural log of 2D binned density of states
    ln_g::BinnedArrayND{Float64, Float64}

    win_bounds::Vector{Float64}

    # Modification factor for accumulating ln_g
    ln_f::Float64
    
    # Replica MPI rank
    rank::Int64

    # MPI communicator for unifying data within windows
    intra_comm_rank::Int64
    intra_win_comm::MPI.Comm

    # MPI communicator for replica exchanges between windows
    inter_comm_ranks::Vector{Int64}
    inter_win_comms::Vector{MPI.Comm}

    # Number of replica exchanges attempted in each dim, dir
    N_rex::Vector{Int64}

    # Alternating direction of exchanges within a dimension
    rex_dir::Int64

    # stats used for round-trip times for bins within window
    lims_curr::Vector{Float64}  # min, max energies ever visited  
    trip_dir::Int64             # (0=unassigned), (-1=down trip), (+1=up trip)
    t_dir_curr::Vector{Int64}   # current mc steps for each trip dir
    t_dir_avg::Vector{Int64}    # average mc steps for each trip dir
    N_dir::Vector{Int64}        # number of trips completed in each dir
end

# Constructor for WL replica that's called internally in 'run_REWL!()'
function WLReplica(
    sys::SpinSystem, 
    win_bounds::Vector{Float64},
    bin_size::Float64,
    rank::Int64,
    N_wins::Int64,
    rex_dir::Int64,
    per_spin::Bool
)
    norm = (per_spin ? length(sys) : 1)

    return WLReplica(
        sys,
        energy(sys)/norm, 
        norm,
        BinnedArrayND{Float64,   Int64}([win_bounds]; bin_sizes=[bin_size]),        
        BinnedArrayND{Float64, Float64}([win_bounds]; bin_sizes=[bin_size]),
        win_bounds,
        1.0,
        rank,
        -1,
        MPI.COMM_NULL,
        -ones(Int64, 2),
        [MPI.COMM_NULL for _ in 1:2],
        zeros(Int64, 2),
        rex_dir,
        Float64[Inf, -Inf],
        0,
        zeros(Int64, 2),
        zeros(Int64, 2),
        zeros(Int64, 2)
    )
end

# setup MPI comms for inter-window exchanges and intra-window communication
function init_comms!(replica::WLReplica, ranks_per_win::Int64, N_wins::Int64)
    N = ranks_per_win * N_wins

    # Ranks from COMM_WORLD that are used to create groups - not new comm. ranks
    intra_ranks = [ collect(i:i+ranks_per_win-1) for i in 0:ranks_per_win:N-1 ]

    w = div(replica.rank, ranks_per_win)
    group_world = MPI.Comm_group(MPI.COMM_WORLD)

    # Create intra-window communicator
    g = MPI.Group_incl(group_world, Int32.(intra_ranks[w+1]))
    replica.intra_win_comm = MPI.Comm_create(MPI.COMM_WORLD, g)
    replica.intra_comm_rank = MPI.Comm_rank(replica.intra_win_comm)

    # Create inter-window communicators
    inter_comms = MPI.Comm[]
    for i in 1:N_wins-1
        comm_ranks = Int32.(vcat(intra_ranks[i], intra_ranks[i+1]))
        g = MPI.Group_incl(group_world, comm_ranks)
        push!(inter_comms, MPI.Comm_create(MPI.COMM_WORLD, g))
    end

    for j in [-1, 0]
        i = j + 2
        if (w+1+j > 0) && (w+1+j < N_wins)
            replica.inter_win_comms[i] = inter_comms[w+1 + j]
            replica.inter_comm_ranks[i] = MPI.Comm_rank(inter_comms[w+1 + j])
        end
    end
    
    return nothing
end

""" 
Check histogram to determine when to advance WL iteration. 
check_type==1 uses the average flatness criterion while 
check_type==2 uses the min histogram > 1/√lnf
"""
function check_hist(replica::WLReplica; p::Float64=0.6, check_type::Int64=1)
    flat = 1

    if check_type == 1 
        # calculate average of visited bins
        avg = 0.0
        vacancies = 0
        for i in 1:replica.hist.size
            if replica.hist.visited[i]
                avg += replica.hist.vals[i]
            else
                vacancies += 1
            end
        end
        avg /= (replica.hist.size - vacancies)

        # check flatness -- all replica hists in window have to be flat
        for i in 1:replica.hist.size
            if replica.hist.visited[i] && replica.hist.vals[i] < p*avg
                flat = 0
                break
            end
        end
    elseif check_type == 2
        Hmin = 1.0 / sqrt(replica.ln_f)

        for i in 1:replica.hist.size
            if replica.hist.visited[i] && replica.hist.vals[i] < Hmin
                flat = 0
                break
            end
        end
    end

    all_flat = MPI.Allreduce(flat, MPI.PROD, replica.intra_win_comm)
    return ((all_flat > 0) ? true : false)
end

""" 
For new bins, shift ln_g to minimum existing and reset histogram
"""
function add_new!(replica::WLReplica, key::Float64)
    ln_g_min = Inf
    for i in 1:replica.hist.size
        if replica.ln_g.visited[i]
            # find min of ln_g
            if replica.ln_g.vals[i] < ln_g_min
                ln_g_min = replica.ln_g.vals[i]
            end
            # reset histogram
            replica.hist.vals[i] = 0
        end
    end
    # shift new ln_g to min value and add to histogram
    # these will be updated again after acceptance
    replica.ln_g[key] = ln_g_min - replica.ln_f
    replica.hist[key] = 0

    return nothing
end

"""
Check whether a key 'x' is within the area described by 'bounds'
"""
function bounds_check(x::Float64, bounds::Vector{Float64})
    r = round(x, digits=10)
    if (r < bounds[1]) || (r > bounds[2])
        return false
    end

    return true
end

"""
Choose random exchange partner from neighboring window.
Idea for this function from Guangjie Shi.
"""
function get_swap_partner(replica::WLReplica)
    swap_partner = [-1]
    ranks_per_win = MPI.Comm_size(replica.intra_win_comm)
    pairs = -ones(Int64, 2*ranks_per_win)

    comm_ID = replica.rex_dir

    if replica.inter_comm_ranks[comm_ID] == 0
        lim = ranks_per_win
        libre = ranks_per_win .+ collect(1:lim)

        for i in 1:ranks_per_win
            select = rand(1:lim)
            pairs[i] = libre[select]-1
            pairs[libre[select]] = i-1
            lim -= 1

            for j in select:lim
                libre[j] = libre[j+1]
            end
        end

        MPI.Scatter!(MPI.UBuffer(pairs,1), swap_partner, 0, replica.inter_win_comms[comm_ID])
    else
        MPI.Scatter!(nothing, swap_partner, 0, replica.inter_win_comms[comm_ID])
    end    

    return swap_partner[1]
end

"""
Perform a replica exchange between niehgboring windows and accept/reject
according to the density of states.
"""
function replica_exchange!(replica::WLReplica)
    comm_ID = replica.rex_dir 

    if replica.inter_comm_ranks[comm_ID] == -1
        return false
    end
    replica.N_rex[comm_ID] += 1
    
    rex_rank = get_swap_partner(replica)
    E_exch = [0.0]

    # Exchange energies to see if in bounds
    MPI.Sendrecv!([replica.E], rex_rank, 1, E_exch, rex_rank, 1, replica.inter_win_comms[comm_ID])

    # A ln_g difference of -Inf signals out of bounds
    Δln_g = [-Inf]
    if bounds_check(E_exch[1], replica.win_bounds)
        Δln_g[1] = replica.ln_g[replica.E] - replica.ln_g[E_exch[1]]
    end

    rex_accept = [false]

    if replica.inter_comm_ranks[comm_ID] < rex_rank
        # Send ln_g difference to neighbor
        MPI.Send(Δln_g, rex_rank, 2, replica.inter_win_comms[comm_ID])

        # Receive accept/reject info from neighbor
        MPI.Recv!(rex_accept, rex_rank, 3, replica.inter_win_comms[comm_ID])
    else
        # Receive lng_g difference from neighbor
        Δln_g_exch = [-Inf]
        MPI.Recv!(Δln_g_exch, rex_rank, 2, replica.inter_win_comms[comm_ID])

        # Acceptance criterion
        if !isinf(Δln_g[1]) && !isinf(Δln_g_exch[1])
            ln_p = Δln_g[1] + Δln_g_exch[1]

            if (ln_p >= 0.0) || (rand() <= exp(ln_p))
                rex_accept[1] = true
            end
        end
        # Send accept/reject info to neighbor
        MPI.Send(rex_accept, rex_rank, 3, replica.inter_win_comms[comm_ID])
    end

    # Replica exchange rejected
    if !rex_accept[1]
        return false
    end

    # Accept replica exch.: swap configurations
    MPI.Sendrecv!(
        deepcopy(replica.sys._dipoles), rex_rank, 4,
                 replica.sys._dipoles , rex_rank, 4,
        replica.inter_win_comms[comm_ID]
    )
    # Recalculate energy of new state
    replica.E = E_exch[1]
    return true
end

""" 
Initialize system to bounded range of states using throw-away WL sampling run.

# Arguments
-`replica::WLReplica`: A WLReplica type to initialize

-`max_mcs::Int64`: The max number of MC sweeps to allow

-`mc_move_type::String`: "flip" for spin flip (use with IsingSampler); "gaussian" for Gaussian perturbation of spin; 
"spherical_cap" for uniform spin displacemnt within spherical cap

-`mc_step_size::Float64`: Set the update magnitude for gaussian or spherical cap spin updates

-`limit_pad::Float64`: Energy distance for 'padding' against lower/upper limit during initialization
"""
function init_REWL!(
    replica::WLReplica;
    max_mcs::Int64 = 100_000,
    mc_move_type::String = "gaussian",
    mc_step_size::Float64 = 0.1,
    limit_pad::Float64 = 0.0
)
    init_output = open(@sprintf("init%06d.dat", replica.rank), "w")

    bounds = replica.hist.bounds[1]

    println(init_output, "# begin init")
    println(init_output, "# bounds = ", replica.hist.bounds, "\n")

    replica.E = energy(replica.sys) / replica.norm

    lim = replica.E
    fac = 0

    # Make bounds for temp. ln_g that include current state 
    init_space = deepcopy(bounds)
    for i in 1:2
        k = 2*i - 3
        if (k*lim) > (k*bounds[i])
            init_space[i] = lim + k*limit_pad
            fac = k
        end
    end

    # WL is already in bounds
    if fac == 0
        println(init_output, "\n# finish init with E = ", replica.E)
        close(init_output)
        return :SUCCESS
    end

    ln_g_tmp = BinnedArrayND{Float64, Int64}([init_space]; bin_sizes=replica.hist.bin_sizes)
    ln_g_tmp[replica.E] = 1

    # start init with finite length
    for mcs in 1:max_mcs
        for pos in CartesianIndices(replica.sys._dipoles)
            # propose single spin move
            if mc_move_type == "flip"
                new_spin = -replica.sys._dipoles[pos]
            elseif mc_move_type == "gaussian"
                new_spin = gaussian_spin_update(replica.sys._dipoles[pos], mc_step_size)
            elseif mc_move_type == "spherical_cap"
                new_spin = spherical_cap_update(replica.sys._dipoles[pos], mc_step_size)
            end

            # Calculate observables
            E_next = replica.E + local_energy_change(replica.sys, pos, new_spin) / replica.norm

            if bounds_check(E_next, init_space)
 
                if (E_next < bounds[1]) && (E_next-limit_pad > init_space[1])
                    init_space[1] = E_next - limit_pad
                end
                if (E_next > bounds[2]) && (E_next+limit_pad < init_space[2])
                    init_space[2] = E_next + limit_pad
                end

                Δln_g = ln_g_tmp[replica.E] - ln_g_tmp[E_next]

                if (Δln_g >= 0) || ( rand() <= exp(Δln_g) )
                    replica.sys._dipoles[pos] = new_spin
                    replica.E = E_next

                    if (fac*replica.E) < lim
                        lim = fac * replica.E
                        println(init_output, replica.E)
                        flush(init_output)
                    end

                    if bounds_check(replica.E, bounds)
                        println(init_output, "\n# finish init with E = ",replica.E)
                        close(init_output)
                        return :SUCCESS
                    end
                end
            end
            ln_g_tmp[replica.E] += 1
        end
    end

    println(init_output, "# init failed.")

    close(init_output)
    return :MCSLIMIT
end

"""
Calculate and return {min, max} bounds for energy windows of equal size

# Arguments
-`global_bounds::Vector{Float64}`: {min, max} bouns of total energy range

-`N_wins::Int64`: Number of windows

-`overlap::Float64`: Overlap ratio (recommend > 0.5) between windows
"""
function get_windows1D(
    global_bounds::Vector{Float64}, 
    N_wins::Int64,
    overlap::Float64
)
    Δ = abs(global_bounds[2] - global_bounds[1])

    n = 1 / (1.0 - overlap)

    width = Δ * n / (N_wins + n - 1)

    wins = Vector{Float64}[]
    pos = global_bounds[1]

    for i in 1:N_wins
        w_beg = round(pos, digits=3)
        w_end = round(pos + width, digits=3)
        push!(wins, [w_beg, w_end])
        pos += width / n
    end

    return wins
end

"""
Update min/max energies found and update stats for round-trip times
"""
function update_stats(replica::WLReplica)
    for j in 1:2
        k = 2*j - 3
        if (k*replica.E) > (k*replica.lims_curr[j])
            replica.lims_curr[j] = replica.E
        end 

        if replica.trip_dir == 0
            if abs(replica.E - replica.lims_curr[j]) <= replica.hist.bin_sizes[1]
                replica.trip_dir = -k
                return nothing
            end
        elseif replica.trip_dir == k
            replica.t_dir_curr[j] += 1

            if abs(replica.E - replica.lims_curr[j]) <= replica.hist.bin_sizes[1]
                replica.trip_dir *= -1
                replica.t_dir_avg[j] += replica.t_dir_curr[j]
                replica.t_dir_curr[j] = 0
                replica.N_dir[j] += 1
            end
        end
    end
    return nothing
end

"""
Run a 1D replica-exchange Wang-Landau simulation.

# Arguments
-`sys::SpinSystem`: Sunny system for simulation

-`ranks_per_window::Int64`: Number of WL replicas per window

-`windows_bounds::Vector{Vector{Float64}}`: Vector of {min, max} bounds for each window 

-`bin_size::Float64`: Energy bin size 

-`max_mcs_init::Int64`: Maximum number of MC sweeps to allow during initialization

-`max_mcs::Int64`: Maximum number of MC sweeps to allow during REWL simulation

-`hcheck_interval::Int64`: Number of MC sweeps between histogram checks

-`hcheck_type::Int64`: 1 = min hist is above some factor of avg hist; 2 = hist min is above 1/√lnf

-`hist_flatness::Float64`: if hcheck_type == 1, then min hist entry must be above hist_flatness * (avg hist) 

-`exch_interval::Int64`: Number of MC sweeps between replica exchanges

-`tune_exch_interval::Bool`: Adapt the replica exchange interval to be >= the largest mean round-trip time if true

-`ln_f_final::Float64`: Cutoff threshold for modification factor

-`ln_f_sched::Function`: Function to reduce modification factor at each iteration

-`per_spin::Bool`: Use energy and magnetization per spin if True

-`mc_move_type::String`: "flip" for spin flip (use with IsingSampler); "gaussian" for Gaussian perturbation of spin; 
"spherical_cap" for uniform spin displacemnt within spherical cap

-`mc_step_size::Float64`: Set the update magnitude for gaussian or spherical cap spin updates
"""
function run_REWL!(
    sys::SpinSystem,
    ranks_per_window::Int64,
    windows_bounds::Vector{Vector{Float64}},
    bin_size::Float64;
    max_init_mcs::Int64 = 100_000,
    max_mcs::Int64 = 1_000_000,
    hcheck_interval::Int64 = 10_000,
    hcheck_type::Int64 = 1,
    hist_flatness::Float64 = 0.6,
    exch_interval::Int64 = 1_000,
    tune_exch_interval::Bool = false,
    ln_f_final::Float64 = 1e-6,
    ln_f_sched::F = (ln_f, i)->(0.5*ln_f),
    per_spin::Bool = true,
    mc_move_type::String = "gaussian",
    mc_step_size::Float64 = 0.1
) where {F <: Function}

    # Initialize MPI if not already and get rank
    if !MPI.Initialized()
        rank, total_ranks = init_MPI()
    else
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        total_ranks = MPI.Comm_size(MPI.COMM_WORLD)
    end

    N_wins = length(windows_bounds)
    win_ID = div(rank, ranks_per_window)
    rex_dir = (win_ID % 2 == 0) ? 1 : 2
    bounds = windows_bounds[win_ID+1]

    # Create WL replica and set MPI communicators for intra- and inter-window communication
    replica = WLReplica(sys, bounds, bin_size, rank, N_wins, rex_dir, per_spin)
    init_comms!(replica, ranks_per_window, N_wins)

    # Initialize sys if bounded - must supply [min, max]
    if init_REWL!(
        replica; 
        max_mcs = max_init_mcs, 
        mc_move_type = mc_move_type, 
        mc_step_size = mc_step_size,
        limit_pad = 0.5 * (bounds[2] - bounds[1])
    ) == :MCSLIMIT
        return :INITFAIL
    end

    output = open(@sprintf("R%06d_out.dat", replica.rank), "w")
    println(output, "begin REWL sampling.")

    rex_accepts = zeros(Int64, 2)
    iteration = 1
    N_exch = 0

    trip_output = open(@sprintf("R%06d_mcs_down-up.dat", replica.rank), "w")
    println(trip_output, "# iter \t N_down \t <mcs_down> \t N_up \t <mcs_up>")

    # Initial state
    replica.E = energy(replica.sys) / replica.norm

    # Record initial state
    replica.ln_g[replica.E] += replica.ln_f
    replica.hist[replica.E] += 1

    total_mcs = 0

    for mcs in 1:max_mcs
        for pos in CartesianIndices(replica.sys._dipoles)
            # propose single spin move
            if mc_move_type == "flip"
                new_spin = -replica.sys._dipoles[pos]
            elseif mc_move_type == "gaussian"
                new_spin = gaussian_spin_update(replica.sys._dipoles[pos], mc_step_size)
            elseif mc_move_type == "spherical_cap"
                new_spin = spherical_cap_update(replica.sys._dipoles[pos], mc_step_size)
            end

            # Calculate observables
            E_next = replica.E + local_energy_change(replica.sys, pos, new_spin) / replica.norm

            if bounds_check(E_next, bounds)
                # Add new bin to ln_g, histogram
                if replica.ln_g[E_next] <= eps()
                    add_new!(replica, E_next)
                end

                Δln_g = (replica.ln_g[replica.E] - replica.ln_g[E_next])

                if (Δln_g >= 0) || ( rand() <= exp(Δln_g) )
                    replica.sys._dipoles[pos] = new_spin
                    replica.E = E_next
                end
            end
            # Update ln_g, hist
            replica.ln_g[replica.E] += replica.ln_f 
            replica.hist[replica.E] += 1

            update_stats(replica)
        end
        total_mcs += 1
    
        # Attempt replica exchange
        if mcs % exch_interval == 0
            if replica_exchange!(replica)
                rex_accepts[replica.rex_dir] += 1
            end
            N_exch += 1

            # Alternate exch. direction
            replica.rex_dir = 3 - replica.rex_dir
        end

        # Check histogram criterion - start new iteration if satisfied
        if mcs % hcheck_interval == 0
            # Print current histograms and ln_g to file
            fn = open(@sprintf("R%06d_hist-current.dat", replica.rank), "w")
            print(fn, replica.hist)
            close(fn)

            fn = open(@sprintf("R%06d_ln_g-current.dat", replica.rank), "w")
            print(fn, replica.ln_g)
            close(fn)

            # Print replica exchange rates
            print(output, "replica exchanges: ")
            for i in 1:2
                if replica.inter_comm_ranks[i] >= 0
                    print(output, ((replica.N_rex[i] > 0) ? rex_accepts[i]/replica.N_rex[i] : 0), " ")
                end
            end
            println(output,"")
            flush(output)

            if tune_exch_interval 
                round_trip_time = 0
                for i in 1:2
                    round_trip_time += ((replica.N_dir[i] > 0) ? round(Int64, replica.t_dir_avg[i]/replica.N_dir[i]/length(replica.sys)) : 0)
                end

                all_round_trip_times = MPI.Allgather([round_trip_time], MPI.COMM_WORLD)
                max_round_trip_time = maximum(all_round_trip_times)

                if max_round_trip_time > exch_interval            
                    exch_interval = max_round_trip_time
                    println(output, "change replica exchange interval to ", exch_interval, " mcs")
                end
            end            

            # Histogram is flat
            if check_hist(replica; p=hist_flatness, check_type=hcheck_type)

                @printf(output, "iteration %d complete: mcs = %d, ln_f = %.8f\n", iteration, mcs, replica.ln_f)
                flush(output)
    
                # Average the ln_g and sum up histograms within each window 
                ln_g_avg = MPI.Allreduce(replica.ln_g.vals, MPI.SUM, replica.intra_win_comm) ./ ranks_per_window
                hist_tot = MPI.Allreduce(replica.hist.vals, MPI.SUM, replica.intra_win_comm)

                for i in 1:replica.ln_g.size
                    replica.ln_g.vals[i] = ln_g_avg[i]
                    replica.hist.vals[i] = hist_tot[i]

                    if (!replica.ln_g.visited[i]) && (ln_g_avg[i] > eps())
                        replica.ln_g.visited[i] = true
                        replica.hist.visited[i] = true
                    end
                end 

                # print window combined histogram and average ln_g
                if (replica.intra_comm_rank == 0) && (replica.ln_f > ln_f_final)
                    fn = open(@sprintf("W%05d_hist-iteration%02d.dat", win_ID, iteration), "w")
                    print(fn, replica.hist)
                    close(fn)

                    fn = open(@sprintf("W%05d_ln_g-iteration%02d.dat", win_ID, iteration), "w")
                    print(fn, replica.ln_g)
                    close(fn)
                end

                # print round trip stats and update exchange interval
                if replica.ln_f > ln_f_final
                    print(trip_output, iteration)
                    for i in 1:2
                        t_avg = ((replica.N_dir[i] > 0) ? round(Int64, replica.t_dir_avg[i]/replica.N_dir[i]/length(replica.sys)) : 0)
                        print(trip_output, "\t", replica.N_dir[i], "\t", t_avg)
                    end
                    println(trip_output, "")
                    flush(trip_output)
                end

                # Reset histogram
                reset!(replica.hist)

                rex_accepts .= 0
                replica.N_rex .= 0

                replica.trip_dir = 0
                replica.t_dir_curr .= 0
                replica.t_dir_avg .= 0
                replica.N_dir .= 0

                # Reduce modification factor by some schedule
                replica.ln_f = ln_f_sched(replica.ln_f, iteration)

                iteration += 1

                # Check for termination
                stop = ((replica.ln_f <= ln_f_final) ? 1 : 0)
                if MPI.Allreduce(stop, MPI.PROD, replica.intra_win_comm) > 0
                    break
                end
            end
        end
    end

    println("total mcs = ", total_mcs)

    close(trip_output)
    close(output)
    return :SUCCESS
end
