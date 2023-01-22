""" 
    mutable struct WLReplica2D

Contains data for a Wang-Ladau 'replica' used in 
a 2D replica-exchange simulation
"""
mutable struct WLReplica2D
    # Spin system
    sys::SpinSystem

    # Running energy and magnetization
    E::Float64
    m⃗::Vec3

	# Per-spin normalization
	norm::Int64

    # Binned 2D histogram
    hist::BinnedArrayND{Float64, Int64}

    # Natural log of 2D binned density of states
    ln_g::BinnedArrayND{Float64, Float64}

    win_bounds::Vector{Vector{Float64}}

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

    # Alternating dimension that exchanges are made in
    rex_dim::Int64
end

# Constructor for WL replica that's called internally in 'run_REWL2D!()'
function WLReplica2D(
    sys::SpinSystem, 
    win_bounds::Vector{Vector{Float64}},
    bin_sizes::Vector{Float64},
    rank::Int64,
    N_wins::Int64,
    rex_dir::Int64,
	per_spin::Bool
)
    dim = length(bin_sizes)

    return WLReplica2D(
        sys,
        energy(sys), 
        sum(sys._dipoles),
		(per_spin ? length(sys) : 1),
        BinnedArrayND{Float64,   Int64}(win_bounds; bin_sizes=bin_sizes),        
        BinnedArrayND{Float64, Float64}(win_bounds; bin_sizes=bin_sizes),
        win_bounds,
        1.0,
        rank,
        -1,
        MPI.COMM_NULL,
        -ones(Int64, 2*dim),
        [MPI.COMM_NULL for _ in 1:(2*dim)],
        zeros(Int64, 2*dim),
        rex_dir,
        1
    )
end

# setup MPI comms for inter-window exchanges and intra-window communication
function init_comms!(replica::WLReplica2D, ranks_per_win::Int64, wins_per_dim::Vector{Int64})

    N_wins = prod(wins_per_dim)
    N = ranks_per_win * N_wins

    # Ranks from COMM_WORLD that are used to create groups - not new comm. ranks
    intra_ranks = [ collect(i:i+ranks_per_win-1) for i in 0:ranks_per_win:N-1 ]
	inter_ranks = Vector{Int64}[]
    IDs = [ -ones(Int64, 4) for _ in 1:N_wins ]
	pos = 1

    # Sweep the grid of processes and store neighbors for each dimension
    # -- Replica exchange directions are: (right=1, 2=left, 3=down, 4=up)
    for w in 0:N_wins-1
        # Window Cartesian coordinates
        x = w % wins_per_dim[1]
        y = div(w, wins_per_dim[1])

        # Right neighbor to window 'w'
        if x < wins_per_dim[1]-1
            w_right = w + 1
            push!(inter_ranks, sort(vcat(intra_ranks[w+1], intra_ranks[w_right+1])))
			IDs[w+1][1] = pos
			IDs[w_right+1][2] = pos
			pos += 1
        end
        # Down neighbor to window 'w'
        if y < wins_per_dim[2]-1
            w_down = w + wins_per_dim[1]
            push!(inter_ranks, sort(vcat(intra_ranks[w+1], intra_ranks[w_down+1])))
			IDs[w+1][3] = pos
			IDs[w_down+1][4] = pos
			pos += 1
        end
    end

    w = div(replica.rank, ranks_per_win)
    group_world = MPI.Comm_group(MPI.COMM_WORLD)

    # Create intra-window communicator
    g = MPI.Group_incl(group_world, Int32.(intra_ranks[w+1]))
    replica.intra_win_comm = MPI.Comm_create(MPI.COMM_WORLD, g)
    replica.intra_comm_rank = MPI.Comm_rank(replica.intra_win_comm)

    # Create inter-window communicators
	N_comms = length(inter_ranks)
	inter_comms = MPI.Comm[]

	for i in 1:length(inter_ranks)
		g = MPI.Group_incl(group_world, Int32.(inter_ranks[i]))
		push!(inter_comms, MPI.Comm_create(MPI.COMM_WORLD, g) )
	end

	for j in 1:4
		if IDs[w+1][j] > 0
			replica.inter_win_comms[j] = inter_comms[IDs[w+1][j]]
			replica.inter_comm_ranks[j] = MPI.Comm_rank(replica.inter_win_comms[j])
		end
	end

    return nothing
end

""" 
Check histogram to determine when to advance WL iteration. 
check_type==1 uses the average flatness criterion while 
check_type==2 uses the min histogram > 1/√lnf
"""
function check_hist(replica::WLReplica2D; p::Float64=0.6, check_type::Int64=1)
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
function add_new!(replica::WLReplica2D, keys::Float64...)
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
    replica.ln_g[keys...] = ln_g_min - replica.ln_f
    replica.hist[keys...] = 0

    return nothing
end

function bounds_check(x::Vector{Float64}, bounds::Vector{Vector{Float64}})
    for i in 1:length(x)
        r = round(x[i], digits=10)
        if (r < bounds[i][1]) || (r > bounds[i][2])
            return false
        end
    end

    return true
end

"""
Choose random exchange partner from neighboring window.
Idea for this function from Guangjie Shi.
"""
function get_swap_partner(replica::WLReplica2D)
    ranks_per_win = MPI.Comm_size(replica.intra_win_comm)
    pairs = -ones(Int64, 2*ranks_per_win)

    comm_ID = 2*(replica.rex_dim-1) + replica.rex_dir

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
    end    
	swap_partner = MPI.Scatter(pairs, 1, 0, replica.inter_win_comms[comm_ID])

    return swap_partner[1]
end

"""
Perform a replica exchange between niehgboring windows and accept/reject
according to the joint density of states.
"""
function replica_exchange!(replica::WLReplica2D)
    comm_ID = 2*(replica.rex_dim-1) + replica.rex_dir 

    if replica.inter_comm_ranks[comm_ID] == -1
        return false
    end
    replica.N_rex[comm_ID] += 1
    
    rex_rank = get_swap_partner(replica)

    x_curr = [replica.E, replica.m⃗[3]] ./ replica.norm
    x_exch = zeros(Float64, 2)

    # Exchange energies to see if in bounds
    MPI.Sendrecv!(x_curr, rex_rank, 1, x_exch, rex_rank, 1, replica.inter_win_comms[comm_ID])

    # A ln_g difference of -Inf signals out of bounds
    Δln_g = [-Inf]
    if bounds_check(x_exch, replica.win_bounds)
        Δln_g[1] = replica.ln_g[x_curr...] - replica.ln_g[x_exch...] 
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
    replica.E = x_exch[1] * replica.norm
    replica.m⃗ = sum(replica.sys._dipoles)
    return true
end

""" 
Initialize system to bounded range of states using throw-away WL sampling run.
"""
function init_REWL2D!(
	replica::WLReplica2D;
    max_mcs::Int64 = 100_000,
    mc_move_type::String = "gaussian",
    mc_step_size::Float64 = 0.1
)
    init_output = open(@sprintf("init%06d.dat", replica.rank), "w")
    println(init_output, "# begin init")
	println(init_output, "# bounds = ", replica.hist.bounds, "\n")
	flush(init_output)

    # The state 'x' is a pair of energy and magnetization values
    replica.E = energy(replica.sys) 
    replica.m⃗ = sum(replica.sys._dipoles) 
    dim = 2

    x = [replica.E, replica.m⃗[3]] ./ replica.norm
    x_next = copy(x)

    # Make bounds for temp. ln_g that include current state 
	bounds = replica.hist.bounds
    init_space = deepcopy(bounds)
    x_lims = copy(x)
    fac = zeros(Int64, dim)

    for i in 1:dim
        for j in 1:2
            k = 2*j - 3
            if (k*x_lims[i]) > (k*bounds[i][j])
                init_space[i][j] = x_lims[i]
                fac[i] = k
            end
        end
    end

    # WL is already in bounds
    if norm(fac) <= eps()
        println(init_output, "\n# skip init with E = ", x[1], ", m = ", x[2])
        close(init_output)
        return :SUCCESS
    end

    ln_g_tmp = BinnedArrayND{Float64, Int64}(init_space; bin_sizes=replica.hist.bin_sizes)
    ln_g_tmp[x...] = 1

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
            E_next = replica.E + local_energy_change(replica.sys, pos, new_spin)
            m⃗_next = replica.m⃗ + (new_spin - replica.sys._dipoles[pos])
            x_next .= [E_next, m⃗_next[3]] ./ replica.norm

            if bounds_check(x_next, init_space)
                Δln_g = ln_g_tmp[x...] - ln_g_tmp[x_next...]

                if (Δln_g >= 0) || ( rand() <= exp(Δln_g) )
                    replica.sys._dipoles[pos] = new_spin
                    replica.E = E_next
                    replica.m⃗ = m⃗_next
                    x .= x_next

                    print_flag = false
                    for i in 1:dim
                        if (fac[i] != 0) && ((fac[i]*x[i]) < x_lims[i])
                            x_lims[i] = fac[i] * x[i]
                            print_flag = true
                        end
                    end

                    if print_flag
                        println(init_output, mcs," ",x[1]," ",x[2])
						flush(init_output)
                    end

                    if bounds_check(x, bounds)
                        println(init_output, "\n# finish init with E = ",x[1], ", m = ",x[2])
                        close(init_output)
                        return :SUCCESS
                    end
                end
            end
            ln_g_tmp[x...] += 1
        end
    end
    println(init_output, "init failed.")

    close(init_output)
    return :MCSLIMIT
end

"""
Calculate and return vector of 2D window boundaries for equal-sized windows.

# Arguments
-`global_bounds::Vector{Vector{Float64}}`: {min, max} bounds on observable for each dimension 

-`wins_per_dim::Vector{Int64}`: Number of windows for each dimension

-`overlap_per_dim::Vector{Float64}`: Overlap ratio (recommend > 0.5) for each dimension
"""
function get_windows2D(
    global_bounds::Vector{Vector{Float64}}, 
    wins_per_dim::Vector{Int64},
    overlap_per_dim::Vector{Float64}
)
    Δ = [abs(b[2]-b[1]) for b in global_bounds]

    n = 1 ./ (1.0 .- overlap_per_dim)

    widths = Δ .* n ./ (wins_per_dim .+ n .- 1)

    wins = Vector{Vector{Float64}}[]
    pos = [b[1] for b in global_bounds]

    for i in 1:wins_per_dim[1]
        for j in 1:wins_per_dim[2]
            w_beg = round.(pos, digits=3)
            w_end = round.(pos .+ widths, digits=3)
            push!(wins, [[w_beg[1], w_end[1]], [w_beg[2], w_end[2]]])
            pos[2] += widths[2] / n[2]
        end
        pos[1] += widths[1] / n[1]
        pos[2] = global_bounds[2][1]
    end

    return wins
end

"""
Run a 2D replica-exchange Wang-Landau simulation.

# Arguments
-`sys::SpinSystem`: Sunny system for simulation

-`ranks_per_window::Int64`: Number of WL replicas per window

-`windows_per_dim::Vector{Int64}`: Number of windows in each dim

-`windows_bounds::Vector{Vector{Vector{Float64}}}`: Vector of {min, max} bounds for each dim of each window 

-`bin_sizes::Vector{Float64}`: Bin size for each dim

-`max_mcs::Int64`: Maximum number of MC sweeps to allow

-`hcheck_interval::Int64`: Number of MC sweeps between histogram checks

-`hcheck_type::Int64`: 1 = min hist is above some factor of avg hist; 2 = hist min is above 1/√lnf

-`hist_flatness::Float64`: if hcheck_type == 1, then min hist entry must be above hist_flatness * (avg hist) 

-`exch_interval::Int64`: Number of MC sweeps between replica exchanges

-`ln_f_final::Float64`: Cutoff threshold for modification factor

-`ln_f_sched::Function`: Function to reduce modification factor at each iteration

-`per_spin::Bool`: Use energy and magnetization per spin if True

-`mc_move_type::String`: "flip" for spin flip (use with IsingSampler); "gaussian" for Gaussian perturbation of spin; 
"spherical_cap" for uniform spin displacemnt within spherical cap

-`mc_step_size::Float64`: Set the update magnitude for gaussian or spherical cap spin updates
"""
function run_REWL2D!(
    sys::SpinSystem,
    ranks_per_window::Int64,
    windows_per_dim::Vector{Int64},
    windows_bounds::Vector{Vector{Vector{Float64}}},
    bin_sizes::Vector{Float64};
    max_mcs::Int64 = 1_000_000,
    max_init_mcs::Int64 = 100_000,
    hcheck_interval::Int64 = 10_000,
    hcheck_type::Int64 = 1,
    hist_flatness::Float64 = 0.6,
    exch_interval::Int64 = 1_000,
    ln_f_final::Float64 = 1e-6,
    ln_f_sched::F = (ln_f, i)->(0.5*ln_f),
    per_spin::Bool = true,
    mc_move_type::String = "gaussian",
    mc_step_size::Float64 = 0.1,
) where {F <: Function}

    # Initialize MPI if not already and get rank
    if !MPI.Initialized()
        rank, total_ranks = init_MPI()
    else
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        total_ranks = MPI.Comm_size(MPI.COMM_WORLD)
    end

    # Implemented for dim=2 for now
    dim = length(windows_per_dim)
    N_wins = prod(windows_per_dim)
    # Index and Cartesian coordinates (both 0-based indexing) for window
    win_ID = div(rank, ranks_per_window)
    win_coords = [ win_ID % windows_per_dim[1], div(win_ID, windows_per_dim[1]) ]
    rex_dir = sum(win_coords)%2 + 1
	bounds = windows_bounds[win_ID+1]

    # Create WL replica and set MPI communicators for intra- and inter-window communication
    replica = WLReplica2D(sys, bounds, bin_sizes, rank, N_wins, rex_dir, per_spin)
    init_comms!(replica, ranks_per_window, windows_per_dim)

    # Initialize system if bounded - must supply [min, max]
    if init_REWL2D!(replica; max_mcs=max_init_mcs, mc_move_type=mc_move_type, mc_step_size=mc_step_size) == :MCSLIMIT
        return :INITFAIL
    end
	MPI.Barrier(MPI.COMM_WORLD)

    output = open(@sprintf("R%06d_out.dat", replica.rank), "w")
    println(output, "begin REWL sampling.")

    rex_accepts = zeros(Int64, 2*dim)
    iteration = 1
    N_exch = 0

    # The state 'x' is a pair of energy and magnetization values
    replica.E = energy(replica.sys)
    replica.m⃗ = sum(replica.sys._dipoles)
    x = [replica.E, replica.m⃗[3]] ./ replica.norm
    x_next = copy(x)

    # Record initial state
	println(output, "x = ", x[1],", ",x[2],", bounds = ", replica.win_bounds)
	flush(output)
    replica.ln_g[x...] = replica.ln_f
    replica.hist[x...] = 1

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
            E_next = replica.E + local_energy_change(replica.sys, pos, new_spin)
            m⃗_next = replica.m⃗ + (new_spin - replica.sys._dipoles[pos])
            x_next .= [E_next, m⃗_next[3]] ./ replica.norm

            if bounds_check(x_next, bounds)
				# Add new bin to ln_g, histogram
                if (iteration > 1) && (replica.ln_g[x_next...] <= eps())
                    add_new!(replica, x_next...)
                end

                Δln_g = replica.ln_g[x...] - replica.ln_g[x_next...]

                if (Δln_g >= 0) || ( rand() <= exp(Δln_g) )
                    replica.sys._dipoles[pos] = new_spin
                    replica.E = E_next
                    replica.m⃗ = m⃗_next
                    x .= x_next
                end
            end
            # Update ln_g, hist
            replica.ln_g[x...] += replica.ln_f
            replica.hist[x...] += 1
        end
    
        # Attempt replica exchange
        if mcs % exch_interval == 0
            if replica_exchange!(replica)
                rex_accepts[2*(replica.rex_dim-1) + replica.rex_dir] += 1
            end
            N_exch += 1

            # Alternate exch. direction after cycling through each dimension
            if N_exch % dim == 0
                replica.rex_dir = 3 - replica.rex_dir
            end
            replica.rex_dim = 3 - replica.rex_dim
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
			for i in 1:(2*dim)
				if replica.inter_comm_ranks[i] >= 0
					print(output, rex_accepts[i]/replica.N_rex[i], " ")
				end
			end
			println(output,"")
			flush(output)

            # Histogram is flat
            if check_hist(replica; p=hist_flatness, check_type=hcheck_type)

                @printf(output, "iteration %d complete: mcs = %d, ln_f = %.8f\n", iteration, mcs, replica.ln_f)
				flush(output)
				
                # Average the ln_g within each window 
                ln_g_avg = MPI.Allreduce(replica.ln_g.vals, MPI.SUM, replica.intra_win_comm) ./ ranks_per_window

				for i in 1:replica.ln_g.size
					replica.ln_g.vals[i] = ln_g_avg[i]

					if (!replica.ln_g.visited[i]) && (ln_g_avg[i] > eps())
						replica.ln_g.visited[i] = true
					end
				end	

                # only print hist and ln_g for each window
                if replica.intra_comm_rank == 0
                    fn = open(@sprintf("W%05d_hist-iteration%02d.dat", win_ID, iteration), "w")
                    print(fn, replica.hist)
                    close(fn)

                    fn = open(@sprintf("W%05d_ln_g-iteration%02d.dat", win_ID, iteration), "w")
                    print(fn, replica.ln_g)
                    close(fn)
                end

                # Reset histogram
                reset!(replica.hist)

				rex_accepts .= 0
				replica.N_rex .= 0

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

    close(output)
    return :SUCCESS
end
