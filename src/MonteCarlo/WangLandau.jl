# data for multicanonical ensemble -- either MUCA or WL
mutable struct WangLandau{PR}
    ln_g::BinnedArray{Float64, Float64} # log density of states, ln[g(E)]
    hist::BinnedArray{Float64, Float64} # energy histogram
    bounds::Vector{Float64} # energy bounds to sample [E_min, E_max]
    propose::PR     # MC move proposal (function)
    per_spin::Bool  # choose whether to use energy per spin 
    ln_f::Float64   # factor for insstantaneous modification of ln[g(E)] 
    E::Float64      # current energy state

    function WangLandau(; bin_size=1.0, bounds=Float64[], propose=propose_uniform, per_spin=true, ln_f=0.0)
        new{typeof(propose)}(
            BinnedArray{Float64, Float64}(bin_size=bin_size), 
            BinnedArray{Float64, Float64}(bin_size=bin_size), 
            bounds, propose, per_spin, ln_f, 0.0
        )
    end
end

# copy constructor
function Base.copy(WL::WangLandau{PR}) where{PR}
    return WangLandau(; WL.hist.bin_size, WL.bounds, WL.propose, WL.per_spin, WL.ln_f)
end

# check average flatness of histogram
function check_flat(H::BinnedArray{Float64, Float64}; p::Float64=0.8)
    avg = sum(get_vals(H)) / sum(H.visited)
    return !any(x -> x < p*avg, get_vals(H))
end

# check if energy is in specified bounds
function bounds_check(x::Float64, bounds::Vector{Float64})
    r = round(x, digits=10)
    return (r >= bounds[1]) && (r <= bounds[2])
end

# add new energy bin by shifting ln_g to current min. and reset hist
function add_new!(WL::WangLandau, E_new::Float64)
    !any(WL.ln_g.visited) && (WL.ln_g[E_new] = WL.ln_f)
    
    WL.ln_g[E_new] = minimum(get_vals(WL.ln_g)) - WL.ln_f
    WL.hist[E_new] = 0.0
    reset!(WL.hist)
end

# use WL sampling to get system in bounded energy range -- **will reset MUCA data**
function init_system!(sys::System{N}, WL::WangLandau, nsteps::Int64; limit_pad::Float64=0.0) where N

    fac = WL.per_spin ? 1.0/length(sys.dipoles) : 1.0
    WL.E = fac * energy(sys)

    # adjust current bounds to contain initial state
    init_space = copy(WL.bounds)
    lim = WL.E
    dir = 0
    for (i, j) in enumerate([-1, 1])
        if (j * lim) > (j * WL.bounds[i])
            init_space[i] = lim + (j * limit_pad)
            dir = j
        end
    end
    # already in bounds
    (dir == 0) && (return 0)

    reset!(WL.ln_g)
    WL.ln_g[WL.E] = 1.0

    # try to enter energy bounds with Wang-Landau
    mcs = 0
    for _ in 1:nsteps*length(sys.dipoles)
        # perform single-spin update and calculate proposal energy
        site = rand(sys.rng, all_sites(sys))
        state = WL.propose(sys, site)
        E_next = WL.E + fac * local_energy_change(sys, site, state)
        mcs += 1

        # stay within initial bounds 
        if bounds_check(E_next, init_space)
            # raise/lower boundary edge with some 'pad' to avoid trapping
            if (E_next < WL.bounds[1]) && (E_next-limit_pad > init_space[1])
                init_space[1] = E_next - limit_pad
            end
            if (E_next > WL.bounds[2]) && (E_next+limit_pad < init_space[2])
                init_space[2] = E_next + limit_pad
            end

            Δln_g = WL.ln_g[WL.E] - WL.ln_g[E_next]

            # accept spin update  
            if (Δln_g >= 0) || (rand(sys.rng) <= exp(Δln_g))
                setspin!(sys, state, site)
                WL.E = E_next

                # update limit if new energy found above/below pad
                if (fac * WL.E) < lim
                    lim = fac * WL.E
                    println(mcs, " ", WL.E)
                end
                # initialization complete
                if bounds_check(WL.E, WL.bounds)
                    reset!(WL.ln_g, reset_visited=true)
                    return mcs
                end
            end
        end
        WL.ln_g[WL.E] += 1.0
    end
    # initialization failed
    return mcs
end

# perform  MUCA or WL sampling for specified number of sweeps
function sample!(sys::System{N}, WL::WangLandau, nsteps::Int64) where N
    # choose whether to use energy per spin
    fac = WL.per_spin ? 1.0/length(sys.dipoles) : 1.0
    WL.E = fac * energy(sys)

    # prompt user if system is out of specified bounds
    bounded = (length(WL.bounds) == 2) 
    if bounded && !bounds_check(WL.E, WL.bounds)
        println("initialize system before bounded MC")
        return 0
    end

    # try to fulfill the histogram criterion in set number of MC sweeps
    accepts = 0
    for _ in 1:nsteps*length(sys.dipoles)
        # perform single-spin update and calculate proposal energy
        site = rand(sys.rng, all_sites(sys))
        state = WL.propose(sys, site)
        E_next = WL.E + fac * local_energy_change(sys, site, state)
        
        # stay within bounds if specified
        if !bounded || bounds_check(E_next, WL.bounds)
            # shift ln_g for new state up to current minimum
            if WL.ln_g[E_next] <= eps()
                add_new!(WL, E_next)
            end
            
            Δln_g = WL.ln_g[WL.E] - WL.ln_g[E_next]

            # accept spin update  
            if (Δln_g >= 0) || (rand(sys.rng) <= exp(Δln_g))
                setspin!(sys, state, site)
                WL.E = E_next
                accepts += 1
            end
        end
        WL.ln_g[WL.E] += WL.ln_f
        WL.hist[WL.E] += 1.0
    end
    return accepts
end
