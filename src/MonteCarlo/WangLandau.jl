mutable struct WangLandau{N, PR}
    sys::System{N}
    ln_g::BinnedArray{Float64, Float64} # log density of states, ln[g(E)]
    hist::BinnedArray{Float64, Float64} # energy histogram
    bounds::Tuple{Float64, Float64} # energy bounds to sample [E_min, E_max]
    propose::PR     # MC move proposal (function)
    ln_f::Float64   # visiting a site with energy E sets: ln[g(E)] += ln_f
    E::Float64      # current energy state (per spin)

    function WangLandau(; sys::System{N}, bin_size, bounds, propose::PR, ln_f=1.0, max_sweeps_relax=100) where {N, PR}
        drive_system_to_energy_bounds!(sys, bounds, propose, max_sweeps_relax)
        E = energy(sys)/length(all_sites(sys))
        return new{N, PR}(
            sys,
            BinnedArray{Float64, Float64}(; bin_size),
            BinnedArray{Float64, Float64}(; bin_size),
            bounds, propose, ln_f, E
        )
    end
end

# use MCMC sampling to get system into bounded energy range
function drive_system_to_energy_bounds!(sys::System{N}, bounds, propose::PR, max_sweeps) where {N, PR}
    nsites = length(all_sites(sys))

    E0 = energy(sys)

    inbounds(E0/nsites, bounds) && return

    kT = E0/nsites < bounds[1] ? Inf : 0.0

    sampler = LocalSampler(; kT, nsweeps=0.01, propose)
    for _ in 1:max_sweeps
        step!(sys, sampler)
        E = E0+sampler.ΔE
        inbounds(E/nsites, bounds) && return
    end

    if iszero(kT)
        error("Failed to drive system to energy bounds. Try initializing the system to a low-energy configuration.")
    else
        error("Failed to drive system to energy bounds. The energy scale may be unphysically high.")
    end
end

# check average flatness of histogram
function check_flat(H::BinnedArray{Float64, Float64}; p::Float64)
    avg = sum(get_vals(H)) / sum(H.visited)
    return !any(x -> x < p*avg, get_vals(H))
end

# check if energy is in specified bounds
function inbounds(x, bounds)
    return bounds[1] <= x <= bounds[2]
end

# add new energy bin by shifting ln_g to current min. and reset hist
function add_new!(WL::WangLandau, E_new::Float64)
    ln_gs = get_vals(WL.ln_g)
    WL.ln_g[E_new] = isempty(ln_gs) ? 0.0 : minimum(ln_gs)
    reset!(WL.hist)
end

# sample for specified number of sweeps
function sample!(WL::WangLandau, nsweeps::Int64)
    (; sys) = WL
    nsites = length(all_sites(sys))
    WL.E = energy(sys) / nsites
    @assert inbounds(WL.E, WL.bounds)

    # try to fulfill the histogram criterion in set number of MC sweeps
    accepts = 0
    for _ in 1:nsweeps*nsites
        # perform single-spin update and calculate proposal energy
        site = rand(sys.rng, all_sites(sys))
        state = WL.propose(sys, site)
        E_next = WL.E + local_energy_change(sys, site, state) / nsites
        
        if inbounds(E_next, WL.bounds)
            iszero(WL.ln_g[E_next]) && add_new!(WL, E_next)

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
