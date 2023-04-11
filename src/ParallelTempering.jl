# vanilla parallel tempering using Julia multithreading
mutable struct ParallelTempering
    n_replicas::Int64
    # temperatures
    kT_sched::Vector{Float64}
    # sampler that evolves/updates systems at temperature
    samplers::Vector{LocalSampler}
    # system (replica) containing state information 
    systems::Vector{System}
    # ids that tells which system is at each sampler
    system_ids::Vector{Int64}
    # acceptance rate for replica exchanges btw (rank, rank+1)
    n_accept::Vector{Int64}
    # number of attempted exchanges
    n_exch::Vector{Int64}
end

function ParallelTempering(system::System, kT_sched::Vector{Float64}, propose)
    n_replicas = length(kT_sched)
    samplers = [LocalSampler(; kT, propose) for kT in kT_sched]
    systems = [clone_system(system) for _ in 1:n_replicas]
    system_ids = collect(1:n_replicas)

    #return ParallelTempering(n_replicas, kT_sched, samplers, systems, system_ids, zeros(n_replicas), 0)
    return ParallelTempering(n_replicas, kT_sched, samplers, systems, system_ids, zeros(n_replicas), zeros(n_replicas))
end

# attempt a replica exchange
function replica_exchange!(PT::ParallelTempering, exch_start::Int64)
    # let one replica in the pair handle the exchange
    @Threads.threads for rank in exch_start : 2 : PT.n_replicas-1
        id₁, id₂ = rank, rank+1 

        # action for current thread and its action under exchange
        S₁  = PT.samplers[id₁].ΔE / PT.samplers[id₁].kT
        S₁′ = PT.samplers[id₂].ΔE / PT.samplers[id₁].kT
        ΔS₁ = S₁ - S₁′

        # action for partner thread and its action under exchange
        S₂  = PT.samplers[id₂].ΔE / PT.samplers[id₂].kT
        S₂′ = PT.samplers[id₁].ΔE / PT.samplers[id₂].kT
        ΔS₂ = S₂ - S₂′

        ln_P = ΔS₁ + ΔS₂

        # acceptance criterion -- exchange labels and energy/magnetization for samplers
        if ln_P >= 0 || rand(PT.systems[PT.system_ids[id₁]].rng) < exp(ln_P)
            PT.system_ids[id₁], PT.system_ids[id₂] = PT.system_ids[id₂], PT.system_ids[id₁]

            PT.samplers[id₁].ΔE, PT.samplers[id₂].ΔE = PT.samplers[id₂].ΔE, PT.samplers[id₁].ΔE
            PT.samplers[id₁].Δs, PT.samplers[id₂].Δs = PT.samplers[id₂].Δs, PT.samplers[id₁].Δs

            PT.n_accept[id₁] += 1
        end
        PT.n_exch[id₁] += 1
    end
    #PT.n_exch += 1
end

# run a parallel tempering simulation for 'nsweeps' MC sweeps
function sample!(PT::ParallelTempering, nsteps::Int64, exch_interval::Int64)
    # set number of sweeps btw replica exchanges
    for rank in 1:PT.n_replicas
        PT.samplers[rank].nsweeps = exch_interval
    end
    n_exch = cld(nsteps, exch_interval)

    # start PT simulation
    for i in 1:n_exch
        # sample the systems at each kT in parallel
        @Threads.threads for rank in 1:PT.n_replicas
            step!(PT.systems[PT.system_ids[rank]], PT.samplers[rank])
        end

        # attempt a replica exchange - alternate exchange direction
        replica_exchange!(PT, (i%2)+1)
    end
end

# start PT simulation to measure average energy
function internal_energy!(PT::ParallelTempering, n_measure, measure_interval, exch_interval)
    U = zeros(PT.n_replicas)
    for i in 1:n_measure
        sample!(PT, measure_interval, exch_interval)
        U  .+= [sampler.ΔE for sampler in PT.samplers]
    end
    return (U ./ n_measure)
end

# using the method described in numpy's gradient documentation
function finite_diff(x, y)
    N = length(x)
    dydx = zeros(N)

    for i in 2:N-1
        dx₋ = x[i] - x[i-1]
        dx₊ = x[i+1] - x[i]
        dydx[i] = (dx₋^2 * y[i+1] + (dx₊^2 - dx₋^2) * y[i] - dx₊^2 * y[i-1]) / (dx₋ * dx₊ * (dx₊ + dx₋))
    end
    dydx[1] = (y[2] - y[1]) / (x[2] - x[1])
    dydx[N] = (y[N] - y[N-1]) / (x[N] - x[N-1])

    return dydx
end

# choose new temperature points so that exchange acceptance rates are equal for all replicas
function update_kT(PT::ParallelTempering, U::Vector{Float64}, dkT′::Float64=1e-3)
    # temperature bounds will be preserved
    kT_min, kT_max = extrema(PT.kT_sched)

    # integration variable
    n_int = round(Int64, (kT_max - kT_min)/dkT′, RoundUp)
    kT′ = collect(range(kT_min, kT_max, length=n_int))

    # get heat capacity w/ finite difference and then fit spline to Cᵥ curve
    spline = interpolate(PT.kT_sched, finite_diff(PT.kT_sched, U), FiniteDifferenceMonotonicInterpolation())
    Cᵥ = abs.(spline.(kT′))

    # density (λ) for reassigning temperatures is related to heat capacity
    λ = .√Cᵥ ./ kT′
    λ ./= sum(λ) * dkT′

    # find new temperature schedule 
    kT_new = [kT_min]
    𝓂 = 0.0
    for i in 2:n_int-1
        𝓂 += λ[i] * dkT′

        if 𝓂 >= (length(kT_new) + 1) / PT.n_replicas
            push!(kT_new, kT′[i])
        end
    end
    push!(kT_new, kT_max)

    return kT_new
end

# iteratively adapt kT_sched so that exchange acceptance rates are constant for all replicas
# result will be stored in PT.kT_sched
function optimize_kT!(PT::ParallelTempering, n_iters, n_measure, measure_interval, exch_interval)
    for i in 1:n_iters
        # estimate average energy with parallel tempering
        U = internal_energy!(PT, n_measure, measure_interval, exch_interval)

        # reset statistics
        PT.n_accept .= PT.n_exch .= 0

        # choose updated temperature schedule
        kT_new = update_kT(PT, U)

        # copy system's configuration from the sampler with nearest temperature for each replica
        closest_sampler_ids = argmin.([abs.(PT.kT_sched .- kT) for kT in kT_new])

        for (i, j) in enumerate(closest_sampler_ids)
            PT.systems[PT.system_ids[i]] = Sunny.clone_system(PT.systems[PT.system_ids[j]])
            PT.samplers[i].ΔE = PT.samplers[j].ΔE
            PT.samplers[i].Δs = PT.samplers[j].Δs
        end

        # reassign new temperatures to replicas
        for i in 1:PT.n_replicas
            PT.samplers[i].kT = kT_new[i]
        end
        PT.kT_sched .= kT_new
    end
    return (PT.n_accept ./ PT.n_exch)
end
