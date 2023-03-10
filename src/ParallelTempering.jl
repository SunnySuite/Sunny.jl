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
    n_exch::Int64
end

function ParallelTempering(system::System, kT_sched::Vector{Float64}, propose)
    n_replicas = length(kT_sched)
    samplers = [LocalSampler(kT=kT, propose=propose) for kT in kT_sched]
    systems = [clone_system(system) for _ in 1:n_replicas]
    system_ids = collect(1:n_replicas)

    return ParallelTempering(n_replicas, kT_sched, samplers, systems, system_ids, zeros(n_replicas), 0)
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
    end
    PT.n_exch += 1
end

# run a parallel tempering simulation for 'nsweeps' MC sweeps
function sample!(PT::ParallelTempering, n_sweeps::Int64, exch_interval::Int64)
    # set number of sweeps btw replica exchanges
    for rank in 1:PT.n_replicas
        PT.samplers[rank].nsweeps = exch_interval
    end
    n_exch = cld(n_sweeps, exch_interval)

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
