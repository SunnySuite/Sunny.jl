# vanilla parallel tempering using Julia multithreading
mutable struct ParallelTempering{N, SMP}
    n_replicas::Int64
    # temperatures
    kT_sched::Vector{Float64}
    #sampler that evolves/updates systems at temperature
    samplers::Vector{SMP}
    # system (replica) containing state information 
    systems::Vector{System{N}}
    # ids that tells which system is at each sampler
    system_ids::Vector{Int64}
    # acceptance rate for replica exchanges btw (rank, rank+1)
    n_accept::Vector{Int64}
    # number of attempted exchanges
    n_exch::Vector{Int64}
end

function ParallelTempering(system::System{N}, sampler::SMP, kT_sched::Vector{Float64}) where {N, SMP}
    n_replicas = length(kT_sched)

    samplers = SMP[copy(sampler) for _ in 1:n_replicas]
    setproperty!.(samplers, :kT, kT_sched)
 
    systems = [clone_system(system) for _ in 1:n_replicas]
    system_ids = collect(1:n_replicas)

    # set initial energies/magnetizations and set number of sweeps to 1
    if samplers[1] isa LocalSampler
        for rank in 1:n_replicas
            samplers[rank].ΔE = energy(systems[rank])
            samplers[rank].ΔS = sum(systems[rank].dipoles)
            samplers[rank].nsweeps = 1.0
        end
    end

    return ParallelTempering(n_replicas, kT_sched, samplers, systems, system_ids, zeros(Int64, n_replicas), zeros(Int64, n_replicas))
end

# attempt a replica exchange
function replica_exchange!(PT::ParallelTempering, exch_start::Int64)
    # let one replica in the pair handle the exchange
    @Threads.threads for rank in exch_start : 2 : PT.n_replicas-1
        id₁, id₂ = rank, rank+1 

        E₁ = (PT.samplers[1] isa LocalSampler) ? PT.samplers[id₁].ΔE : energy(PT.systems[PT.system_ids[id₁]])
        E₂ = (PT.samplers[1] isa LocalSampler) ? PT.samplers[id₂].ΔE : energy(PT.systems[PT.system_ids[id₂]])

        # action for current thread and its action under exchange
        S₁  = E₁ / PT.samplers[id₁].kT
        S₁′ = E₂ / PT.samplers[id₁].kT
        ΔS₁ = S₁ - S₁′

        # action for partner thread and its action under exchange
        S₂  = E₂ / PT.samplers[id₂].kT
        S₂′ = E₁ / PT.samplers[id₂].kT
        ΔS₂ = S₂ - S₂′

        ln_P = ΔS₁ + ΔS₂

        # acceptance criterion -- exchange labels and energy/magnetization for samplers
        if ln_P >= 0 || rand(PT.systems[PT.system_ids[id₁]].rng) < exp(ln_P)
            PT.system_ids[id₁], PT.system_ids[id₂] = PT.system_ids[id₂], PT.system_ids[id₁]

            if PT.samplers[1] isa LocalSampler
                PT.samplers[id₁].ΔE, PT.samplers[id₂].ΔE = PT.samplers[id₂].ΔE, PT.samplers[id₁].ΔE
                PT.samplers[id₁].ΔS, PT.samplers[id₂].ΔS = PT.samplers[id₂].ΔS, PT.samplers[id₁].ΔS
            end

            PT.n_accept[id₁] += 1
        end
        PT.n_exch[id₁] += 1
    end
end

"""
"""
function step_ensemble! end

# run a parallel tempering simulation for 'nsweeps' MC sweeps
function step_ensemble!(PT::ParallelTempering, nsteps::Int64, exch_interval::Int64)
    n_exch = cld(nsteps, exch_interval)

    # start PT simulation
    for i in 1:n_exch
        # sample the systems at each kT in parallel
        @Threads.threads for rank in 1:PT.n_replicas
            for _ in 1:exch_interval
                step!(PT.systems[PT.system_ids[rank]], PT.samplers[rank])
            end
        end

        # attempt a replica exchange - alternate exchange direction
        replica_exchange!(PT, (i%2)+1)
    end
end
