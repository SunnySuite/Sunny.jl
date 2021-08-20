struct ParallelTempSampler <: AbstractSampler
    ...
    ...
    ...
end

function ParallelTempSampler(sys::SpinSystem, temps, nsteps)
    ...
    num_replicas = length(temps)
    replicas = [deepcopy(sys) for _ in num_replicas]

    return ParallelTempSampler(replicas, temps, nsteps, ...)
end

function sample!(sampler::ParallelTempSampler)
    ...
end

function thermalize!(sampler::ParallelTempSampler)
    ...
end

(Any other functions)