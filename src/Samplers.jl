"""
    temperature_schedule(kT1, length; curvature=4.0)

Returns an iterator that includes `length` temperatures that decay from `kT1` to
`kT2`. The optional parameter `curvature` controls the initial slope of the
temperature decay.

See also [`anneal!`](@ref).
"""
function temperature_schedule(kT1, kT2, length; curvature=4)
    kT1 >= 0 && kT2 >= 0 || error("Temperature must be positive.")
    α = exp(-curvature)
    c = log((kT2+α)/(kT1+α))
    return (max((kT1+α)*exp(c*x)-α, 0.0) for x in range(0, 1, length))
end

"""
    anneal!(sys::System, sampler, kTs)

Anneals the system by repeated calling `step!(sys, sampler)` using a provided
temperature schedule. Frequently the `kTs` will be constructed using
[`temperature_schedule`](@ref).
"""
function anneal!(sys::System{N}, sampler, kTs) where N
    Es = Float64[]
    for kT in kTs
        sampler.kT = kT
        step!(sys, sampler)
        push!(Es, energy(sys))
    end
    return Es
end


"""
    propose_uniform

Function to propose a uniformly random spin update. For use with
[`LocalSampler`](@ref).
"""
const propose_uniform = randspin

"""
    propose_flip

Function to propose Ising spin flip updates. For use with
[`LocalSampler`](@ref).
"""
propose_flip(sys::System{N}, idx) where N = flip(getspin(sys, idx))

"""
    propose_delta(magnitude)

Generate a proposal function that adds a Gaussian perturbation to the existing
spin. The `magnitude` is typically order one or smaller. For use with
[`LocalSampler`](@ref).
"""
function propose_delta(magnitude)
    function ret(sys::System{N}, idx) where N
        κ = sys.κs[idx]
        if N == 0
            s = sys.dipoles[idx] + magnitude * κ * randn(sys.rng, Vec3)
            s = normalize_dipole(s, κ)
            return SpinState(s, CVec{0}())        
        else
            Z = sys.coherents[idx] + magnitude * sqrt(κ) * randn(sys.rng, CVec{N})
            Z = normalize_ket(Z, κ)
            s = expected_spin(Z)
            return SpinState(s, Z)
        end
    end
    return ret
end

"""
    propose_mix(options, weights)

Generate a proposal function that randomly selects among all `options`,
according to the probability `weights`. For use with [`LocalSampler`](@ref).

# Example
```julia
propose_mix([propose_flip, propose_delta(0.2)], [0.5, 0.5])
```
"""
function propose_mix(options, weights)
    if any(<(0), weights)
        error("Weights must be positive.")
    end
    probabilities = collect(weights) / sum(weights)
    accum = accumulate(+, probabilities)
    @assert accum[end] ≈ 1.0
    accum[end] = 1.0
    function ret(sys::System{N}, idx) where N
        r = rand(sys.rng)
        i = findfirst(>(r), accum)::Int
        return options[i](sys, idx)
    end
    return ret
end

"""
    LocalSampler(; kT, nsweeps=1.0, propose=propose_uniform)

Monte Carlo simulation involving Metropolis updates to individual spins. Use
with the [`step!`](@ref) function.

 - `kT` is the target temperature, and can be updated dynamically.
 - `nsweeps` is the number of full-system MCMC sweeps, and may be fractional.
   The default value of `1.0` means that `step!` performs, on average, one trial
   update for every spin.
 - `propose` is a function to generate new candidate spin states, which may be
   accepted or rejected. Options include [`propose_uniform`](@ref),
   [`propose_flip`](@ref), and [`propose_delta`](@ref). Multiple proposals can
   be mixed with [`propose_mix`](@ref).
"""
mutable struct LocalSampler{Propose}
    kT      :: Float64   # Temperature
    nsweeps :: Float64   # Number of sweeps per step!()
    propose :: Propose   # Function: (System, Site) -> SpinState
    ΔE      :: Float64   # Cumulative energy change
    Δs      :: Vec3      # Cumulative net dipole change

    function LocalSampler(; kT, nsweeps=1.0, propose=propose_uniform)
        new{typeof(propose)}(kT, nsweeps, propose, 0.0, 0.0)
    end
end


function step!(sys::System{N}, sampler::LocalSampler) where N
    niters = round(Int, sampler.nsweeps*length(sys.dipoles), RoundUp)
    for _ in 1:niters
        idx = rand(sys.rng, all_sites(sys))
        state = sampler.propose(sys, idx)
        ΔE = local_energy_change(sys, idx, state)

        # Metropolis acceptance probability
        if rand(sys.rng) < exp(-ΔE/sampler.kT)
            sampler.ΔE += ΔE
            sampler.Δs += state.s - sys.dipoles[idx]
            setspin!(sys, state, idx)
        end
    end
end
