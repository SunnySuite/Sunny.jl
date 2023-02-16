"""
    temperature_schedule(kT1, length; curvature=4.0)

Returns an iterator that includes `length` temperatures that exponentially
decays from `kT1` to `kT2`. The optional parameter `curvature` is a multiplier
on the initial slope of the temperature decay.
"""
function temperature_schedule(kT1, kT2, length; curvature=4)
    kT1 >= 0 && kT2 >= 0 || error("Temperature must be positive.")
    error("FIXME")
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
    @propose_mix weight1 propose1 weight2 propose2 ...

Macro to generate a proposal function that randomly selects among the provided
functions according to the provided probability weights. For use with
[`LocalSampler`](@ref).

# Example
```julia
# A proposal function that proposes a spin flip 40% of the time, and a
# Gaussian perturbation 60% of the time.
@propose_mix 0.4 propose_flip 0.6 propose_delta(0.2)
```
"""
macro mix_proposals(terms...)
    isodd(length(terms)) && error("Alternate weights and proposal functions.")
    terms = reshape(collect(terms), 2, :)
    nterms = size(terms, 2)
    return quote
        let
            # Calculative cumulative probabilities for use in branching
            weights = SVector{$nterms, Float64}($(terms[1,:]...))
            probs = weights / sum(weights)
            cumprobs = accumulate(+, probs)
            @assert cumprobs[end] ≈ 1.0

            # Storing the functions in a tuple preserves type information to avoid
            # dynamic dispatch.
            fns = ($(terms[2,:]...),)

            function ret(sys::System{N}, idx) where N
                r = rand(sys.rng)

                $([:(r < cumprobs[$i] && return fns[$i](sys, idx)) for i in 1:nterms-1]...)
                # It is always true that r ≤ cumprobs[end]
                return fns[end](sys, idx)
            end
            ret
        end
    end
end

"""
    LocalSampler(; kT, nsweeps=1.0, propose=propose_uniform)

Monte Carlo simulation involving Metropolis updates to individual spins. Use
with the [`step!`](@ref) function.

 - `kT` is the target temperature, and can be updated mutably.
 - `nsweeps` is the number of full-system MCMC sweeps, and may be fractional.
   The default value of `1.0` means that `step!` performs, on average, one trial
   update for every spin.
 - `propose` is a function to generate new candidate spin states, which may be
   accepted or rejected. Options include [`propose_uniform`](@ref),
   [`propose_flip`](@ref), and [`propose_delta`](@ref). Multiple proposals can
   be mixed with the macro [`@mix_proposals`](@ref).

The returned object stores fields `ΔE` and `Δs`, which represent the cumulative
change to the net energy and dipole, respectively.
"""
mutable struct LocalSampler{Propose}
    kT      :: Float64   # Temperature
    nsweeps :: Float64   # Number of sweeps per step!()
    propose :: Propose   # Function: (System, Site) -> SpinState
    ΔE      :: Float64   # Cumulative energy change
    Δs      :: Vec3      # Cumulative net dipole change

    function LocalSampler(; kT, nsweeps=1.0, propose=propose_uniform)
        new{typeof(propose)}(kT, nsweeps, propose, 0.0, zero(Vec3))
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
