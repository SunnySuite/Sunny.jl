"""
    propose_uniform

Function to propose a uniformly random spin update in the context of a
[`LocalSampler`](@ref). In `:dipole` mode, the result is a random three-vector
with appropriate normalization. In `:SUN` mode, the result is a random SU(_N_)
coherent state with appropriate normalization.

For low-temperature Monte Carlo simulations, uniform spin proposals can be very
inefficient due to a high probability of rejection in the Metropolis
accept/reject step. Consider also [`Langevin`](@ref) sampling, which is
rejection free.
"""
const propose_uniform = randspin

"""
    propose_flip

Function to propose pure spin flip updates in the context of a
[`LocalSampler`](@ref). Dipoles are flipped as ``𝐒 → -𝐒``. SU(_N_) coherent
states are flipped using the time-reversal operator.
"""
propose_flip(sys::System{N}, site) where N = flip(getspin(sys, site))

"""
    propose_delta(magnitude)

Generate a proposal function that adds a Gaussian perturbation to the existing
spin state. In `:dipole` mode, the procedure is to first introduce a random
three-vector perturbation ``𝐒′ = 𝐒 + |𝐒| ξ`` and then return the properly
normalized spin ``|𝐒| (𝐒′/|𝐒′|)``. Each component of the random vector ``ξ``
is Gaussian distributed with a standard deviation of `magnitude`; the latter is
dimensionless and typically smaller than one. 

In `:SUN` mode, the procedure is analogous, but now involving Gaussian
perturbations to each of the ``N`` complex components of an SU(_N_) coherent
state.

In the limit of very large `magnitude`, this function coincides with
[`propose_uniform`](@ref).

Consider also [`Langevin`](@ref) sampling, which is rejection free.
"""
function propose_delta(magnitude)
    function ret(sys::System{N}, site) where N
        κ = sys.κs[site]
        if N == 0
            S = sys.dipoles[site] + magnitude * κ * randn(sys.rng, Vec3)
            S = normalize_dipole(S, κ)
            return SpinState(S, CVec{0}())        
        else
            Z = sys.coherents[site] + magnitude * sqrt(κ) * randn(sys.rng, CVec{N})
            Z = normalize_ket(Z, κ)
            S = expected_spin(Z)
            return SpinState(S, Z)
        end
    end
    return ret
end

"""
    @mix_proposals weight1 propose1 weight2 propose2 ...

Macro to generate a proposal function that randomly selects among the provided
functions according to the provided probability weights. For use with
[`LocalSampler`](@ref).

# Example
```julia
# A proposal function that proposes a spin flip 40% of the time, and a
# Gaussian perturbation 60% of the time.
@mix_proposals 0.4 propose_flip 0.6 propose_delta(0.2)
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

            function ret(sys::System{N}, site) where N
                r = rand(sys.rng)

                $([:(r < cumprobs[$i] && return fns[$i](sys, site)) for i in 1:nterms-1]...)
                # It is always true that r ≤ cumprobs[end]
                return fns[end](sys, site)
            end
            ret
        end
    end
end

"""
    LocalSampler(; kT, nsweeps=1.0, propose=propose_uniform)

Monte Carlo simulation involving Metropolis updates to individual spins. One
call to the [`step!`](@ref) function will perform `nsweeps` of MCMC sampling for
a provided [`System`](@ref). The default value of `1.0` means that `step!`
performs, on average, one trial update per spin.

Assuming ergodicity, the `LocalSampler` will sample from thermal equilibrium for
the target temperature `kT`. 

The trial spin updates are sampled using the `propose` function. Options include
[`propose_uniform`](@ref), [`propose_flip`](@ref), and [`propose_delta`](@ref).
Multiple proposals can be mixed with the macro [`@mix_proposals`](@ref).

The returned object stores fields `ΔE` and `ΔS`, which represent the cumulative
change to the net energy and dipole, respectively.

!!! warning "Efficiency considerations"  
    Prefer [`Langevin`](@ref) sampling in most cases. Langevin dynamics will
    usually be much more efficient for sampling Heisenberg-like spins that vary
    continuously. `LocalSampler` is most useful for sampling from discrete spin
    states. In particular, [`propose_flip`](@ref) may be required for sampling
    Ising-like spins that arise due to a strong easy-axis anisotropy. For strong
    but finite single-ion anisotropy, consider alternating between `Langevin`
    and `LocalSampler` update steps.
"""
mutable struct LocalSampler{F}
    kT      :: Float64   # Temperature
    nsweeps :: Float64   # Number of MCMC sweeps per `step!`
    propose :: F         # Function: (System, Site) -> SpinState
    ΔE      :: Float64   # Cumulative energy change
    ΔS      :: Vec3      # Cumulative net dipole change

    function LocalSampler(; kT, nsweeps=1.0, propose=propose_uniform)
        new{typeof(propose)}(kT, nsweeps, propose, 0.0, zero(Vec3))
    end
end

function Base.copy(sampler::LocalSampler{F}) where F
    LocalSampler(; sampler.kT, sampler.nsweeps, sampler.propose)
end

function step!(sys::System{N}, sampler::LocalSampler) where N
    niters = round(Int, sampler.nsweeps*length(sys.dipoles), RoundUp)
    for _ in 1:niters
        site = rand(sys.rng, eachsite(sys))
        state = sampler.propose(sys, site)
        ΔE = local_energy_change(sys, site, state)

        # Metropolis acceptance probability
        if iszero(sampler.kT)
            accept = ΔE <= 0
        else
            accept = rand(sys.rng) <= exp(-ΔE/sampler.kT)
        end
        
        if accept
            sampler.ΔE += ΔE
            sampler.ΔS += state.S - sys.dipoles[site]
            setspin!(sys, state, site)
        end
    end
end
