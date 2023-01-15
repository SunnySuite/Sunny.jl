"Types for arbitrary samplers, and a simple Metropolis sampling algorithm."

"""
    AbstractSampler

Samplers should provide a field `sys::System` and implement the following methods
    `set_temp!(sampler::MySampler, kT::Float64)`
    `get_temp(sampler::MySampler) :: Float64`
    `sample!(sampler::MySampler)`

    Optionally, can override behavior of these, which by default do full-system
     recalculations on each call.
    `running_energy(sampler::MySampler)`
    `running_mag(sampler::MySampler)`
"""
abstract type AbstractSampler end

# These should be deprecated/rewritten. Won't have a SpinSystem in the sampler going forward.
running_energy(sampler::S) where {S <: AbstractSampler} = energy(sampler.sys)
running_mag(sampler::S) where {S <: AbstractSampler} = sum(sampler.sys) # Does this apply g-factor?
reset_running_energy!(sampler::S) where {S <: AbstractSampler} = nothing
reset_running_mag!(sampler::S) where {S <: AbstractSampler} = nothing

"""
    thermalize!(sampler, num_samples)

`sample!` a sampler a given number of times.
"""
function thermalize!(sampler::S, num_samples::Int) where {S <: AbstractSampler}
    for _ in 1:num_samples
        sample!(sampler)
    end
end


"""
    anneal!(sampler, temp_schedule, step_schedule)

`sample!` a sampler at a series of temperatures, staying at each temperature
  for the number of steps in `step_schedule`.

    anneal!(sampler, temp_function::Function, num_samples)

`sample!` a sampler `num_samples` times, with the sample at timestep `n`
 drawn at a temperature `temp_function(n)`.
"""
function anneal!(sampler::S,
                         temp_schedule,
                         step_schedule) where {S <: AbstractSampler}
    for (temp, num_steps) in zip(temp_schedule, step_schedule)
        set_temp!(sampler, temp)
        thermalize!(sampler, num_steps)
    end
end

function anneal!(sampler::S,
                         temp_function::Function,
                         num_samples::Int) where {S <: AbstractSampler}
    for t in 1:num_samples
        set_temp!(sampler, temp_function(t))
        sample!(sampler) 
    end
end

"""
    MetropolisSampler(sys::SpinSystem, kT::Float64, nsweeps::Int)

A sampler which performs the standard Metropolis Monte Carlo algorithm to sample
 a `SpinSystem` at the requested temperature.

Each single-spin update attempts to completely randomize the spin. One call to
 `sample!` will attempt to flip each spin `nsweeps` times.
"""
mutable struct MetropolisSampler{N} <: AbstractSampler
    sys        :: SpinSystem{N}
    β          :: Float64
    nsweeps    :: Int
    E          :: Float64
    M          :: Vec3
    function MetropolisSampler(sys::SpinSystem{N}, kT::Float64, nsweeps::Int) where N
        @assert kT != 0. "Temperature must be nonzero!"
        new{N}(sys, 1.0 / kT, nsweeps, energy(sys), sum(sys.dipoles))
    end
end

"""
    IsingSampler(sys::SpinSystem, kT::Float64, nsweeps::Int)

A sampler which performs the standard Metropolis Monte Carlo algorithm to
sample a `SpinSystem` at the requested temperature.

This version differs from `MetropolisSampler` in that each single-spin update
only attempts to completely flip the spin. One call to `sample!` will attempt
to flip each spin `nsweeps` times.

Before constructing, be sure that your `SpinSystem` is initialized so that each
spin points along its "Ising-like" axis.
"""
mutable struct IsingSampler{N} <: AbstractSampler
    sys     :: SpinSystem{N}
    β          :: Float64
    nsweeps    :: Int
    E          :: Float64
    M          :: Vec3
    function IsingSampler(sys::SpinSystem{N}, kT::Float64, nsweeps::Int) where N
        @assert kT != 0. "Temperature must be nonzero!"
        new{N}(sys, 1.0 / kT, nsweeps, energy(sys), sum(sys.dipoles))
    end
end


"""
    set_temp!(sampler, kT)

Changes the temperature of the sampler to `kT`.
"""
function set_temp!(sampler::MetropolisSampler, kT)
    sampler.β = 1 / kT
end
function set_temp!(sampler::IsingSampler, kT)
    sampler.β = 1 / kT
end

"""
    get_temp(sampler) :: Float64

Returns the temperature of the sampler, as `kT`.
"""
get_temp(sampler::MetropolisSampler) = 1 / sampler.β
get_temp(sampler::IsingSampler) = 1 / sampler.β


# KBTODO: Unify with functions in System...

struct SpinState{N}
    s::Vec3
    Z::CVec{N}
end

@inline function random_state(sys::SpinSystem{0}, idx::CartesianIndex{4})
    s = sys.κs[idx[4]] * normalize(randn(sys.rng, Vec3))
    return SpinState(s, CVec{0}())
end

@inline function random_state(sys::SpinSystem{N}, idx::CartesianIndex{4}) where N
    Z = normalize(randn(sys.rng, CVec{N}))
    s = sys.κs[idx[4]] * expected_spin(Z)
    return SpinState(s, Z)
end

@inline function flipped_state(sys::SpinSystem{N}, idx) where N
    SpinState(-sys.dipoles[idx], flip_ket(sys.coherents[idx]))
end

function local_energy_change(sys::SpinSystem{N}, idx, state::SpinState) where N
    energy_local_delta(sys.dipoles, sys.coherents, sys.hamiltonian, sys.κs, idx, state.s, state.Z)
end


"""
    sample!(sampler)

Samples `sampler.sys` to a new state, under the Boltzmann distribution
 as defined by `sampler.sys.hamiltonian`.
"""
function sample!(sampler::MetropolisSampler{N}) where N
    sys = sampler.sys
    for _ in 1:sampler.nsweeps
        for idx in CartesianIndices(sys.dipoles)
            # Try to rotate this spin to somewhere randomly on the unit sphere
            state = random_state(sys, idx)
            ΔE = local_energy_change(sys, idx, state)

            if rand(sys.rng) < exp(-sampler.β * ΔE)
                sampler.E += ΔE
                sampler.M += state.s - sys.dipoles[idx]

                sys.dipoles[idx] = state.s
                sys.coherents[idx] = state.Z
            end
        end
    end
end


function sample!(sampler::IsingSampler{N}) where N
    sys = sampler.sys
    for _ in 1:sampler.nsweeps
        for idx in CartesianIndices(sys.dipoles)
            # Try to completely flip this spin
            state = flipped_state(sys, idx)
            ΔE = local_energy_change(sys, idx, state)

            if rand(sys.rng) < exp(-sampler.β * ΔE)
                sampler.E += ΔE
                sampler.M += 2 * state.s

                sys.dipoles[idx] = state.s 
                sys.coherents[idx] = state.Z
            end
        end
    end
end


@inline running_energy(sampler::MetropolisSampler) = sampler.E
@inline running_mag(sampler::MetropolisSampler) = sampler.M
@inline reset_running_energy!(sampler::MetropolisSampler) = (sampler.E = energy(sampler.sys); nothing)
@inline reset_running_mag!(sampler::MetropolisSampler) = (sampler.M = sum(sampler.sys.dipoles); nothing)
@inline running_energy(sampler::IsingSampler) = sampler.E
@inline running_mag(sampler::IsingSampler) = sampler.M
@inline reset_running_energy!(sampler::IsingSampler) = (sampler.E = energy(sampler.sys); nothing)
@inline reset_running_mag!(sampler::IsingSampler) = (sampler.M = sum(sampler.sys.dipoles); nothing)
