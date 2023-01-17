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


# These should be deprecated/rewritten. Won't have a System in the sampler going forward.
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


################################################################################
# Langevin Sampler 
################################################################################

"""
    LangevinSampler(integrator::LangevinHeunP, nsteps::Int)

Creates a sampler from a Langevin integrator. `nsteps` determines how many
times `step!` is called using the integrator. `nsteps` should be selected large
enough to ensure that the state of the System after integration
is decorrelated with respect to its initial state.
"""
mutable struct LangevinSampler <: AbstractSampler
    integrator :: LangevinHeunP
    nsteps     :: Int
end

set_temp!(sampler::LangevinSampler, kT) = sampler.integrator.kT = kT
get_temp(sampler::LangevinSampler) = sampler.integrator.kT

function sample!(sys, sampler::LangevinSampler)
    (; integrator, nsteps) = sampler
    for _ in 1:nsteps
        step!(sys, integrator)
    end
    return nothing
end


################################################################################
# Metropolis Sampler 
################################################################################


"""
    MetropolisSampler(sys::System, kT::Float64, nsweeps::Int)

A sampler which performs the standard Metropolis Monte Carlo algorithm to sample
 a `System` at the requested temperature.

Each single-spin update attempts to completely randomize the spin. One call to
 `sample!` will attempt to flip each spin `nsweeps` times.
"""
mutable struct MetropolisSampler{N} <: AbstractSampler
    sys        :: System{N}
    β          :: Float64
    nsweeps    :: Int
    E          :: Float64
    M          :: Vec3
    function MetropolisSampler(sys::System{N}, kT::Float64, nsweeps::Int) where N
        @assert kT != 0. "Temperature must be nonzero!"
        new{N}(sys, 1.0 / kT, nsweeps, energy(sys), sum(sys.dipoles))
    end
end

"""
    IsingSampler(sys::System, kT::Float64, nsweeps::Int)

A sampler which performs the standard Metropolis Monte Carlo algorithm to
sample a `System` at the requested temperature.

This version differs from `MetropolisSampler` in that each single-spin update
only attempts to completely flip the spin. One call to `sample!` will attempt
to flip each spin `nsweeps` times.

Before constructing, be sure that your `System` is initialized so that each
spin points along its "Ising-like" axis.
"""
mutable struct IsingSampler{N} <: AbstractSampler
    sys        :: System{N}
    β          :: Float64
    nsweeps    :: Int
    E          :: Float64
    M          :: Vec3
    function IsingSampler(sys::System{N}, kT::Float64, nsweeps::Int) where N
        @assert kT != 0. "Temperature must be nonzero!"
        new{N}(sys, 1.0 / kT, nsweeps, energy(sys), sum(sys.dipoles))
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


"""
    sample!(sampler)

Samples `sampler.sys` to a new state, under the Boltzmann distribution
 as defined by `sampler.sys.interactions`.
"""
function sample!(sampler::MetropolisSampler{N}) where N
    sys = sampler.sys
    for _ in 1:sampler.nsweeps
        for idx in CartesianIndices(sys.dipoles)
            # Try to rotate this spin to somewhere randomly on the unit sphere
            state = randspin(sys, idx)
            ΔE = local_energy_change(sys, idx, state)

            if rand(sys.rng) < exp(-sampler.β * ΔE)
                sampler.E += ΔE
                sampler.M += state.s - sys.dipoles[idx]

                setspin!(sys, state, idx)
            end
        end
    end
end


function sample!(sampler::IsingSampler{N}) where N
    sys = sampler.sys
    for _ in 1:sampler.nsweeps
        for idx in CartesianIndices(sys.dipoles)
            # Try to completely flip this spin
            state = flip(getspin(sys, idx))
            ΔE = local_energy_change(sys, idx, state)

            if rand(sys.rng) < exp(-sampler.β * ΔE)
                sampler.E += ΔE
                sampler.M += 2 * state.s

                setspin!(sys, state, idx)
            end
        end
    end
end
