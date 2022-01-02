"Types for arbitrary samplers, and a simple Metropolis sampling algorithm."

"""
    AbstractSampler

All samplers should subtype this, and implement the following methods
    `set_temp!(sampler::MySampler, kT::Float64)`
    `get_system(sampler::MySampler)
    `sample!(sampler::MySampler)`

    Optionally, can override behavior of these, which by default do full-system
     recalculations on each call.
    `running_energy(sampler::MySampler)`
    `running_mag(sampler::MySampler)`
"""
abstract type AbstractSampler end

running_energy(sampler::S) where {S <: AbstractSampler} = energy(get_system(sampler))
running_mag(sampler::S) where {S <: AbstractSampler} = sum(get_system(sampler))
reset_energy!(sampler::S) where {S <: AbstractSampler} = nothing
reset_mag!(sampler::S) where {S <: AbstractSampler} = nothing

"""
    thermalize!(sampler, num_samples)

`sample!` a sampler a given number of times.
"""
@inline function thermalize!(sampler::S, num_samples::Int) where {S <: AbstractSampler}
    for _ in 1:num_samples
        sample!(sampler)
    end
end


"""
    anneal!(sampler, temp_schedule, step_schedule)

`sample!` a sampler at a series of temperatures, staying at each temperature
  for the number of steps in `step_schedule`.
"""
@inline function anneal!(sampler::S,
                         temp_schedule,
                         step_schedule) where {S <: AbstractSampler}
    for (temp, num_steps) in zip(temp_schedule, step_schedule)
        set_temp!(sampler, temp)
        thermalize!(sampler, num_steps)
    end
end

"""
    anneal!(sampler, temp_function, num_samples)

`sample!` a sampler `num_samples` times, with the sample at timestep `n`
 drawn at a temperature `temp_function(n)`.
"""
@inline function anneal!(sampler::S,
                         temp_function::Function,
                         num_samples::Int) where {S <: AbstractSampler}
    for t in 1:num_steps
        set_temp!(sampler, temp_function(t))
        sample!(sampler) 
    end
end

"""
    MetropolisSampler(sys::SpinSystem, kT::Float64, nsweeps::Int)

A sampler which performs the standard Metropolis Monte Carlo algorithm to
 sample a `SpinSystem` at the requested temperature.

Each single-spin update attempts to move the spin to a random position on
 the unit sphere. One call to `sample!` will attempt to flip each spin
 `nsweeps` times.
"""
mutable struct MetropolisSampler{D, L, Db} <: AbstractSampler
    system     :: SpinSystem{D, L, Db}
    β          :: Float64
    nsweeps    :: Int
    E          :: Float64
    M          :: Vec3
    function MetropolisSampler(sys::SpinSystem{D,L,Db}, kT::Float64, nsweeps::Int) where {D,L,Db}
        @assert kT != 0. "Temperature must be nonzero!"
        new{D, L, Db}(sys, 1.0 / kT, nsweeps, energy(sys), sum(sys))
    end
end

"""
    IsingSampler(sys::SpinSystem, kT::Float64, nsweeps::Int)

A sampler which performs the standard Metropolis Monte Carlo algorithm to
 sample a `SpinSystem` at the requested temperature.

This version differs from `MetropolisSampler` in that each single-spin update
only attempts to completely flip the spin. One call to `sample!` will attempt
to flip each spin `nsweeps` times.

Before construting, be sure that your `SpinSystem` is initialized so that each
spin points along its "Ising-like" axis.
"""
mutable struct IsingSampler{D, L, Db} <: AbstractSampler
    system     :: SpinSystem{D, L, Db}
    β          :: Float64
    nsweeps    :: Int
    E          :: Float64
    M          :: Vec3
    function IsingSampler(sys::SpinSystem{D,L,Db}, kT::Float64, nsweeps::Int) where {D,L,Db}
        @assert kT != 0. "Temperature must be nonzero!"
        new{D, L, Db}(sys, 1.0 / kT, nsweeps, energy(sys), sum(sys))
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
    get_system(sampler)

Returns the `SpinSystem` being updated by the `sampler`.
"""
get_system(sampler::MetropolisSampler) = sampler.system
get_system(sampler::IsingSampler) = sampler.system


@inline function _random_spin() :: Vec3
    n = randn(Vec3)
    return n / norm(n)
end

"""
    sample!(sampler)

Samples `sampler.system` to a new state, under the Boltzmann distribution
 as defined by `sampler.system.hamiltonian`.
"""
function sample!(sampler::MetropolisSampler)
    for _ in 1:sampler.nsweeps
        for idx in CartesianIndices(sampler.system)
            # Try to rotate this spin to somewhere randomly on the unit sphere
            newspin = _random_spin()
            ΔE = local_energy_change(sampler.system, idx, newspin)

            if rand() < exp(-sampler.β * ΔE)
                sampler.E += ΔE
                sampler.M += (newspin - sampler.system[idx])
                sampler.system[idx] = newspin
            end
        end
    end
end

function sample!(sampler::IsingSampler)
    for _ in 1:sampler.nsweeps
        for idx in CartesianIndices(sampler.system)
            # Try to completely flip this spin
            newspin = -sampler.system[idx]
            ΔE = local_energy_change(sampler.system, idx, newspin)

            if rand() < exp(-sampler.β * ΔE)
                sampler.E += ΔE
                sampler.M += 2 * newspin
                sampler.system[idx] = newspin
            end
        end
    end
end

@inline running_energy(sampler::MetropolisSampler) = sampler.E
@inline running_mag(sampler::MetropolisSampler) = sampler.M
@inline reset_running_energy!(sampler::MetropolisSampler) = (sampler.E = energy(sampler.system); nothing)
@inline reset_running_mag!(sampler::MetropolisSampler) = (sampler.M = sum(sampler.system); nothing)
@inline running_energy(sampler::IsingSampler) = sampler.E
@inline running_mag(sampler::IsingSampler) = sampler.M
@inline reset_running_energy!(sampler::IsingSampler) = (sampler.E = energy(sampler.system); nothing)
@inline reset_running_mag!(sampler::IsingSampler) = (sampler.M = sum(sampler.system); nothing)


"""
    local_energy_change(sys, idx, newspin)

Computes the change in energy if we replace the spin at `sys[idx]`
  with `newspin`.
"""
function local_energy_change(sys::SpinSystem{D}, idx, newspin::Vec3) where {D}
    ℋ = sys.hamiltonian
    ΔE = 0.0
    oldspin = sys[idx]
    spindiff = newspin - oldspin
    (i, cell) = splitidx(idx)

    if !isnothing(ℋ.ext_field)
        ΔE -= ℋ.ext_field.effBs[i] ⋅ spindiff
    end
    for heisen in ℋ.heisenbergs
        J = heisen.effJ
        for (bond, _) in sublat_bonds(heisen.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += J * (newspin⋅newspin - oldspin⋅oldspin)
            else
                Sⱼ = sys[bond.j, offset(cell, bond.n, sys.lattice.size)]
                ΔE += J * (spindiff ⋅ Sⱼ)
            end
        end
    end
    for diag_coup in ℋ.diag_coups
        for (bond, J) in sublat_bonds(diag_coup.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += newspin⋅(J.*newspin) - oldspin⋅(J.*oldspin)
            else
                Sⱼ = sys[bond.j, offset(cell, bond.n, sys.lattice.size)]
                ΔE += (J .* spindiff) ⋅ Sⱼ
            end
        end
    end
    for gen_coup in ℋ.gen_coups
        for (bond, J) in sublat_bonds(gen_coup.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += dot(newspin, J, newspin) - dot(oldspin, J, oldspin)
            else
                Sⱼ = sys[bond.j, offset(cell, bond.n, sys.lattice.size)]
                ΔE += dot(spindiff, J, Sⱼ)
            end
        end
    end
    if !isnothing(ℋ.dipole_int)
        throw("Local energy changes not implemented yet for dipole interactions")
    end
    return ΔE
end
