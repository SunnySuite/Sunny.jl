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
    sys     :: SpinSystem{N}
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


""" The following functions have been added or modified to permit code reuse
of the basic MC framework. The newly sampled or flipped spin may be either 
a Vec3 or a CVec{N}, and the rest of the code dispatches accordingly with little
alteration to the prior code.

It may be more sensible to more completely separate the functions for classical
LL-type systems and the new SU(N) stuff. Let me know if you prefer this approach.
"""

## For Metropolis sampler
@inline function _random_spin(sys::SpinSystem{0}, idx::CartesianIndex{4}) :: Vec3
    n = randn(sys.rng, Vec3)
    return sys.κs[idx[4]] * normalize(n)
end
@inline function _random_spin(sys::SpinSystem{N}, idx::CartesianIndex{4}) :: CVec{N} where N
    n = randn(sys.rng, CVec{N})
    return normalize(n)  # Kets are always normalized to 1
end

# For Ising sampler. Non-mutating. Hence "flipped" not "flip"
@inline function _flipped_spin(sys::SpinSystem{0}, idx) :: Vec3
    -sys.dipoles[idx]
end
# Uses time-reversal approach to spin flip -- see Sakurai (3rd ed.), eq. 4.176.
@inline function _flipped_spin(sys::SpinSystem{N}, idx) ::CVec{N} where N
    flip_ket(sys.coherents[idx])
end


##### KBTODO: UNIFY
## Since we don't know what type the spin is in the main functions below,
## always generate both a ket and a dipole from the spin. (The is trivial
## except in the case when getting the dipole from and SU(N) spin, since in
## that case it is necessary to coordinate the dipole values when the coherents
## state is updated.)
@inline function ket(::Vec3) :: CVec{0}
    CVec{0}()
end
@inline function ket(spin::CVec{N}) :: CVec{N} where N
    spin
end
@inline function dipole(spin::Vec3, sys::SpinSystem, idx) :: Vec3
    spin
end
@inline function dipole(spin::CVec{N}, sys::SpinSystem, idx) :: Vec3 where N
    sys.κs[idx[4]] * expected_spin(spin)
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
            new_spin = _random_spin(sys, idx[4])
            ΔE = local_energy_change(sys, idx, new_spin)

            if rand(sys.rng) < exp(-sampler.β * ΔE)
                new_ket = ket(new_spin)
                new_dipole = dipole(new_spin, sys, idx)

                sampler.E += ΔE
                sampler.M += (new_dipole - sys.dipoles[idx])

                sys.dipoles[idx] = new_dipole 
                sys.coherents[idx] = new_ket
            end
        end
    end
end


function sample!(sampler::IsingSampler{N}) where N
    sys = sampler.sys
    for _ in 1:sampler.nsweeps
        for idx in CartesianIndices(sys.dipoles)
            # Try to completely flip this spin
            new_spin = _flipped_spin(sys, idx)
            ΔE = local_energy_change(sys, idx, new_spin)

            if rand(sys.rng) < exp(-sampler.β * ΔE)
                new_ket = ket(new_spin)
                new_dipole = dipole(new_spin, sys, idx)

                sampler.E += ΔE
                sampler.M += 2 * new_dipole   # check this is still sensible...

                sys.dipoles[idx] = new_dipole 
                sys.coherents[idx] = new_ket
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

# KBTODO: move this somewhere else
"""
    local_energy_change(sys, idx, newspin)

Computes the change in energy if we replace the spin at `sys[idx]`
  with `newspin`.
"""
function local_energy_change(sys::SpinSystem{N}, idx, newspin) where N
    ℋ = sys.hamiltonian
    ΔE = 0.0
    new_dipole = dipole(newspin, sys, idx)
    new_ket = ket(newspin)
    old_dipole = sys.dipoles[idx]
    spindiff = new_dipole - old_dipole
    cell, i = splitidx(idx)

    if !isnothing(ℋ.ext_field)
        ΔE -= ℋ.ext_field.effBs[i] ⋅ spindiff
    end
    for heisen in ℋ.heisenbergs
        J = first(heisen.bondtable.data)
        for (bond, _) in sublat_bonds(heisen.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += J * (new_dipole⋅new_dipole - old_dipole⋅old_dipole)
            else
                Sⱼ = sys.dipoles[offsetc(cell, bond.n, sys.size), bond.j]
                ΔE += J * (spindiff ⋅ Sⱼ)
            end
        end
    end
    for diag_coup in ℋ.diag_coups
        for (bond, J) in sublat_bonds(diag_coup.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += new_dipole⋅(J.*new_dipole) - old_dipole⋅(J.*old_dipole)
            else
                Sⱼ = sys.dipoles[offsetc(cell, bond.n, sys.size), bond.j]
                ΔE += (J .* spindiff) ⋅ Sⱼ
            end
        end
    end
    for gen_coup in ℋ.gen_coups
        for (bond, J) in sublat_bonds(gen_coup.bondtable, i)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += dot(new_dipole, J, new_dipole) - dot(old_dipole, J, old_dipole)
            else
                Sⱼ = sys.dipoles[offsetc(cell, bond.n, sys.size), bond.j]
                ΔE += dot(spindiff, J, Sⱼ)
            end
        end
    end
    if !isnothing(ℋ.dipole_aniso)
        # TODO: Add tests!
        aniso = ℋ.dipole_aniso
        for (site, J) in zip(aniso.sites, aniso.Js)
            if site == i
                c2, c4, c6 = aniso.coeff_2[i], aniso.coeff_4[i], aniso.coeff_6[i]
                E_new, _ = energy_and_gradient_for_classical_anisotropy(new_dipole, c2, c4, c6)
                E_old, _ = energy_and_gradient_for_classical_anisotropy(old_dipole, c2, c4, c6)
                ΔE += E_new - E_old
            end
        end
    end
    if N > 0
        aniso = ℋ.sun_aniso
        old_ket = sys.coherents[idx]
        Λ = @view(aniso[:,:,i])
        ΔE += real(new_ket' * Λ * new_ket) - real(old_ket' * Λ * old_ket)
    end
    if !isnothing(ℋ.ewald)
        error("Local energy changes not implemented yet for dipole interactions")
    end
    return ΔE
end
