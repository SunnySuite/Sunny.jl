# Sunny interface refactor notes

## Miscellaneous

* Rename `SpinSystem` to `System`

* Replace `sys.lattice` with `positions(sys)`.

## Instantiate `System`

Construct a spin system using

```julia
sys = System(crystal, dims, siteinfos; SUN=true)
```

SU(N) mode is enabled by default. One can optionally restrict to dipole-only model by setting `SUN=false`.

The parameter `siteinfos` is a list of `SiteInfo` objects,

```julia
SiteInfo(site::Int; S=1/2, g=2)
```

one for each symmetry-distinct `site`. The parameter `S` gives the total angular momentum (units of $\hbar$), and `g` is a scalar or tensor that determines conversion from angular momentum to magnetic moment. Information about form factor corrections has been moved to the StructureFactor portion of the code.

If `SUN` is false then `siteinfo.S` determines the norm of dipoles in `sys.dipoles` corresponding to site `siteinfo.i`. If `SUN=true`, then all $S$ must be the same, and the size of each local Hilbert space is $N = 2S + 1$.


## Setting interactions

Terms in the Hamiltonian can be added mutably to an existing `System`. These are applied globally,

```julia
# Zeeman coupling to external field.
set_external_field!(sys::System, h::Vec3)

# Enable or disable long-range dipole-dipole interactions.
enable_dipole_dipole!(sys)
disable_dipole_diople!(sys)
```

Exchange interactions and single-ion anisotropy are propagated to symmetry equivalent bonds and sites,
```julia
# general 3x3 exchange matrix on a Bond(i, j, cell_offset)
set_exchange!(sys::System, b::Bond, J::Mat3)

# biquadratic scalar interaction
set_biquadratic_scalar!(sys::System, b::Bond, J::Float64)

# single-ion anisotropy where op is an abstract operator
set_anisotropy!(sys::System, site::Int, op::SymbolicOperator)
```

## Inhomogeneity and chemical disorder

Inhomogeneity can be introduced using,

```julia
# Set external field at single spin `idx`
set_external_field!(sys, idx::SpinIndex, h::Vec3)

# Introduce a vacancy at single spin `idx`, which effectively rescales the
# spin magnitude to zero.
set_vacancy!(sys, idx::SpinIndex)
```

As a shorthand, we use the type alias
```julia
# The first three indices label the unit cell, and the last index labels
# the site within a unit cell.
const SpinIndex = CartesianIndex{4}
```

It is also possible to set inhomogeneous exchange interactions,

```julia
# First, we must switch to "inhomogeneous" mode. This will slow down
# simulations, and cannot be undone for a given System.
enable_inhomogeneous_exchange!(sys)

# The result is inherited from the original, "homogeneous" system.
get_exchange!(sys, idx1::SpinIndex, idx2::SpinIndex, J::Mat3)

# Can override exchange matrices on individual bonds.
set_exchange!(sys, idx1::SpinIndex, idx2::SpinIndex, J::Mat3)
```


## Structure factor


```julia
add_trajectory!(sf, sys)
```

This will run a dynamical trajectory of a copy of the system `sys`, and accumulate data into `sf`. Allocations are avoided by using buffer space in `sys`.

_**David: Please describe the rest of the new API**_


# Internal details

The following internal changes will be made. These should not be user facing.

## Dynamics

New types:

```julia
struct LangevinHeunP
    T::Float64
    Δt::Float64
    λ::Float64
end

struct ImplicitMidpint
    Δt::Float64
    atol :: Float64
end
```

One integration time-step:

```julia
 step!(sys, dynamics)
```

## Sampling

This interface is still being evolve. There is currently also a notion of `AbstractSampler` which provides the `sample!(...)` method. For example,

```julia
mutable struct LangevinSampler <: AbstractSampler
    integrator :: LangevinHeunP
    nsteps     :: Int
end
```

is distinguished from `LangevinHeunp` in that one call to `sample!(sys, ::LangevinSampler)` corresponds to `nsteps` calls to `step!(sys,::LangevinHeunP)`.

Something similar will be needed for Monte Carlo samplers, e.g.

```julia
struct MonteCarlo
    # ???
end
```

Here `sample!(...)` will track global energy and magnetization quantities in `sys`.
