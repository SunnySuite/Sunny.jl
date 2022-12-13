# Sunny interface refactor notes

## SpinSystem

```julia
struct SpinSystem
    crystal::Crystal
    dims::(Int,Int,Int)
    ...
end

sys.lattice -> sys.position
```

Basic idea: adopt a mutable interface with interactions added successively rather than all at once.

### Instantiation
The information necessary to instantiate a SpinSystem will then consist of a `Crystal`, lattice dimensions, and local Hilbert space information. It is also necessary to decalare at this point whether to use SU(_N_) mode or not. Creating a `SpinSystem` would then look like this:

```julia
sys = SpinSystem(crystal, dims, siteinfos; SUN_mode=false)
```

The local Hilbert space information, what we have been calling, `siteinfos` contains: the g-factor, S (the total angular momentum). Note that any information about the ion for form factor corrections will be moved to the StructureFactor portion of the code.

```julia
SiteInfo(site::Int; S=1/2, g=2I)
```

If `SUN_mode=false`, then the `S` values are used as `spin_rescaling` values. If `SUN_mode=true`, then the largest `S` in the list is converted into an `N`, with `N = 2S + 1` (probably converting to an `Int` with a tolerance check). The resulting `N` will be used for all sites, but the appropriate `spin_rescaling` will be calculated for the remaining sites in case not all sites.

For example, suppose there are two symmetry inequivalent sites, one with `S=1` and the other with `S=1/2`. In dipole mode, the `spin_rescaling` value for each site will simply be set to `1` and `1/2` respectively. In SU(_N_) mode, `N` will be set to `3` for all sites, but `spin_rescaling` will be set to `1` for the `S=1` site and to `1/2` for the `S=1/2` site.


## Lattice



### Adding interactions

The Hamiltonian would then be specified by successively adding interactions after creation of the `SpinSystem`:

```julia
add_exchange!(sys::SpinSystem, b::Bond, J)      # where J is a 3x3 matrix
add_anisotropy!(sys::SpinSystem, atom::Int, op) # where op is an abstract operator
```

Global Zeeman coupling
```julia
set_external_field!(sys::SpinSystem, H)         # where H is a 3-vector
```

Enable dipole-dipole interactions globally
```julia
enable_dipole_dipole!(sys::SpinSystem; kmax=???, eta=???)
disable_dipole_diople!(sys)
```

### Preparing for inhomogeneity

Local applied field,
```julia
set_local_field!(sys, i, j, k, atom, H)
```

Local exchange interaction
```julia
set_local_exchange!(sys, i, j, k, atom1, ni, nj, nk, atom2, J)
```


## Integrators/Samplers -> Dynamics


abstract type Dynamics

struct LangevinIntegrator <: Dynamics
    Δt::Float64
    λ::Float64
    T::Float64
end

struct SphericalIntegrator <: Dynamics
    Δt::Float64
end

struct MonteCarlo <: Dynamics
    ...
end


```julia
 step!(sys, dynamics)
```

In the case of MonteCarlo dynamics, the global energy and magnetization will be updated in sys, and a flag will be set to true indicating that that data is reliable. For Langevin samples, the flag will be set to false.

## StructureFactor

```julia
add_trajectory!(sf, sys)
```

This will modify sf but not sys. It will create a "clone" of sys which contains a new configuration array for updating the spins along the trajectory. In the implementation, some buffers in sf should be allocated ahead of time so that add_trajectory!() need not allocate memory.

