# Library API

This page describes the public types and functions exported by Sunny. This documentation can be also be accessed using the Julia help system (enter `?` at the Julia command prompt).

Typical Sunny usage will involve the following steps:

1. Create a [`Crystal`](@ref), either by providing explicit geometry information or by loading a `.cif` file.
2. Perform using methods such as [`print_bond_table`](@ref) and [`print_allowed_anisotropy`](@ref).
3. Define a list of [Interactions](@ref), i.e., terms to be included in the Hamiltonian.
4. Specify information for each site through a [`SiteInfo`](@ref) object that specifies, e.g., local spin magnitude and ``g``-tensor.
5. Assemble a [`SpinSystem`](@ref) using the crystal, the interactions, the dimensions of the simulation box (in unit cells), and the site information.
6. Perform some flavor of Monte Carlo simulation, which is used to sample equilibrated spin configurations.
7. Measure the static or dynamical structure factor. For this, Sunny
   provides high-level helper functions [`dynamic_structure_factor`](@ref) and
   [`static_structure_factor`](@ref). For more documentation, see [Structure factor
   calculations](@ref).


## Crystal definition

```@docs
Crystal
Crystal(lat_vecs, positions; types, symprec)
Crystal(::AbstractString; symprec)
subcrystal
nbasis
cell_volume
lattice_vectors
lattice_params
```

## Symmetry analysis

```@docs
Bond
displacement
distance
coordination_number
print_bond
print_bond_table
reference_bonds
basis_for_symmetry_allowed_couplings
all_symmetry_related_bonds
all_symmetry_related_bonds_for_atom
all_symmetry_related_couplings
all_symmetry_related_couplings_for_atom
print_allowed_anisotropy
```

## Interactions

```@docs
easy_axis
easy_plane
quadratic_anisotropy
heisenberg
dm_interaction
exchange
external_field
dipole_dipole
```

## System definition

```@docs
SpinSystem
SpinSystem(::Crystal, ::Vector{<:Sunny.AbstractInteraction}, latsize, ::Vector{SiteInfo}; μB, μ0)
SiteInfo
rand!(::SpinSystem{N}) where N
randflips!
energy
field
field!
```

## Sampling

```@docs
LangevinSampler
MetropolisSampler
IsingSampler
set_temp!
sample!
thermalize!
anneal!
```

## Structure factor calculations

For extended details on what these functions compute, and how they do it,
see the page [Structure Factor Calculations](@ref)

```@docs
StructureFactor
Sunny.update!
apply_dipole_factor
dynamic_structure_factor
static_structure_factor
```

## Plotting

To reduce package load times, Sunny plotting functions are initially hidden, and only become available when the user explicitly executes "`using GLMakie`". It is a good idea to check that the [GLMakie](https://github.com/JuliaPlots/Makie.jl/tree/master/GLMakie) installation is working correctly (execute "`] test GLMakie`" from the Julia REPL).


```@docs
plot_lattice
plot_spins
plot_bonds
plot_all_bonds
anim_integration
live_integration
live_langevin_integration
```

## Integrators

These functions are not intended to be used by typical users, who instead
should instead perform dynamics either using [`LangevinSampler`](@ref) or implicitly in [Structure Factor Calculations](@ref). However, advanced users and developers may want direct access to an interface to perform dynamics
integrations.

```@docs
HeunP
LangevinHeunP
SphericalMidpoint
evolve!
```