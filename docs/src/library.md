# Library

Here, we document all publically exposed types and methods in our Module.
Developers may be interested in further documentation of the [Internals](@ref).

Our package makes extensive usage of
[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl).
In particular, throughout the documentation we make use of aliases
`Vec3 = SVector{3, Float64}`, `Mat3 = SMatrix{3, 3, Float64, 9}`. Additionally, some features
and internals utilize [OffsetArrays.jl](https://github.com/JuliaArrays/OffsetArrays.jl) for
pleasant indexing. In particular, structure factors are often returned as these.

## Geometry definition

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
```

## Interactions

```@docs
easy_axis
easy_plane
single_ion_anisotropy
heisenberg
dm_interaction
exchange
external_field
dipole_dipole
```

## System definition

```@docs
SpinSystem
SpinSystem(::Crystal, ::Vector{<:Sunny.Interaction}, latsize, sites_info)
rand!(::SpinSystem)
energy
field
field!
```

## Sampling

```@docs
LangevinSampler(::SpinSystem, ::Float64, ::Float64, ::Float64, ::Int)
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
should instead perform dynamics either using [`LangevinSampler`](@ref) or implicitly in [Structure factor calculations](@ref). However, advanced users and developers may want direct access to an interface to perform dynamics
integrations.

```@docs
HeunP
LangevinHeunP
evolve!
```