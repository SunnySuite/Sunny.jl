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
lattice_vectors
lattice_params
```

## Symmetry analysis

```@docs
Bond
allowed_J
print_bond_table
all_symmetry_related_bonds
all_symmetry_related_interactions
```

## Interactions

```@docs
ExternalField
onsite_anisotropy
heisenberg
dm_interaction
exchange
QuadraticInteraction
OnSiteQuadratic
DipoleDipole
```

## System definition

```@docs
SpinSystem
SpinSystem(::Crystal, ::Vector{<:Interaction}, latsize, ::Rational{Int})
nbasis
eachcellindex
rand!(::SpinSystem)
energy
field
field!
```

## Sampling

```@docs
LangevinSampler(::SpinSystem, ::Float64, ::Float64, ::Float64, ::Int)
MetropolisSampler
set_temp!
sample!
thermalize!
anneal!
```

## Structure factor calculations

For extended details on what these functions compute, and how they do it,
see the page [Structure Factor Calculations](@ref)

```@docs
DynStructFactor
update!
apply_dipole_factor
dynamic_structure_factor
static_structure_factor
```

## Plotting

All plotting functions rely on a successful user installation of [GLMakie](https://github.com/JuliaPlots/Makie.jl/tree/master/GLMakie), which is dependent on the correct video drivers being installed on your system. To ensure your installation is working correctly, please press `]` in a Julia REPL to access the package manager, then execute `test GLMakie`.

```@docs
plot_lattice
plot_spins
plot_bonds
plot_all_bonds
anim_integration
live_integration
live_langevin_integration
```

<!-- ## Ewald summation

These functions are not intended to be used by typical users, who instead should utilize dipole interactions purely through [`DipoleDipole`](@ref).
However, developers may find the following documentation of the internals useful.

```@docs
ewald_sum_monopole
ewald_sum_dipole
precompute_monopole_ewald_compressed
precompute_dipole_ewald_compressed
contract_monopole_compressed
contract_dipole_compressed
``` -->

## Integrators

These functions are not intended to be used by typical users, who instead
should instead perform dynamics either using [`LangevinSampler`](@ref) or implicitly in [Structure factor calculations](@ref). However, advanced users and developers may want direct access to an interface to perform dynamics
integrations.

```@docs
HeunP
LangevinHeunP
evolve!
```