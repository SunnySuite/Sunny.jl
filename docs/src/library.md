# Library API

This page describes the public types and functions exported by Sunny. This documentation can be also be accessed using the Julia help system (enter `?` at the Julia command prompt).

Typical Sunny usage will involve the following steps:

1. Create a [`Crystal`](@ref), either by providing explicit geometry information
   or by loading a `.cif` file.
2. Perform symmetry analysis using [`print_symmetry_table`](@ref).
3. Define a list of [`Interactions`](@ref), i.e., terms to be included in the
   Hamiltonian.
4. Specify information for each site through a [`SiteInfo`](@ref) object that
   specifies, e.g., local spin magnitude and ``g``-tensor.
5. Assemble a [`SpinSystem`](@ref) using the crystal, the interactions, the
   dimensions of the simulation box (in unit cells), and the site information.
6. Perform some flavor of Monte Carlo simulation, which is used to sample
   equilibrated spin configurations. The [`LangevinSampler`](@ref) uses a
   continuous dynamics that can very efficiently handle long-range dipole-dipole
   interactions. The [`MetropolisSampler`](@ref) may be more effective in the
   presence of strong anisotropy (e.g. the Ising limit) because it employs local
   moves.
7. Measure the static or dynamical structure factor. For details, see the page
   [Structure Factor Calculations](@ref)


## Function list

```@index
```

```@autodocs
Modules = [Sunny]
```


<!-- 
## Plotting

To reduce package load times, Sunny plotting functions are initially hidden, and only become available when the user explicitly executes "`using GLMakie`". It is a good idea to check that the [GLMakie](https://github.com/JuliaPlots/Makie.jl/tree/master/GLMakie) installation is working correctly (execute "`] test GLMakie`" from the Julia REPL).

```_AT_docs
plot_lattice
plot_spins
plot_bonds
plot_all_bonds
anim_integration
live_integration
live_langevin_integration
```
-->
