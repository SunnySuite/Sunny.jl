# Library Reference

This page describes the public types and functions exported by Sunny. This documentation can be also be accessed using the Julia help system (enter `?` at the Julia command prompt).

## Index

```@index
Modules = [Sunny]
```

## Core Sunny functions

```@docs
BinningParameters
Bond
Crystal
FormFactor
ImplicitMidpoint
Langevin
LocalSampler
Moment
SampledCorrelations
SampledCorrelationsStatic
SCGA
Site
SpinWaveTheory
SpinWaveTheoryKPM
SpinWaveTheorySpiral
System
Units
add_sample!
clone_correlations
clone_system
copy_spins!
dispersion
dmvec
domain_average
eachsite
enable_dipole_dipole!
energy
energy_per_site
excitations
excitations!
gaussian
global_position
intensities
intensities_bands
intensities_static
lattice_params
lattice_vectors
load_nxs
lorentzian
magnetic_moment
merge_correlations
minimize_energy!
minimize_spiral_energy!
modify_exchange_with_truncated_dipole_dipole!
polarize_spins!
position_to_site
powder_average
primitive_cell
print_bond
print_irreducible_bz_paths
print_site
print_stevens_expansion
print_suggested_frame
print_symmetry_table
print_wrapped_intensities
propose_delta
propose_flip
propose_uniform
q_space_grid
q_space_path
randomize_spins!
reference_bonds
remove_periodicity!
repeat_periodically
repeat_periodically_as_spiral
reshape_supercell
resize_supercell
rotate_operator
set_coherent!
set_dipole!
set_dipoles_from_mcif!
set_exchange!
set_exchange_at!
set_field!
set_field_at!
set_onsite_coupling!
set_onsite_coupling_at!
set_pair_coupling!
set_pair_coupling_at!
set_spin_rescaling!
set_spin_rescaling_for_static_sum_rule!
set_spin_s_at!
set_vacancy_at!
spin_label
spin_matrices
spiral_energy
spiral_energy_per_site
ssf_custom
ssf_custom_bm
ssf_perp
ssf_trace
standardize
step!
stevens_matrices
subcrystal
suggest_magnetic_supercell
suggest_timestep
symmetry_equivalent_bonds
to_inhomogeneous
to_product_space
@mix_proposals
```

## Optional Makie extensions

Load a Makie graphics package (`GLMakie`, `WGLMakie`, or `CairoMakie`) to enable
the following extensions:

```@docs
view_bz
view_crystal
plot_intensities
plot_intensities!
plot_spins
plot_spins!
```

## Optional WriteVTK extensions

Load the `WriteVTK` package to enable the following extensions:

```@docs
export_vtk
```
