# Version History

## v0.5.6

* General pair couplings are now supported in [`set_pair_coupling!`](@ref) and
  [`set_pair_coupling_at!`](@ref). `:SUN` supports interactions of any order,
  but `:dipole` mode is limited to bilinear and biquadratic coupling of the
  spin.
* To perform a calculation with dipoles in the large-$S$ limit, use the new mode
  `:dipole_large_S` when constructing a [`System`](@ref).
* Deprecate the option `biquad` to [`set_exchange!`](@ref). Use instead
  `set_pair_coupling!`, which generalizes beyond the scalar biquadratic.
* Deprecate `spin_operators`, `stevens_operators`, `large_S_spin_operators` and
  `large_S_stevens_operators`. Use instead [`spin_matrices`](@ref) and
  [`stevens_matrices`](@ref), which require a specific spin-``S`` label. To
  infer this, one can use [`spin_label`](@ref).
* Remove unused option `energy_tol` in [`SpinWaveTheory`](@ref).

## v0.5.5

* [`reshape_supercell`](@ref) now allows reshaping to multiples of the primitive
  unit cell, which can speed up certain calculations. This is illustrated in the
  CoRh₂O₄ powder averaging tutorial.
* [`resize_supercell`](@ref) now allows all resizings.
* Added [`energy_per_site`](@ref).
* [`set_spiral_order_on_sublattice!`](@ref) cannot work on reshaped systems.
* Various bug fixes. In particular, an `intensity_formula` with `:full` will now
  uniformly calculate a `3x3` matrix of complex numbers.

# v0.5.4

* Various enhancements to [`view_crystal`](@ref). Atoms are now labeled by
  index, and bonds support interactive inspection (GLMakie only). Font sizes
  work correctly on Makie v0.20-beta. If using Makie v0.19 on a high-resolution
  display, pass `rescale=1.5` to enlarge font sizes.
* The function [`suggest_magnetic_supercell`](@ref) now requires only a list of
  wavevectors, and will return a $3×3$ matrix that can be programmatically
  passed to [`reshape_supercell`](@ref). The new tolerance parameter `tol`
  allows `suggest_magnetic_supercell` to approximate incommensurate wavevectors
  with nearby commensurate ones.
* New functions [`set_spiral_order!`](@ref) and
  [`set_spiral_order_on_sublattice!`](@ref) can be used to initialize a spiral,
  single-$Q$ order.
* Sunny now retains all constant energy shifts that have been introduced by
  anisotropy operators.
* Fix `export_vtk` functionality.

# v0.5.3

* Add `large_S_spin_operators` and `large_S_stevens_operators`
  to support single-ion anisotropies in dipole mode without renormalization. Set
  `large_S=true` in [`set_exchange!`](@ref) to avoid renormalization of
  biquadratics.
* [`view_crystal`](@ref) has been rewritten in Makie.
* [`plot_spins`](@ref) now expects `ghost_radius` in physical length units.
* [`SpinWaveTheory`](@ref) will (currently) error if provided a system with
  [`enable_dipole_dipole!`](@ref).

# v0.5.2

* Form factors for 5d transition ions.
* Download links for notebooks and scripts on each doc example
* Various bug fixes.

# v0.5.1

* Fix binning edge cases.
* [`plot_spins`](@ref) accepts resolution argument.

# v0.5.0

**New features**.

Support for Linear Spin Wave Theory in `:dipole` and `:SUN` modes. (Thanks Hao
Zhang!)

New function [`minimize_energy!`](@ref) to efficiently find an optimal
configuration of spin dipoles or SU(_N_) coherent states.

Major refactors and enhancements to intensity calculations. This new interface
allows unification between LSWT and classical spin dynamics calculations. This
interface allows: Custom observables as local quantum operators, better support
for linebroadening, and automatic binning to facilitate comparison with
experimental data. See [`intensity_formula`](@ref) for documentation. Use
[`load_nxs`](@ref) to load experimental neutron scattering data.

**Breaking changes**.

Require Julia 1.9.

Replace `set_anisotropy!` with a new function [`set_onsite_coupling!`](@ref)
(and similarly [`set_onsite_coupling_at!`](@ref)). The latter expects an
explicit matrix representation for the local Hamiltonian. This can be
constructed, e.g., as a linear combination of `stevens_operators`, or as a
polynomial of `spin_operators`. To understand the mapping between these two, the
new function [`print_stevens_expansion`](@ref) acts on an arbitrary local
operator.

Remove `set_biquadratic!`. Instead, use an optional keyword argument `biquad` to
[`set_exchange!`](@ref).

Rename `DynamicStructureFactor` to [`dynamical_correlations`](@ref).
Similarly, replace `InstantStructureFactor` with [`instant_correlations`](@ref).
The return type has been renamed [`SampledCorrelations`](@ref) to emphasize that
the object may be holding thermodynamic samples, which are collected using
[`add_sample!`](@ref).

Remove `intensities` function. Instead, use one of
[`intensities_interpolated`](@ref) or [`intensities_binned`](@ref). These will
require an [`intensity_formula`](@ref), which defines a calculator (e.g., LSWT).

Rename `connected_path` to [`reciprocal_space_path`](@ref), which now returns an
`xticks` object that can be used in plotting. Replace `spherical_shell` with
[`reciprocal_space_shell`](@ref) that functions similarly.

Rename `polarize_spin!` to [`set_dipole!`](@ref) for consistency with
[`set_coherent!`](@ref). The behavior of the former function is unchanged: the
spin at a given site will still be polarized along the provided direction.

Rename `all_sites` to `eachsite` consistent with Julia convention for iterators.

Rename `reshape_geometry` to [`reshape_supercell`](@ref), which is the
fundamental reshaping function. Rename `resize_periodically` to
[`resize_supercell`](@ref).

The constructor [`SpinInfo`](@ref) now requires a $g$-factor or tensor as a
named argument.

The constructor [`FormFactor`](@ref) no longer accepts an atom index. Instead,
the form factors are associated with site-symmetry classes in order of
appearance.

# v0.4.3

**Experimental** support for linear [`SpinWaveTheory`](@ref), implemented in
SU(_N_) mode. This module may evolve rapidly.

Implement renormalization of single-ion anisotropy and biquadratic interactions
when in `:dipole` mode. This makes the model more faithful to the quantum
mechanical Hamiltonian, but is also a **breaking change**.

Various improvements and bugfixes for [`to_inhomogeneous`](@ref). Setting
inhomogeneous interactions via [`set_exchange_at!`](@ref) should now infer the
correct bond offset direction, or will report an ambiguity error. Ambiguities
can be resolved by passing an explicit `offset`.

The function [`remove_periodicity!`](@ref) disables periodicity along specified
dimensions.

Rename `StaticStructureFactor` to `InstantStructureFactor`.


# v0.4.2

Introduce [`LocalSampler`](@ref), a framework for MCMC sampling with local spin
updates.

Rename `print_dominant_wavevectors` to [`print_wrapped_intensities`](@ref) to
reduce confusion with the physical instantaneous intensities.

The function `spherical_shell` now takes a radius in physical units of inverse
Å.

New exported functions [`global_position`](@ref), [`magnetic_moment`](@ref), `all_sites`.

Remove all uses of
[`Base.deepcopy`](https://docs.julialang.org/en/v1/base/base/#Base.deepcopy)
which [resolves crashes](https://github.com/SunnySuite/Sunny.jl/issues/65).

# v0.4.1

The function [`to_inhomogeneous`](@ref) creates a system that supports
inhomogeneous interactions, which can be set using [`set_exchange_at!`](@ref),
etc.

`set_biquadratic!` replaces `set_exchange_with_biquadratic!`.


# v0.4.0

This update includes many breaking changes, and is missing some features of
0.3.0.

### Creating a spin `System`

Rename `SpinSystem` to [`System`](@ref). Its constructor now has the form,

```julia
System(crystal, latsize, infos, mode)
```

The parameter `infos` is now a list of [`SpinInfo`](@ref) objects. Each defines
spin angular momentum $S = \frac{1}{2}, 1, \frac{3}{2}, …$, and an optional
$g$-factor or tensor.

The parameter `mode` is one of `:SUN` or `:dipole`.

### Setting interactions

Interactions are now added mutably to an existing `System` using the following
functions: [`set_external_field!`](@ref), [`set_exchange!`](@ref),
[`set_onsite_coupling!`](@ref), [`enable_dipole_dipole!`](@ref).

As a convenience, one can use [`dmvec(D)`](@ref) to convert a DM vector to a
$3×3$ antisymmetric exchange matrix.

Fully general single-ion anisotropy is now possible. The function
[`set_onsite_coupling!`](@ref) expects the single ion anisotropy to be expressed as a
polynomial in symbolic spin operators `𝒮`, or as a linear combination
of symbolic Stevens operators `𝒪`. For example, an easy axis anisotropy
in the direction `n` may be written `D*(𝒮⋅n)^2`.

Stevens operators `𝒪[k,q]` admit polynomial expression in spin operators
`𝒮[α]`. Conversely, a polynomial of spin operators can be expressed as a linear
combination of Stevens operators. To see this expansion use
`print_anisotropy_as_stevens`.


### Inhomogeneous field

An external field can be applied to a single site with
[`set_external_field_at!`](@ref). 


### Structure factor rewrite

The calculation of structure factors has been completely rewritten. For the new
interface, see the [FeI₂ at Finite Temperature](@ref) page.


### Various

* The "Sampler" interface is in flux. [`Langevin`](@ref) replaces both
  `LangevinHeunP` and `LangevinSampler`. Local spin-flip Monte Carlo sampling
  methods are temporarily broken.

* [`repeat_periodically`](@ref) replaces `extend_periodically`.

Additional related functions include `resize_periodically` and
`reshape_geometry`, the latter being fundamental.

* [`print_symmetry_table`](@ref) replaces `print_bond_table()`.

The new function includes the list of symmetry-allowed single ion anisotropies
in addition to exchange interactions.

* When reading CIF files, the field `_atom_site_label` is now used in place of
  the field `_atom_site_type_symbol`.

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to [`subcrystal`](@ref)) will need to be updated to use the new label.
