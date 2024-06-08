# Version History

## v0.6.0
(In progress)

* Various correctness fixes. The magnetic moment is now anti-aligned with the
  spin dipole ([Issue 190](https://github.com/SunnySuite/Sunny.jl/issues/190)),
  and the wavevector $𝐪$ in structure factor intensities $\mathcal{S}(𝐪,ω)$
  now consistently represents momentum transfer _to_ the sample ([Issue
  270](https://github.com/SunnySuite/Sunny.jl/issues/270)). The new [Example
  8](@ref "8. Momentum transfer conventions") demonstrates a model system where
  momentum transfers $±𝐪$ are inequivalent.
* Dynamical structure factor intensities now have a [precisely defined
  scale](@ref "Conventions for the Sunny-calculated structure factor"),
  independent of the calculator ([Issue
  264](https://github.com/SunnySuite/Sunny.jl/issues/264)). Consequently, color
  ranges in plots may need to be rescaled.

## v0.5.11
(June 2, 2024)

* Fixes for Makie 0.21.

## v0.5.10
(May 27, 2024)

* [`view_crystal`](@ref) called on a [`System`](@ref) now shows interactions,
  and optionally the spin or magnetic dipoles.
* Interactions for [`enable_dipole_dipole!`](@ref) are now supported in linear
  spin wave theory, with proper Ewald summation. For a faster alternative, the
  experimental function [`modify_exchange_with_truncated_dipole_dipole!`](@ref)
  will accept a real-space cutoff.
* Intensities calculated with [`dynamical_correlations`](@ref) now avoid
  "bleeding artifacts" at low-energy (long-timescale) modes. See [PR
  246](https://github.com/SunnySuite/Sunny.jl/pull/246) for details. This
  eliminates the need for `process_trajectory=:symmetrize`.
* When passed to `intensity_formula`, the special value `zero(FormFactor)` can
  now be used to disable contributions from a given site. For an example, see
  the ported [SpinW tutorial 19](@ref "SW19 - Different magnetic ions").
* Broadening kernels [`gaussian`](@ref) and [`lorentzian`](@ref) now expect a
  full width at half maximum (`fwhm`) keyword argument.
* Experimental support for calculations on generalized spiral phases. For an
  example, see the ported [SpinW tutorial 18](@ref "SW18 - Distorted kagome").
* Correctness fix for the case where spin-$S$ varies between sites in
  dipole-mode. In SU($N$) mode, however, there is still no support for varying
  the Hilbert space dimension $N$ between sites.
* Correctness fix in long-range dipole-dipole interactions for systems with
  multiple cells.
* Correctness fix in general biquadratic interactions (beyond scalar) for spin
  wave theory in dipole-mode.
* Correctness fix for reading Mantid `.nxs` files.


## v0.5.9
(Mar 25, 2024)

* **Correctness fixes**: Structure factor conventions are now uniform across
  modes and [precisely specified](@ref "Structure Factor Calculations"). The
  g-tensor is applied by default (disable with `apply_g = false`). The intensity
  is additive with increasing number of magnetic ions in the chemical cell,
  consistent with SpinW. [Issue
  #235](https://github.com/SunnySuite/Sunny.jl/issues/235).
* Enhancements to [`view_crystal`](@ref). If a bond allows a DM interaction, its
  orientation will be shown visually. If a [`System`](@ref) argument is
  supplied, its exchange interactions will be shown..
* New function [`suggest_timestep`](@ref) to assist in performing accurate and
  efficient simulation of classical spin dynamics. [Issue
  #149](https://github.com/SunnySuite/Sunny.jl/issues/149).
* Scalar biquadratic interactions can again be set in `:dipole_large_S` mode via
  the keyword argument `biquad` of [`set_exchange!`](@ref).
* Significantly speed up [`dynamical_correlations`](@ref) for crystals with many
  atoms in the unit cell. [Issue
  #204](https://github.com/SunnySuite/Sunny.jl/issues/204).
* Renamings: `dt` replaces `Δt` and `damping` replaces `λ`. This affects
  [`Langevin`](@ref), [`ImplicitMidpoint`], and [`dynamical_correlations`](@ref)
  functions.

## v0.5.8
(Jan 4, 2024)

* Many bugs in the WGLMakie backend have become apparent, and are being tracked
  at [Issue #211](https://github.com/SunnySuite/Sunny.jl/issues/211). Emit a
  warning if WGLMakie is detected, suggesting that GLMakie is preferred.
* Various improvements to [`view_crystal`](@ref). A distance parameter is no
  longer expected. Cartesian axes now appear as "compass" in bottom-left. Custom
  list of reference bonds can be passed. Toggle to view non-magnetic atoms in
  root crystal. Atoms now colored using [CPK/JMol
  conventions](https://en.wikipedia.org/wiki/CPK_coloring).

## v0.5.7
(Nov 26, 2023)

* Update form factor coefficients, which now include `Mn5`.
* Fix [`merge_correlations`](@ref) and the [Parallelizing Calculations](@ref)
  tutorial.
* Remove internal functions `*_primitive_crystal`. Instead, it is recommended to
  use the conventional unit cell, and later call [`reshape_supercell`](@ref).
* Require Makie 0.20. An important new feature is resolution-independent scaling
  of font sizes. New figures expect `size` instead of `resolution`, and no
  longer accept `rescale`.

## v0.5.6
(Nov 8, 2023)

This release initiates some **major enhancements** to the user interface in support
of generalized SU(_N_) spin models. See [this documentation
page](https://sunnysuite.github.io/Sunny.jl/dev/renormalization.html) for an
illustration of the new features. Most existing Sunny 0.5 models will continue
to work with deprecation warnings, but these will become hard errors Sunny v0.6.

* General pair couplings are now supported in [`set_pair_coupling!`](@ref) and
  [`set_pair_coupling_at!`](@ref). `:SUN` mode supports interactions of any
  order, but `:dipole` mode is limited to bilinear and biquadratic coupling of
  the spin.
* To perform a calculation with dipoles in the large-$S$ limit, use the new mode
  `:dipole_large_S` when constructing a [`System`](@ref).
* Deprecate the option `biquad` to [`set_exchange!`](@ref). Use instead
  `set_pair_coupling!`, which generalizes beyond the scalar biquadratic.
* Deprecate `spin_operators`, `stevens_operators`, `large_S_spin_operators` and
  `large_S_stevens_operators`. Use instead [`spin_matrices`](@ref) and
  [`stevens_matrices`](@ref), which require a specific spin-``S`` label. To
  infer this, one can use [`spin_label`](@ref).
* Remove unused option `energy_tol` in [`SpinWaveTheory`](@ref).
* Animated spin dynamics is now possible. Call `notify` on the result of
  [`plot_spins`](@ref) to trigger redrawing of the frame. The argument `colorfn`
  to `plot_spins` supports animation of colors. See [example usage for a
  Heisenberg
  ferromagnetic.](https://github.com/SunnySuite/Sunny.jl/blob/main/examples/extra/heisenberg_animation.jl)
* Add [`set_spin_rescaling!`](@ref) feature, which supports improved spectral
  measurements at finite-$T$. This follows the method proposed in [Dahlbom et
  al., [arXiv:2310.19905]](https://arxiv.org/abs/2310.19905).

## v0.5.5
(Sept 29, 2023)

* [`reshape_supercell`](@ref) now allows reshaping to multiples of the primitive
  unit cell, which can speed up certain calculations. This is illustrated in the
  CoRh₂O₄ powder averaging tutorial.
* [`resize_supercell`](@ref) now allows all resizings.
* Added [`energy_per_site`](@ref).
* [`set_spiral_order_on_sublattice!`](@ref) cannot work on reshaped systems.
* Various bug fixes. In particular, an `intensity_formula` with `:full` will now
  uniformly calculate a `3x3` matrix of complex numbers.

## v0.5.4
(Sept 11, 2023)

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

## v0.5.3
(Sept 8, 2023)

* Add `large_S_spin_operators` and `large_S_stevens_operators`
  to support single-ion anisotropies in dipole mode without renormalization. Set
  `large_S=true` in [`set_exchange!`](@ref) to avoid renormalization of
  biquadratics.
* [`view_crystal`](@ref) has been rewritten in Makie.
* [`plot_spins`](@ref) now expects `ghost_radius` in physical length units.
* [`SpinWaveTheory`](@ref) will (currently) error if provided a system with
  [`enable_dipole_dipole!`](@ref).

## v0.5.2
(Aug 30, 2023)

* Form factors for 5d transition ions.
* Download links for notebooks and scripts on each doc example
* Various bug fixes.

## v0.5.1
(Aug 23, 2023)

* Fix binning edge cases.
* [`plot_spins`](@ref) accepts resolution argument.

## v0.5.0
(Aug 21, 2023)

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

Rename `DynamicStructureFactor` to [`dynamical_correlations`](@ref). Similarly,
replace `InstantStructureFactor` with [`instant_correlations`](@ref). The return
type has been renamed [`SampledCorrelations`](@ref) to emphasize that the object
may be holding thermodynamic samples, which are collected using
[`add_sample!`](@ref). Upon construction, the `SampledCorrelations` object will
be empty (no initial sample).

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

## v0.4.3
(Jun 23, 2023)

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


## v0.4.2
(Feb 27, 2023)

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

## v0.4.1
(Feb 13, 2023)

The function [`to_inhomogeneous`](@ref) creates a system that supports
inhomogeneous interactions, which can be set using [`set_exchange_at!`](@ref),
etc.

`set_biquadratic!` replaces `set_exchange_with_biquadratic!`.


## v0.4.0
(Feb 10, 2023)

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
interface see the documentation tutorials.


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
