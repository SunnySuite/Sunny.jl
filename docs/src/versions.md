# Release Notes

## v0.8.1
(Dec 23, 2025)

* Fix plotting error with [`SCGA`](@ref)-calculated intensities ([#464](@ref)).
* Add [`copy_spins!`](@ref). Expose `jitter` parameter in
  [`minimize_energy!`](@ref) ([#465](@ref)).

## v0.8.0
(Nov 30, 2025)

**Breaking changes** in this release.

* Significant enhancements to crystal construction and symmetry analysis, which
  may **reorder certain atom and site indices** ([#421](@ref)). Crystal lattice
  vectors and positions are now idealized according to spacegroup data. In rare
  cases, this will cause re-indexing of atoms in the crystal. Site indexing
  conventions for certain reshaped systems may also change; use
  [`position_to_site`](@ref) for robust indexing of reshaped systems. Loading an
  mCIF as a chemical cell now employs a standardized Cartesian coordinate
  system, as would be obtained from [`lattice_vectors`](@ref). Attempting to
  perform symmetry analysis on a crystal with a larger-than-standard chemical
  cell will error rather than give a wrong result.
* When loading a CIF or mCIF, the precision parameter `symprec` becomes optional
  ([#413](@ref)).
* Various enhancements to [`minimize_energy!`](@ref). The return value becomes a
  struct that stores optimization statistics ([#430](@ref)). A small
  perturbation to the initial spin state breaks accidental symmetries
  ([#442](@ref)). Convergence to the local minimum becomes faster and more
  robust ([#453](@ref)).
* Fixes to [`load_nxs`](@ref) ([#420](@ref)).
* Add `interpolate` option to [`plot_intensities`](@ref). Selecting
  `interpolate=true` will significantly reduce file sizes of PDF exports
  containing 2D heatmap data ([#411](@ref)).
* [`set_spin_rescaling!`](@ref) now expects a scaling factor for each
  symmetry-distinct sublattice ([#444](@ref)).
* Introduce [`set_spin_s_at!`](@ref) to set the local quantum spin-``s``
  ([#454](@ref)).
* Add missing 3D support for [`q_space_grid`](@ref) ([#457](@ref)).
* If user-provided lattice vectors do not match crystallographic conventions,
  suggest the use of [`standardize`](@ref) ([#461](@ref)).
* Improve robustness of [`SpinWaveTheoryKPM`](@ref) ([#462](@ref)).

## v0.7.8
(Jul 1, 2025)

* Compatibility with Makie v0.24 ([#393](@ref)).

## v0.7.7
(Jun 25, 2025)

* Add [`set_spin_rescaling_for_static_sum_rule!`](@ref) which sets the classical
  dipole magnitude to ``\sqrt{s (s + 1)}`` for each quantum spin-``s`` moment.
* Add module [`SCGA`](@ref) for calculating [`intensities_static`](@ref) within
  the self-consistent Gaussian approximation ([#355](@ref)).
* Extend [`enable_dipole_dipole!`](@ref) to accept a demagnetization factor or
  tensor `demag`. The new default is isotropic demagnetization, `demag = 1/3`,
  appropriate for a spherical sample in vacuum. Set `demag = 0` to disable
  demagnetization ([#380](@ref)).
* Fix heatmaps in [`plot_intensities`](@ref) for very large grids
  ([#379](@ref)).
* Make energy minimization more reliable ([#397](@ref)).

## v0.7.6
(May 1, 2025)

* Extend [`powder_average`](@ref) to support static intensities.
* Vacancies defined by [`set_vacancy_at!`](@ref) are supported in linear spin
  wave theory. Empty sites are modeled using bosons that do not excite.
* The default implementation of [`SpinWaveTheoryKPM`](@ref) now uses Lanczos for
  higher accuracy.
* Fix correctness of [`suggest_magnetic_supercell`](@ref) when multiple
  wavevectors are provided.
* Fix atom indexing when setting interactions for a reshaped system
  ([#359](@ref)).
* Normalize `axis` argument to [`SpinWaveTheorySpiral`](@ref) for correctness.
* Fix thermal prefactor `kT` in spin wave theory ([#370](@ref)).
* In `:dipole` or `:SUN` mode, functions [`set_onsite_coupling!`](@ref) and
  [`set_pair_coupling!`](@ref) throw an error when used with quantum spin 1/2.
  Construct [`System`](@ref) using mode `:dipole_uncorrected` for such models
  ([#376](@ref)).

## v0.7.5
(Jan 20, 2025)

* [`view_crystal`](@ref) shows allowed quadratic anisotropy.
  [`print_site`](@ref) accepts an optional reference atom `i_ref`, with default
  of `i`. The optional reference bond `b_ref` of [`print_bond`](@ref) now
  defaults to `b`.
* The `regularization` parameter in [`SpinWaveTheory`](@ref) is reduced by half,
  and now corresponds to an effective energy shift. This may affect intensities,
  especially at small excitation energies.
* Update dependencies and in particular fix `color` compile error.

## v0.7.4
(Dec 6, 2024)

* Higher-precision convergence in [`minimize_energy!`](@ref).
* Make [`minimize_energy!`](@ref) compatible with [`set_vacancy_at!`](@ref).
* The [`System`](@ref) constructor now seeds its internal random number
  generator using Julia's task-local random number generator.
* Add [`print_irreducible_bz_paths`](@ref), which builds on
  [Brillouin.jl](https://github.com/thchr/Brillouin.jl) and
  [SeeK-path](http://www.materialscloud.org/tools/seekpath/).
* Add function [`view_bz`](@ref) for visualizing reciprocal-space objects in the
  context of the first Brillouin zone.
* Fix [`load_nxs`](@ref) for compatibility with recent JLD2.
* Fix Makie precompiles for faster time-to-first-plot in Julia 1.11
  ([#329](@ref)).

## v0.7.3
(Nov 12, 2024)

* Fix error in `print_symmetry_table` for slightly-distorted crystal cells
  ([#317](@ref)).
* Stabilize [`SpinWaveTheoryKPM`](@ref). It now automatically selects the
  polynomial order according to an error tolerance.
* Rename mode `:dipole_large_S` to `:dipole_uncorrected` to emphasize that
  corrections are missing.
* The [`Crystal`](@ref) constructor, by default, interprets a spacegroup number
  in its ITA standard setting, e.g., as used by the [Bilbao crystallographic
  server](https://www.cryst.ehu.es/cryst/get_wp.html). The keyword argument
  `setting` becomes `choice`, and can typically be omitted.
* Rename `primitive_cell_shape` to [`primitive_cell`](@ref).

## v0.7.2
(Sep 11, 2024)

* Fix error in `SampledCorrelations` with a coarse ``𝐪``-grid. ([#314](@ref)).
* Fix colorbar in `plot_intensities!` when all data is uniform ([#315](@ref)).
* An explicit `colorrange` can be used for plotting `intensities_bands`.

## v0.7.1
(Sep 3, 2024)

* Correctness fix for scalar biquadratic interactions specified with option
  `biquad` to [`set_exchange!`](@ref).
* Prototype implementation of entangled units.

## v0.7.0
(Aug 30, 2024)

**Breaking changes** in this release.

* The interface for calculating intensities has been revised to unify
  functionality across backends. The functions [`intensities_bands`](@ref),
  [`intensities`](@ref), and [`intensities_static`](@ref) no longer expect a
  "formula", and instead take keyword arguments directly. Pair correlations are
  now specified using [`ssf_perp`](@ref) and related functions. The constructors
  [`SampledCorrelations`](@ref) and [`SampledCorrelationsStatic`](@ref) replace
  `dynamic_correlations` and `static_correlations`, respectively.
* New function [`plot_intensities`](@ref) enables convenient plotting for many
  types of intensities plots. Mutating variant [`plot_intensities!`](@ref)
  enables multi-panel plots.
* One should now specify a range of ``𝐪``-points with [`q_space_path`](@ref) or
  [`q_space_grid`](@ref).
* [`SpinWaveTheorySpiral`](@ref) is available to perform calculations on
  generalized spiral structures, which may be incommensurate.
* [`repeat_periodically_as_spiral`](@ref) replaces
  `set_spiral_order_on_sublattice!` and `set_spiral_order!`.
* New convenience functions [`powder_average`](@ref) and
  [`domain_average`](@ref), which wrap [`intensities`](@ref).
* [`System`](@ref) now expects supercell dimensions as a `dims` keyword
  argument. [`Moment`](@ref) replaces `SpinInfo`. Lower-case `s` now labels
  quantum spin.
* In [`view_crystal`](@ref) and [`plot_spins`](@ref) use `ndims` instead of
  `dims` for the number of spatial dimensions.
* Binning features have been removed. Some functionality may be added back in a
  future release.
* Experimental `SpinWaveTheoryKPM` feature implements a [new
  algorithm](https://arxiv.org/abs/2312.08349) to enable intensities
  calculations at a computational cost that scales linearly in system size.


## v0.6.1
(Aug 2, 2024)

* **Breaking changes**: [`magnetic_moment`](@ref) is now reported in units of
  the Bohr magneton, ``μ_B``. For model systems where the Zeeman coupling aligns
  spin dipole with field (e.g., the Ising model convention), create a `SpinInfo`
  with `g=-1` ([#284](@ref)).
* More flexible [`Units`](@ref) system. `set_external_field!` is deprecated in
  favor of [`set_field!`](@ref), which now expects a field in energy units.
  [`enable_dipole_dipole!`](@ref) now expects a scale parameter ``μ_0 μ_B^2``
  that can be obtained from `units.vacuum_permeability`.

## v0.6.0
(Jun 18, 2024)

* Various correctness fixes. The magnetic moment is now anti-aligned with the
  spin dipole ([#190](@ref)), and the wavevector $𝐪$ in structure factor
  intensities $\mathcal{S}(𝐪,ω)$ now consistently represents momentum transfer
  _to_ the sample ([#270](@ref)). The new [Example 8](@ref "8. Momentum transfer
  conventions") demonstrates a model system where momentum transfers $±𝐪$ are
  inequivalent.
* Dynamical structure factor intensities now have a [precisely defined
  scale](@ref "Conventions for the Sunny-calculated structure factor"),
  independent of the calculator ([#264](@ref)). Consequently, color ranges in
  plots may need to be rescaled.
* [`Crystal`](@ref) can now infer a chemical unit cell from an mCIF file.
  `System` now supports [`set_dipoles_from_mcif!`](@ref). Through spglib, one
  can now [`standardize`](@ref) any `Crystal`, with an option to idealize site
  positions.

## v0.5.11
(Jun 2, 2024)

* Fixes for Makie 0.21.

## v0.5.10
(May 27, 2024)

* [`view_crystal`](@ref) called on a [`System`](@ref) now shows interactions,
  and optionally the spin or magnetic dipoles.
* Interactions for [`enable_dipole_dipole!`](@ref) are now supported in linear
  spin wave theory, with proper Ewald summation. For a faster alternative, the
  experimental function [`modify_exchange_with_truncated_dipole_dipole!`](@ref)
  will accept a real-space cutoff.
* Intensities calculated with `dynamic_correlations` now avoid "bleeding
  artifacts" at low-energy (long-timescale) modes. See [#246](@ref) for details.
  This eliminates the need for `process_trajectory=:symmetrize`.
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
  modes and [precisely specified](@ref "Structure Factor Conventions"). The
  g-tensor is applied by default (disable with `apply_g = false`). The intensity
  is additive with increasing number of magnetic ions in the chemical cell,
  consistent with SpinW ([#235](@ref)).
* Enhancements to [`view_crystal`](@ref). If a bond allows a DM interaction, its
  orientation will be shown visually. If a [`System`](@ref) argument is
  supplied, its exchange interactions will be shown.
* New function [`suggest_timestep`](@ref) to assist in performing accurate and
  efficient simulation of classical spin dynamics ([#149](@ref)).
* Scalar biquadratic interactions can again be set in `:dipole_large_S` mode via
  the keyword argument `biquad` of [`set_exchange!`](@ref).
* Significantly speed up `dynamic_correlations` for crystals with many atoms in
  the unit cell ([#204](@ref)).
* Renamings: `dt` replaces `Δt` and `damping` replaces `λ`. This affects
  [`Langevin`](@ref), [`ImplicitMidpoint`], and `dynamic_correlations`
  functions.

## v0.5.8
(Jan 4, 2024)

* Many bugs in the WGLMakie backend have become apparent, and are being tracked
  at [#211](@ref). Emit a warning if WGLMakie is detected, suggesting that
  GLMakie is preferred.
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

This release initiates some **major enhancements** to the user interface in
support of generalized SU(_N_) spin models. See [this documentation page](@ref
"Interaction Renormalization") for an illustration of the new features. Most
existing Sunny 0.5 models will continue to work with deprecation warnings, but
these will become hard errors Sunny v0.6.

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
  ferromagnet](https://github.com/SunnySuite/Sunny.jl/blob/main/examples/extra/heisenberg_animation.jl).
* Add [`set_spin_rescaling!`](@ref) feature, which supports improved spectral
  measurements at finite-$T$. This follows the method proposed in [Dahlbom et
  al., [arXiv:2310.19905]](https://arxiv.org/abs/2310.19905).

## v0.5.5
(Sep 29, 2023)

* [`reshape_supercell`](@ref) now allows reshaping to multiples of the primitive
  unit cell, which can speed up certain calculations. This is illustrated in the
  CoRh₂O₄ powder averaging tutorial.
* [`resize_supercell`](@ref) now allows all resizings.
* Added [`energy_per_site`](@ref).
* `set_spiral_order_on_sublattice!` cannot work on reshaped systems.
* Various bug fixes. In particular, an `intensity_formula` with `:full` will now
  uniformly calculate a `3x3` matrix of complex numbers.

## v0.5.4
(Sep 11, 2023)

* Various enhancements to [`view_crystal`](@ref). Atoms are now labeled by
  index, and bonds support interactive inspection (GLMakie only). Font sizes
  work correctly on Makie v0.20-beta. If using Makie v0.19 on a high-resolution
  display, pass `rescale=1.5` to enlarge font sizes.
* The function [`suggest_magnetic_supercell`](@ref) now requires only a list of
  wavevectors, and will return a $3×3$ matrix that can be programmatically
  passed to [`reshape_supercell`](@ref). The new tolerance parameter `tol`
  allows `suggest_magnetic_supercell` to approximate incommensurate wavevectors
  with nearby commensurate ones.
* New functions `set_spiral_order!` and `set_spiral_order_on_sublattice!` can be
  used to initialize a spiral, single-$Q$ order.
* Sunny now retains all constant energy shifts that have been introduced by
  anisotropy operators.
* Fix [`export_vtk`](@ref) functionality.

## v0.5.3
(Sep 8, 2023)

* Add `large_S_spin_operators` and `large_S_stevens_operators` to support
  single-ion anisotropies in dipole mode without renormalization. Set
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
experimental data. See `intensity_formula` for documentation. Use
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

Rename `DynamicStructureFactor` to `dynamic_correlations`. Similarly, replace
`InstantStructureFactor` with `instant_correlations`. The return type has been
renamed [`SampledCorrelations`](@ref) to emphasize that the object may be
holding thermodynamic samples, which are collected using [`add_sample!`](@ref).
Upon construction, the `SampledCorrelations` object will be empty (no initial
sample).

Remove `intensities` function. Instead, use one of `intensities_interpolated` or
`intensities_binned`. These will require an `intensity_formula`, which defines a
calculator (e.g., LSWT).

Rename `connected_path` to `reciprocal_space_path`, which now returns an
`xticks` object that can be used in plotting. Replace `spherical_shell` with
`reciprocal_space_shell` that functions similarly.

Rename `polarize_spin!` to [`set_dipole!`](@ref) for consistency with
[`set_coherent!`](@ref). The behavior of the former function is unchanged: the
spin at a given site will still be polarized along the provided direction.

Rename `all_sites` to `eachsite` consistent with Julia convention for iterators.

Rename `reshape_geometry` to [`reshape_supercell`](@ref), which is the
fundamental reshaping function. Rename `resize_periodically` to
[`resize_supercell`](@ref).

The constructor `SpinInfo` now requires a $g$-factor or tensor as a named
argument.

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

New exported functions [`global_position`](@ref), [`magnetic_moment`](@ref),
`all_sites`.

Remove all uses of
[`Base.deepcopy`](https://docs.julialang.org/en/v1/base/base/#Base.deepcopy)
which resolves crashes ([#65](@ref)).

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
System(crystal, dims, infos, mode)
```

The parameter `infos` is now a list of `SpinInfo` objects. Each defines spin
angular momentum $S = \frac{1}{2}, 1, \frac{3}{2}, …$, and an optional
$g$-factor or tensor.

The parameter `mode` is one of `:SUN` or `:dipole`.

### Setting interactions

Interactions are now added mutably to an existing `System` using the following
functions: `set_external_field!`, [`set_exchange!`](@ref),
[`set_onsite_coupling!`](@ref), [`enable_dipole_dipole!`](@ref).

As a convenience, one can use [`dmvec(D)`](@ref) to convert a DM vector to a
$3×3$ antisymmetric exchange matrix.

Fully general single-ion anisotropy is now possible. The function
[`set_onsite_coupling!`](@ref) expects the single ion anisotropy to be expressed
as a polynomial in symbolic spin operators `𝒮`, or as a linear combination of
symbolic Stevens operators `𝒪`. For example, an easy axis anisotropy in the
direction `n` may be written `D*(𝒮⋅n)^2`.

Stevens operators `𝒪[k,q]` admit polynomial expression in spin operators
`𝒮[α]`. Conversely, a polynomial of spin operators can be expressed as a linear
combination of Stevens operators. To see this expansion use
`print_anisotropy_as_stevens`.


### Inhomogeneous field

An external field can be applied to a single site with `set_external_field_at!`. 


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
