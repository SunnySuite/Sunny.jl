# Version 0.4

This update includes many breaking changes.

### Creating a spin `System`

`SpinSystem` has been renamed [`System`](@ref). Its constructor now has the form,

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
[`set_exchange_with_biquadratic!`](@ref), [`set_anisotropy!`](@ref),
[`enable_dipole_dipole!`](@ref).

As a convenience, one can use [`dmvec(D)`](@ref) to convert a DM vector to a
$3×3$ antisymmetric exchange matrix.

Fully general single-ion anisotropy is now possible. The function
[`set_anisotropy!`](@ref) expects the single ion anisotropy to be expressed as a
polynomial in symbolic spin operators [`𝒮`](@ref), or as a linear combination
of symbolic Stevens operators [`𝒪`](@ref). For example, an easy axis anisotropy
in the direction `n` may be written `D*(𝒮⋅n)^2`.

Stevens operators `𝒪[k,q]` admit polynomial expression in spin operators
`𝒮[α]`. Conversely, a polynomial of spin operators can be expressed as a linear
combination of Stevens operators. To see this expansion use
[`print_anisotropy_as_stevens`](@ref).


### Inhomogeneous interactions

Vacancies can be introduced with [`set_vacancy_at!`](@ref).

An external field can be applied to a single site with
[`set_external_field_at!`](@ref). 

Support for inhomogeneous exchange interactions is planned. One will first call
`enable_inhomogeneous_exchange!`, allowing subsequent use of `set_exchange_at!`.


### Structure factor rewrite

The calculation of structure factors has been completely rewritten. For the new
interface, see the [Structure Factor Calculations](@ref) page.


### Various

* The "Sampler" interface is in flux. [`Langevin`](@ref) replaces both
  `LangevinHeunP` and `LangevinSampler`. Local spin-flip Monte Carlo sampling
  methods are temporarily broken.

* [`repeat_periodically`](@ref) replaces `extend_periodically`.

Additional related functions include [`resize_periodically`](@ref) and
[`reshape_geometry`](@ref), the latter being fundamental.

* [`print_symmetry_table`](@ref) replaces `print_bond_table()`.

The new function includes the list of symmetry-allowed single ion anisotropies
in addition to exchange interactions.


* When reading CIF files, the field `_atom_site_label` is now used in place of
  the field `_atom_site_type_symbol`.

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to [`subcrystal`](@ref)) will need to be updated to use the new label.
