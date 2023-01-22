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

The parameter `mode` is one of `:dipole`, `:SUN`, or `:projected`. 

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

### Inhomogeneous interactions (Planned)

Spatially inhomogeneous interactions can be get or set using the following methods:

```julia
set_vacancy_at!(sys, idx)

set_external_field_at!(sys, h, idx)
get_external_field_at!(sys, idx)

enable_inhomogeneous_exchange!(sys) # Once enabled, cannot be disabled

set_exchange_at!(sys, J, idx)
get_exchange_at(sys, idx)
```

The parameter `idx` has the shape `(n1, n2, n3, atom)`, where `(n1,n2,n3)`
labels a unit cell, and `atom` is an index within this unit cell.


### Structure factor rewrite

The calculation of structure factors has been completely rewritten. For the new
interface, see the [Structure Factor Calculations](@ref) page.


### Various

* `print_symmetry_table()` replaces `print_bond_table()`.

The new function includes the list of symmetry-allowed single ion anisotropies
in addition to exchange interactions.


* When reading CIF files, the field `_atom_site_label` is now used in place of the field `_atom_site_type_symbol`.

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to `subcrystal()`) will need to be updated to use the new label.
