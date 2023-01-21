# Sunny v0.4 development

This is a big update with many breaking changes. The example FeI2 notebook
(_**Reference**_) illustrates many of the changes.


## Creating a spin `System`

`SpinSystem` has been renamed `System`. Its constructor is,

```julia
sys = System(crystal, latsize, infos; mode)
```
The parameter `infos` is now a list of `SpinInfo`
objects,

```julia
SpinInfo(site::Int, S; g=2),
```

involving spin angular momentum $S = (1/2, 1, 3/2, …)$ and an optional
``g``-factor or tensor.

`System` now requires an additional parameter `mode` which must be one of
`:dipole`, `:SUN`, or `:projected`. 

## Setting interactions

Interactions are now added mutably to an existing System. The previous
intraction list will be replaced with the following functions:

```julia
set_external_field!(sys, h)
set_exchange!(sys, J, i)
set_exchange_with_biquadratic!(sys, J_quad, J_biquad, i)
set_anisotropy!(sys, Λ, i)
enable_dipole_dipole!(sys)
disable_dipole_diople!(sys)
```

As a convenience, the new function `dmvec(D)` converts a DM vector an
antisymmetric exchange matrix.

## Anisotropy operators

The function `set_anisotropy!(sys, Λ, i)` expects the single ion anisotropy `Λ`
to be expressed as a polynomial in spin operators `𝒮[α]`, or as a linear
combination of Stevens operators `𝒪[k,q]`. For example,
```julia
Λ1 = 20*(𝒮[1]^4 + 𝒮[2]^4 + 𝒮[3]^4)
Λ2 = 𝒪[4,0] + 5𝒪[4,4]
```

Stevens operators `𝒪[k,q]` admit polynomial expression in spin operators
`𝒮[α]`. Conversely, a polynomial of spin operators can be expressed as a linear
combination of Stevens operators. To see this expansion use
```julia
print_anisotropy_as_stevens(Λ1; N)
```
where `N = 2S+1` is the dimension of the spin operators. Alternatively, the
special value `N = 0` indicates the large-_S_ classical limit. In this case, `N
= 0` prints `12X² + 𝒪₄₀ + 5𝒪₄₄`, where `X` is the spin magnitude squared.
Observe that `Λ1` and `Λ2` agree up to a constant shift.

An easy-axis anisotropy may be expressed as `-D*(𝒮⋅n)^2`. A general quadratic
anisotropy with matrix elements `A` may be written `𝒮'*A*𝒮`.

## Inhomogeneous interactions (Planned)

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


## Structure factor


_**David: Please describe the rest of the new API**_

```julia
add_trajectory!(sf, sys)
```

This will run a dynamical trajectory of a copy of the system `sys`, and accumulate data into `sf`. Allocations are avoided by using buffer space in `sys`.


## Various

* `print_symmetry_table()` replaces `print_bond_table()`.

The new function includes the list of symmetry-allowed single ion anisotropies
in addition to exchange interactions.


* When reading CIF files, the field `_atom_site_label` is now used in place of the field `_atom_site_type_symbol`.

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to `subcrystal()`) will need to be updated to use the new label.
