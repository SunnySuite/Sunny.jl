# Sunny v0.4 development

## Breaking changes

**1. The interface for specifying anisotropy operators has changed.**

Anisotropy operators can now be specified as either a polynomial in spin operators, or a linear combination of Stevens operators. For example:
```julia
Sx, Sy, Sz = spin_operators
a1 = Sx^4 + Sy^4 + Sz^4

𝒪₄ = stevens_operators[4]
a2 = 𝒪₄[0] + 5𝒪₄[4]
```

To get the spin-_S_ matrix representation in the eigenbasis of ``\hat{S}_z``, use
```julia
m1 = operator_as_matrix(a1; S=3/2)
m2 = operator_as_matrix(a2; S=3/2)
```

In this example, notice
```julia
20m1 ≈ m2 + 408I
```
Up to a rescaling, both `a1` and `a2` represent the unique quartic anisotropy that is symmetry-allowed for cubic point groups.

To get an `Interaction`, use, e.g., `anisotropy(a1, site_index; label)`. This interaction can be used in either dipole-only mode or SU(_N_) mode.

**2. When reading CIF files, the field `_atom_site_label` is now used in place of the field `_atom_site_type_symbol`**

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to `subcrystal()`) will need to be updated to use the new label.`