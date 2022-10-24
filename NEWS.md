# Sunny v0.4 development

## Breaking changes

**1. The interface for specifying anisotropy operators has changed.**

Anisotropy can now be expressed as a polynomial in spin operators `𝒮[α]`, or as
a linear combination of Stevens operators `𝒪[k,q]`. For example,
```julia
a1 = 20*(𝒮[1]^4 + 𝒮[2]^4 + 𝒮[3]^4)
a2 = 𝒪[4,0] + 5𝒪[4,4]
```

Stevens operators `𝒪[k,q]` admit polynomial expression in spin operators
`𝒮[α]`. Conversely, a polynomial of spin operators can be expressed as a linear
combination of Stevens operators. To see this expansion use
```julia
print_anisotropy_as_stevens(a1; N)
```
where `N = 2S+1` is the dimension of the spin operators. Alternatively, the
special value `N = 0` indicates the large-_S_ classical limit. Setting `N = 0`
prints `12X² + 𝒪₄₀ + 5𝒪₄₄`, where `X` is the spin magnitude squared. Observe
that `a1` and `a2` agree up to a constant shift.

The `anisotropy()` function takes a symbolic expression such as `a1` or `a2` and
produces an `Interaction`, which can be used in either dipole-only mode or
SU(_N_) mode. For example, to specify an easy-axis in the `n` direction with
magnitude `D`, one may use:
```julia
anisotropy(-D*(𝒮⋅n)^2, site_index; label)
```

Another convenient syntax is `𝒮'*J*𝒮` to produce a general quadratic
interaction with matrix-elements `J`.

**2. The function `print_symmetry_table()` replaces `print_bond_table()`.**

This new function describes both symmetry-allowed anisotropies in addition to
symmetry-allowed interactions on bonds.

**3. When reading CIF files, the field `_atom_site_label` is now used in place of the field `_atom_site_type_symbol`.**

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to `subcrystal()`) will need to be updated to use the new label.`
