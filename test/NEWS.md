# Sunny v0.4 development

## Breaking changes

**1. The interface for specifying anisotropy operators has changed.**

Anisotropy operators can now be specified as either a polynomial in spin
operators `𝒮` or a linear combination of Stevens operators `𝒪`. For example:
```julia
a1 = 20*(𝒮[1]^4 + 𝒮[2]^4 + 𝒮[3]^4)
a2 = 𝒪[4,0] + 5𝒪[4,4]
```

In the classical limit, spin operators are replaced with expectation values. In
this limit, Stevens operators retain only the leading order terms in powers of
_S_ and become homogeneous polynomials. 
```julia
print_anisotropy_as_spins(a2) 
# Output: 8𝒮₁⁴ - 24𝒮₁²𝒮₂² - 24𝒮₁²𝒮₃² + 8𝒮₂⁴ - 24𝒮₂²𝒮₃² + 8𝒮₃⁴
```

Conversely, given a classical spin polynomial, Sunny can print the corresponding
expansion in Stevens operators,
```julia
print_anisotropy_as_stevens(a1)
# Output: 12X² + 𝒪₄₀ + 5𝒪₄₄
```

In this case, `a1` and `a2` are the same up to an irrelevant shift. The symbol
`X` indicates spin magnitude squared.

The `anisotropy()` function takes these operators and produces an `Interaction`,
which can be used in either dipole-only mode or SU(_N_) mode. For example, to
specify an easy-axis in the `n` direction with magnitude `D`, one may use:
```julia
anisotropy((-D*(𝒮⋅n)^2, site_index; label)
```

Another convenient syntax is `𝒮'*J*𝒮` to produce a general quadratic
interaction with matrix-elements `J`.

**2. When reading CIF files, the field `_atom_site_label` is now used in place of the field `_atom_site_type_symbol`**

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to `subcrystal()`) will need to be updated to use the new label.`