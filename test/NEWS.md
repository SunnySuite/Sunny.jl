# Sunny v0.4 development

## Breaking changes

**1. The interface for specifying anisotropy operators has changed.**

Anisotropy operators can now be specified as either a polynomial in spin
operators, or a linear combination of Stevens operators. For example:
```julia
Sx, Sy, Sz = spin_operators
a1 = Sx^4 + Sy^4 + Sz^4

𝒪₄ = stevens_operators[4]
a2 = 𝒪₄[0] + 5𝒪₄[4]
```

In the classical limit, spin operators are replaced with expectation values.
Here, Stevens operators retain only the leading order terms in powers of _S_ and
become homogeneous polynomials, e.g.
```julia
print_operator_as_classical_polynomial(a2) 
# Output: 8sx⁴ - 24sx²sy² - 24sx²sz² + 8sy⁴ - 24sy²sz² + 8sz⁴
```

In this example, `a2` corresponds to `a1` up to a rescaling and irrelevant shift
because `(sx²+sy²+sz²)²` is a constant.

To get an `Interaction`, use, e.g., `anisotropy(a1, site_index; label)`. This
interaction can be used in either dipole-only mode or SU(_N_) mode.

**2. When reading CIF files, the field `_atom_site_label` is now used in place of the field `_atom_site_type_symbol`**

This is required for correctness. The field `_atom_site_label` is guaranteed to
be present, and is guaranteed to be a distinct label for each
symmetry-inequivalent site. Code that explicitly referred to site labels (e.g.
in calls to `subcrystal()`) will need to be updated to use the new label.`