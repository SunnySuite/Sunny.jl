# Single-ion anisotropy

A unique feature of Sunny is its support for building classical models where
each quantum spin is represented as a full $N$-level system, rather than just an
expected dipole. This generalized formalism is especially useful for modeling,
e.g., quantum spin Hamiltonians that include a strong effective
single-anisotropy.

## Defining on-site couplings

A single-ion anisotropy can be specified either as a polynomial of spin
operators, or as a linear combination of Stevens operators. A quantum spin-$S$
has $N = 2S + 1$ levels, and local quantum operators can be fully described as
$N×N$ matrices. To get the spin operators $\hat{S}^{\{x,y,z\}}$ as matrices, use
[`spin_operators`](@ref). For example, in the case of spin-$1/2$, this function
returns Pauli matrices. The Stevens operators $\hat{\mathcal{O}}_{k,q}$ are
defined as polynomials of the $\hat{S}^α$, and are accessible using
[`stevens_operators`](@ref).

Although Sunny encourages representing an explicit $N×N$ matrix representation
of local operators, it is also possible work in a "large-$S$ limit", whereby
each spin operator is replaced with its expectation value. To build such a model
in Sunny, use [`large_S_spin_operators`](@ref) and
[`large_S_stevens_operators`](@ref). Each of these will provide operators that
are symbolic polynomials in the commuting spin expectation values.

After constructing the single-ion anisotropy, either as a matrix or as a
symbolic polynomial, use [`set_onsite_coupling!`](@ref) to add it to a system. To
see the expansion of any local operators in terms of Stevens operators, use
[`print_stevens_expansion`](@ref). 

## A note about renormalization

When constructing a [`System`](@ref) in Sunny, one must select from two possible
modes. The mode `:SUN` models each spin as $N$ complex amplitudes, and is the
most variationally accurate. The mode `:dipole` constrains the spin dynamics to
the space of pure dipoles. In either mode, we generally recommend defining the
anisotropy operator as an $N×N$ matrix. In `:dipole` mode, Sunny will
automatically perform a [renormalization of
anisotropy](https://arxiv.org/abs/2304.03874) to achieve maximal consistency
with `:SUN` mode.

Anisotropy renormalization is not desirable if one is attempting to reproduce a
previous study, or to model a micromagnetic systems for which $S \to \infty$ is
physically valid. To _avoid_ renormalization, simply construct the anisotropy
operators using [`large_S_spin_operators`](@ref) in place of `spin_operators`.

Note that Sunny will also perform renormalization on the scalar biquadratic
exchange interactions. See [`set_exchange!`](@ref) for more details.

## Stevens operators

The Stevens operators $\hat{\mathcal{O}}_{k,q}$ are defined as polynomials of
angular momentum operators $\hat{S}_{\{x,y,z\}}$ in some spin-$S$ representation.

Using

```math
\begin{align*}
X &= \mathbf{\hat{S}} \cdot \mathbf{\hat{S}} = S (S+1) \\
\hat{S}_\pm &= \hat{S}_x \pm i \hat{S}_y \\
\phi_+ &= \frac{1}{4},\quad \phi_- = \frac{1}{4 i},
\end{align*}
```

some Stevens operators are defined as,

```math
\begin{align*}
\hat{\mathcal{O}}_{2,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm \hat{S}_{-}^{2})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{2,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm \hat{S}_{-})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{2,0} & =3\hat{S}_{z}^{2}-X\\
\\
\hat{\mathcal{O}}_{4,\pm4} & =\phi_{\pm}(\hat{S}_{+}^{4}\pm \hat{S}_{-}^{4})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm \hat{S}_{-}^{3})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm \hat{S}_{-}^{2})(7\hat{S}_{z}^{2}-(X+5))+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm \hat{S}_{-})(7\hat{S}_{z}^{3}-(3X+1)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,0} & =35\hat{S}_{z}^{4}-(30X-25)\hat{S}_{z}^{2}+(3X^{2}-6X)+\mathrm{h.c.}\\
\\
\hat{\mathcal{O}}_{6,\pm6} & =\phi_{\pm}(\hat{S}_{+}^{6}\pm \hat{S}_{-}^{6})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm5} & =\phi_{\pm}(\hat{S}_{+}^{5}\pm \hat{S}_{-}^{5})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm4} & =\phi_{\pm}(\hat{S}_{+}^{4}\pm \hat{S}_{-}^{4})(11\hat{S}_{z}^{2}-X-38)+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm \hat{S}_{-}^{3})(11\hat{S}_{z}^{3}-(3X+59)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm \hat{S}_{-}^{2})(33\hat{S}_{z}^{4}-(18X+123)\hat{S}_{z}^{2}+X^{2}+10X+102)+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm \hat{S}_{-})(33\hat{S}_{z}^{5}-(30X-15)\hat{S}_{z}^{3}+(5X^{2}-10X+12)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,0} & =231\hat{S}_{z}^{6}-(315X-735)\hat{S}_{z}^{4}+(105X^{2}-525X+294)\hat{S}_{z}^{2}-5X^{3}+40X^{2}-60X+\mathrm{h.c.}
\end{align*}
```

There also exist Stevens operators $\hat{\mathcal{O}}_{k,q}$ for odd $k$.
Because these would violate time-reversal symmetry, Sunny disallows them from
the single-ion anisotropy. Computer-generated tables of Stevens operators with
larger k are available from C. Rudowicz and C. Y. Chung, J. Phys.: Condens.
Matter 16, 5825 (2004).

For each $k$ value, the collection of operators $\{\hat{\mathcal{O}}_{k,q'}\}$
for $q' = -k, \dots, k$ is an irreducible representation of the group of
rotations O(3). That is, a physical rotation will transform
$\hat{\mathcal{O}}_{k,q}$ into a linear combination of
$\hat{\mathcal{O}}_{k,q'}$ where $q'$ varies but $k$ remains fixed. 

In taking the large-$S$ limit, each dipole operator is replaced by its
expectation value $\mathbf{s} = \langle \hat{\mathbf{S}} \rangle$, and only
leading-order terms. That is, $\hat{\mathcal{O}}_{k,q}$ becomes a homogeneous
polynomial $O_{k,q}(\mathbf{s})$ of order $k$ in the spin components. One can
see these polynomials using [`large_S_stevens_operators`](@ref). Due to the
normalization constraint, each dipole is equivalently written as polar angles,
$(\theta, \phi)$. Then the Stevens functions $O_{k,q}(\mathbf{s})$ map almost
exactly to the spherical harmonic functions $Y_l^m(\theta, \phi)$ where $l=k$
and $m=q$ (up to a $k$ and $q$-dependent rescaling factor).
