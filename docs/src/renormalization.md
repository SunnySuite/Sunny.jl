# Interaction Strength Renormalization

A unique feature of Sunny is its support for building classical models where
quantum spin is represented as an $N$-level system, rather than just an expected
dipole. This generalization can be important when modeling quantum spin
Hamiltonians that include, e.g., a single-ion anisotropy, or a biquadratic
coupling between sites. Sunny also supports constraining quantum spin to the
space of pure dipoles; in this case, Sunny will automatically perform an
interaction strength renormalization that maximizes accuracy.

## Local operators

A quantum spin-$S$ state has $N = 2S + 1$ levels. Each local spin operator
$\hat{\mathcal{S}}^{\{x,y,z\}}$ is faithfully represented as an $N×N$ matrix.
These matrices can be accessed using [`spin_matrices`](@ref) for a given label
$S$. For example, the Pauli matrices are associated with $S = 1/2$.

When $S > 1/2$, it is possible to construct multipole moments beyond the
spin-dipole. For example,

```julia
S = spin_matrices(3/2)
@assert S[3] ≈ diagm([3/2, 1/2, -1/2, -3/2])
@assert S[3]^2 ≈ diagm([9/4, 1/4, 1/4, 9/4])
```

If the operator `-S[3]^2` is passed to [`set_onsite_coupling!`](@ref), it would
set an easy-axis anisotropy in the $\hat{z}$ direction.

Any Hermitian operator can be expanded in the basis of Stevens operators
$\hat{\mathcal{O}}_{k,q}$ up to a constant shift. To see this expansion, use
[`print_stevens_expansion`](@ref):
```julia
print_stevens_expansion((S[1]^2 + S[2]^2)) # Prints -(1/3)𝒪₂₀ + 5/2
```

Alternatively, the same operator could have been constructed directly from
[`stevens_matrices`](@ref):

```julia
O = stevens_matrices(3/2)
@assert S[1]^2 + S[2]^2 ≈ -O[2, 0]/3 + (5/2)*I
```

See below for an explicit definition of Stevens operators as polynomials of the
spin operators.

## Renormalization procedure for `:dipole` mode

Sunny will typically operate in one of two modes: `:SUN` or `:dipole`. The
former faithfully represents quantum spin as an SU(_N_) coherent-state which,
for our purposes, is an $N$-component complex vector. In contrast, `:dipole`
mode constrains the coherent-state to the space of pure dipoles. Here, Sunny
will automatically renormalize the magnitude of each Stevens operator to achieve
maximal consistency with `:SUN` mode. This procedure was derived in [D. Dahlbom
et al., [arXiv:2304.03874]](https://arxiv.org/abs/2304.03874).

By way of illustration, consider a quantum operator
$\hat{\mathcal{H}}_{\mathrm{local}}$ giving a single-ion anisotropy for one
site. In Stevens operators,
```math
\hat{\mathcal H}_{\mathrm{local}} = \sum_{k, q} A_{k,q} \hat{\mathcal{O}}_{k,q},
```
for some coefficients $A_{k,q}$.

In `:SUN` mode, Sunny will faithfully represent $\hat{\mathcal
H}_{\mathrm{local}}$ as an $N×N$ matrix. In `:dipole` mode, the expected energy
$\langle \hat{\mathcal H}_{\mathrm{local}} \rangle$ must somehow be approximated
as a function of the expected dipole $\mathbf{s}$.

One possibility is to formally take the $S \to \infty$ limit, whereby each spin
operator $\hat{\mathcal{S}}$ is replaced by its expectation value $\mathbf{s}$.
Correspondingly, each Stevens operator $\hat{\mathcal{O}}_{k,q}$ is replaced by
the Stevens _function_ $\mathcal{O}_{k,q}(\mathbf{s})$, which is a polynomial of
the expected dipole $\mathbf{s}$ rather than of the spin operators
$\hat{\mathcal{S}}$.

In a real magnetic compound, however, the spin magnitude $S$ may not be large,
and to achieve a better approximation one should avoid the large-$S$ limit. The
strategy is to begin with the full dynamics of SU(_N_) coherent states, and then
constrain it to the space of pure dipole states $|\mathbf{s}\rangle$. The latter
are defined such that expectation values,
```math
\langle \mathbf{s}| \hat{\mathcal{S}^\alpha} | \mathbf{s}\rangle = s^\alpha,
```
yield the maximum expected dipole magnitude, $|\mathbf{s}| = S$.

For pure dipole states, it can be demonstrated that
```math
\langle \mathbf{s}| \hat{\mathcal{O}^\alpha} | \mathbf{s}\rangle = c_k \mathcal{O}_{k,q}(\mathbf{s}),
```
where the Stevens functions on the right are scaled by the factors,

```math
\begin{align*}
c_1 &= 1 \\
c_2 &= 1-\frac{1}{2}S^{-1} \\
c_3 &= 1-\frac{3}{2}S^{-1}+\frac{1}{2}S^{-2} \\
c_4 &= 1-3S^{-1}+\frac{11}{4}S^{-2}-\frac{3}{4}S^{-3} \\
c_5 &= 1-5S^{-1}+\frac{35}{4}S^{-2}-\frac{25}{4}S^{-3}+\frac{3}{2}S^{-4} \\
c_6 &= 1-\frac{15}{2}S^{-1}+\frac{85}{4}S^{-2}-\frac{225}{8}S^{-3}+\frac{137}{8}S^{-4}-\frac{15}{4}S^{-5} \\
&\vdots
\end{align*}
```

Collecting results, the SU(_N_) dynamics constrained to the space of dipoles
reduces to the usual Landau-Lifshitz dynamics, but now involving the
_renormalized_ expected energy,
```math
H_{\mathrm{renormalized}}(\mathbf{s}) = \sum_{k, q} c_k A_{k,q} \mathcal{O}_{k,q}(\mathbf{s}).
```

Through these renormalization factors $c_k$, **Sunny avoids the large-$S$
assumption, and gives a more variationally accurate result than traditional
codes**.

Renormalization also applies to the coupling between different sites. In Sunny,
couplings will often be expressed as a polynomial of spin operators using
[`set_pair_coupling!`](@ref), but any such coupling can be decomposed as sum of
tensor products of Stevens operators. Without loss of generality, consider a
single coupling between two Stevens operators $\hat{\mathcal{H}}_\mathrm{biquad}
= \hat{\mathcal{O}}_{k,q} \otimes \hat{\mathcal{O}}_{k',q'}$ along a bond
connecting sites $i$ and $j$. Upon constraining to pure dipole states
$|\mathbf{s}_i\rangle$ and $|\mathbf{s}_j\rangle$, the expected energy takes the
form $c_k c_k' \mathcal{O}_{k,q}(\mathbf{s}_i)
\mathcal{O}_{k',q'}(\mathbf{s}_j)$, which now involves a product of renormalized
Stevens functions. 

## How and when to disable renormalization?

Although we generally recommend the above renormalization procedure, there are
circumstances where it is not desirable. Examples include reproducing a
model-system study, or describing a micromagnetic system for which the
$S\to\infty$ limit is a good approximation. To simulate dipoles without
interaction strength renormalization, construct a [`System`](@ref) using the
mode `:dipole_large_S`. Symbolic operators in the large-$S$ limit can be
constructed by passing `Inf` to either [`spin_matrices`](@ref) or
[`stevens_matrices`](@ref).

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

the relevant Stevens operators are defined as,

```math
\begin{align*}
\hat{\mathcal{O}}_{0,0} & =1\\
\\
\hat{\mathcal{O}}_{1,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm\hat{S}_{-})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{1,0} & =\hat{S}_{z}\\
\\
\hat{\mathcal{O}}_{2,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm\hat{S}_{-}^{2})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{2,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm\hat{S}_{-})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{2,0} & =3\hat{S}_{z}^{2}-X\\
\\
\hat{\mathcal{O}}_{3,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm\hat{S}_{-}^{3})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{3,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm\hat{S}_{-}^{2})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{3,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm\hat{S}_{-})(5\hat{S}_{z}^{2}-X-1/2)+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{3,0} & =5\hat{S}_{z}^{3}-(3X-1)\hat{S}_{z}\\
\\
\hat{\mathcal{O}}_{4,\pm4} & =\phi_{\pm}(\hat{S}_{+}^{4}\pm\hat{S}_{-}^{4})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm\hat{S}_{-}^{3})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm\hat{S}_{-}^{2})(7\hat{S}_{z}^{2}-(X+5))+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm\hat{S}_{-})(7\hat{S}_{z}^{3}-(3X+1)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,0} & =35\hat{S}_{z}^{4}-(30X-25)\hat{S}_{z}^{2}+(3X^{2}-6X)\\
\\
\hat{\mathcal{O}}_{5,\pm5} & =\phi_{\pm}(\hat{S}_{+}^{5}\pm\hat{S}_{-}^{5})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{5,\pm4} & =\phi_{\pm}(\hat{S}_{+}^{4}\pm\hat{S}_{-}^{4})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{5,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm\hat{S}_{-}^{3})(9\hat{S}_{z}^{2}-(X+33/2))+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{5,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm\hat{S}_{-}^{2})(3\hat{S}_{z}^{3}-(X+6)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{5,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm\hat{S}_{-})(21\hat{S}_{z}^{4}-14X\hat{S}_{z}^{2}+(X^{2}-X+3/2))+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{5,0} & =63\hat{S}_{z}^{5}-(70X-105)\hat{S}_{z}^{3}+(15X^{2}-50X+12)\hat{S}_{z}\\
\\
\hat{\mathcal{O}}_{6,\pm6} & =\phi_{\pm}(\hat{S}_{+}^{6}\pm\hat{S}_{-}^{6})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm5} & =\phi_{\pm}(\hat{S}_{+}^{5}\pm\hat{S}_{-}^{5})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm4} & =\phi_{\pm}(\hat{S}_{+}^{4}\pm\hat{S}_{-}^{4})(11\hat{S}_{z}^{2}-X-38)+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm\hat{S}_{-}^{3})(11\hat{S}_{z}^{3}-(3X+59)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm\hat{S}_{-}^{2})(33\hat{S}_{z}^{4}-(18X+123)\hat{S}_{z}^{2}+X^{2}+10X+102)+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm\hat{S}_{-})(33\hat{S}_{z}^{5}-(30X-15)\hat{S}_{z}^{3}+(5X^{2}-10X+12)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,0} & =231\hat{S}_{z}^{6}-(315X-735)\hat{S}_{z}^{4}+(105X^{2}-525X+294)\hat{S}_{z}^{2}-5X^{3}+40X^{2}-60X
\end{align*}
```

Computer-generated tables of Stevens operators with $k > 6$ are available from
[C. Rudowicz and C. Y. Chung, J. Phys.: Condens. Matter 16, 5825
(2004)](https://doi.org/10.1088/0953-8984/16/32/018), but these typically do not
appear in magnetic simulations.

For each $k$ value, the set of operators $\hat{\mathcal{O}}_{k,q'}$ for $q' =
-k, \dots, k$ form an irreducible representation of the group of rotations O(3).
That is, rotation will transform $\hat{\mathcal{O}}_{k,q}$ into a linear
combination of $\hat{\mathcal{O}}_{k,q'}$ where $q'$ varies but $k$ remains
fixed. 

In taking the large-$S$ limit, each dipole operator is replaced by its
expectation value $\mathbf{s} = \langle \hat{\mathbf{S}} \rangle$, and only
leading-order terms are retained. The operator $\hat{\mathcal{O}}_{k,q}$ becomes
a homogeneous polynomial $O_{k,q}(\mathbf{s})$ of order $k$ in the spin
components. One can see these polynomials by constructing
[`stevens_matrices`](@ref) with the argument `S = Inf`. Due to the normalization
constraint, each dipole can be expressed in polar angles, $(\theta, \phi)$. Then
the Stevens functions $O_{k,q}(\mathbf{s})$ correspond to the spherical harmonic
functions $Y_l^m(\theta, \phi)$ where $l=k$ and $m=q$; this correspondence is
valid up to $k$ and $q$-dependent rescaling factors.
