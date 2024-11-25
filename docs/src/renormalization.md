# Interaction Renormalization

A unique feature of Sunny is its support for building classical models where
quantum spin is represented as an $N$-level system, rather than just an expected
dipole. This generalization can be important when modeling quantum spin
Hamiltonians that include, e.g., a single-ion anisotropy, or a biquadratic
coupling between sites. Sunny also supports constraining quantum spin to the
space of pure dipoles; in this case, Sunny will automatically perform an
interaction strength renormalization that enhances accuracy.

## Local operators

A quantum spin-$s$ state has $N = 2s + 1$ levels. Each local spin operator
$\hat{S}^{\{x,y,z\}}$ is faithfully represented as an $N×N$ matrix. These
matrices can be accessed using [`spin_matrices`](@ref) for a given label $s$.
For example, the Pauli matrices are associated with $s = 1/2$.

When $s > 1/2$, it is possible to construct multipole moments beyond the
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

The Stevens operators ``\mathcal{O}_{k, q}`` are [defined below](@ref
"Definition of Stevens operators") as ``k``th order polynomials of the spin
operators.

## Renormalization procedure for `:dipole` mode

Sunny will typically operate in one of two modes: `:SUN` or `:dipole`. The
former faithfully represents quantum spin as an SU(_N_) coherent-state which,
for our purposes, is an $N$-component complex vector. In contrast, `:dipole`
mode constrains the coherent-state to the space of pure dipoles. Here, Sunny
will automatically renormalize the magnitude of each Stevens expectation value
to achieve maximal consistency with `:SUN` mode. This procedure was derived in
[D. Dahlbom et al., [arXiv:2304.03874]](https://arxiv.org/abs/2304.03874).

By way of illustration, consider a quantum operator
$\hat{\mathcal{H}}_{\mathrm{local}}$ giving a single-ion anisotropy for one
site. It can be expanded in [Stevens operators](@ref "Definition of Stevens
operators"),
```math
\hat{\mathcal H}_{\mathrm{local}} = \sum_{k, q} A_{k,q} \hat{\mathcal{O}}_{k,q},
```
for some coefficients $A_{k,q}$.

In `:SUN` mode, Sunny will faithfully represent $\hat{\mathcal
H}_{\mathrm{local}}$ as an $N×N$ matrix. In `:dipole` mode, the expected energy
$\langle \hat{\mathcal H}_{\mathrm{local}} \rangle$ must somehow be approximated
using the expected dipole data.

One approach is to formally take $s \to \infty$, and this yields the traditional
classical limit of a spin system. In this limit spin operators commute and
expectation values of polynomials become polynomials of expectation values. For
example, $\langle \hat{S}^\alpha \hat{S}^\beta\rangle \to \langle \hat{S}^\alpha
\rangle \langle \hat{S}^\beta\rangle$, because any corrections are damped by the
factor $s^{-1} \to 0$. The expectation of a Stevens operator $\langle
\hat{\mathcal{O}}_{k,q} \rangle$ would then become a classical Stevens function
$\mathcal{O}_{k,q}(\langle\hat{\mathbf{S}}\rangle)$, i.e., a polynomial of the
same form, but now applied to the expected dipole. Classical Stevens functions
are constructed as homogeneous polynomials of order $k$, because lower-order
terms would vanish in the limit $s \to \infty$.

For real compounds with finite quantum spin-$s$, one can obtain a better
approximation by avoiding the formal $s \to \infty$ limit. Corrections can be
derived by starting from the full dynamics of SU(_N_) coherent states and then
constraining to the space of pure dipole states $|\boldsymbol{\Omega}\rangle$.
The latter are defined as any states where the expected dipole 3-vector,
```math
\boldsymbol{\Omega} ≡ \langle \boldsymbol{\Omega}| \hat{\mathbf{S}} | \boldsymbol{\Omega}\rangle,
```
has maximal magnitude $|\boldsymbol{\Omega}| = s$ and arbitrary direction.

For a pure dipole state, group theory dictates that expectations of the Stevens
operators can be expressed as a renormalization of the classical Stevens
functions,
```math
\langle \boldsymbol{\Omega}| \hat{\mathcal{O}}_{k,q} | \boldsymbol{\Omega}\rangle = c_k \mathcal{O}_{k,q}(\boldsymbol{\Omega}).
```

At fixed $k$, the two sides must be proportional because they are both spin-$k$
irreducible representations of SO(3). The renormalization factors [can be
calculated explicitly](https://arxiv.org/abs/2304.03874):

```math
\begin{align*}
c_1 &= 1 \\
c_2 &= 1-\frac{1}{2}s^{-1} \\
c_3 &= 1-\frac{3}{2}s^{-1}+\frac{1}{2}s^{-2} \\
c_4 &= 1-3s^{-1}+\frac{11}{4}s^{-2}-\frac{3}{4}s^{-3} \\
c_5 &= 1-5s^{-1}+\frac{35}{4}s^{-2}-\frac{25}{4}s^{-3}+\frac{3}{2}s^{-4} \\
c_6 &= 1-\frac{15}{2}s^{-1}+\frac{85}{4}s^{-2}-\frac{225}{8}s^{-3}+\frac{137}{8}s^{-4}-\frac{15}{4}s^{-5} \\
&\vdots
\end{align*}
```

Constrained to the space of dipoles, the expected local energy becomes
```math
E_{\mathrm{local}}(\boldsymbol{\Omega}) = \langle \boldsymbol{\Omega}| \hat{\mathcal H}_{\mathrm{local}} | \boldsymbol{\Omega}\rangle = \sum_{k, q} c_k A_{k,q} \mathcal{O}_{k,q}(\boldsymbol{\Omega}).
```

It can be shown that SU(_N_) dynamics reduces to the usual Landau-Lifshitz
dynamics of dipoles, but involving $E_{\mathrm{local}}(\boldsymbol{\Omega})$ as
the classical Hamiltonian. The renormalization factors $c_k$ can therefore be
interpreted as a correction to the traditional large-$s$ classical limit.

Renormalization also applies to the coupling between different sites. In Sunny,
couplings will often be expressed as a polynomial of spin operators using
[`set_pair_coupling!`](@ref), but any such coupling can be decomposed as a sum
of tensor products of Stevens operators. Without loss of generality, consider a
single coupling between two Stevens operators
$\hat{\mathcal{H}}_\mathrm{coupling} = \hat{\mathcal{O}}_{k,q} \otimes
\hat{\mathcal{O}}_{k',q'}$ along a bond connecting sites $i$ and $j$. Upon
constraining to pure dipole states $|\boldsymbol{\Omega}_i\rangle$ and
$|\boldsymbol{\Omega}_j\rangle$, the expected energy takes the form
$E_\mathrm{coupling} = c_k c_k' \mathcal{O}_{k,q}(\boldsymbol{\Omega}_i)
\mathcal{O}_{k',q'}(\boldsymbol{\Omega}_j)$, which now involves a product of
renormalized Stevens functions. 

## Use `:dipole_uncorrected` mode to disable renormalization

The above renormalization procedure is valid under the assumption that the
starting model is the true microscopic quantum Hamiltonian. Sometimes, however,
the starting model is instead an effective classical Hamiltonian that has been
fitted to experimental data. In this case, the model parameters will already
incorporate any appropriate renormalizations, and no further renormalization
should be applied. Similarly, renormalization is unneeded for effective models
of micromagnets that are far from the quantum regime.

To specify an effective spin-dipole Hamiltonian with renormalization disabled,
construct a [`System`](@ref) using the mode `:dipole_uncorrected` instead of
`:dipole`. Formally, `:dipole_uncorrected` takes the $s → ∞$ limit, such that
all local operators become infinite dimensional and commute. A symbolic
representation of these operators can be obtained by passing `Inf` to either
[`spin_matrices`](@ref) or [`stevens_matrices`](@ref). Polynomials of such spin
operators can be used, e.g., in [`set_onsite_coupling!`](@ref).

## Definition of Stevens operators

The Stevens operators $\hat{\mathcal{O}}_{k,q}$ are defined as polynomials of
angular momentum operators $\hat{S}_{\{x,y,z\}}$ in some spin-$s$
representation.

Using

```math
\begin{align*}
X &= \mathbf{\hat{S}} \cdot \mathbf{\hat{S}} = s (s+1) \\
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

The case $k=1$ gives the dipole operators,
```math
(\hat{\mathcal{O}}_{1,1}, \hat{\mathcal{O}}_{1,0}, \hat{\mathcal{O}}_{1,-1}) = (\hat{S}_{x}, \hat{S}_{z}, \hat{S}_{y}).
```

The case $k=2$ gives the quadrupole operators,
```math
(\hat{\mathcal{O}}_{2,2}, \dots, \hat{\mathcal{O}}_{2,-2}) = \left(\hat{S}_x^2 - \hat{S}_y^2, \frac{\hat{S}_x \hat{S}_z + \hat{S}_z \hat{S}_x}{2}, 2\hat{S}_z^2-\hat{S}_x^2-\hat{S}_y^2, \frac{\hat{S}_y \hat{S}_z + \hat{S}_z \hat{S}_y}{2}, \hat{S}_x \hat{S}_y + \hat{S}_y \hat{S}_x\right).
```

For each $k$ value, the set of operators $\hat{\mathcal{O}}_{k,q'}$ for $q' =
-k, \dots, k$ form an irreducible representation of the group of rotations O(3).
That is, rotation will transform $\hat{\mathcal{O}}_{k,q}$ into a linear
combination of $\hat{\mathcal{O}}_{k,q'}$ where $q'$ varies but $k$ remains
fixed. 

In taking the large-$s$ limit, each dipole operator is replaced by its
expectation value $\boldsymbol{\Omega} = \langle \hat{\mathbf{S}} \rangle$, and
only leading-order terms are retained. The operator $\hat{\mathcal{O}}_{k,q}$
becomes a homogeneous polynomial $O_{k,q}(\boldsymbol{\Omega})$ of order $k$ in
the spin components $\Omega^\alpha$. One can see these polynomials by
constructing [`stevens_matrices`](@ref) with the argument `s = Inf`. Due to the
normalization constraint, each dipole can be expressed in polar angles,
$(\theta, \phi)$. Then the Stevens functions $O_{k,q}(\boldsymbol{\Omega})$
correspond to the spherical harmonic functions $Y_l^m(\theta, \phi)$ where $l=k$
and $m=q$; this correspondence is valid up to $k$ and $q$-dependent rescaling
factors.
