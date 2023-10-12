# Interaction Strength Renormalization

A unique feature of Sunny is its support for building classical models where
each quantum spin is represented as a full $N$-level system, rather than just an
expected dipole. This formalism enables more accurate modeling of quantum spin
Hamiltonians that include, e.g., a single-ion anisotropy, or a biquadratic
coupling between sites.

## Local operators

A quantum spin-$S$ state has $N = 2S + 1$ levels. Each local spin operator
$\hat{S}^{\{x,y,z\}}$ is faithfully represented as an $N×N$ matrix. Access these
matrices using [`spin_matrices`](@ref) for a given label $S$. For example,
`spin_matrices(1/2)` returns the Pauli matrices divided by 2.

When $S > 1/2$, it is possible to construct multipole moments beyond the
spin-dipole. For example,

XXX TODO XXX
```julia
S = spin_matrices(2)
```

The Stevens operators
$\hat{\mathcal{O}}_{k,q}$ are polynomials of the spin operators, and are
accessed using [`stevens_matrices`](@ref). With these building blocks, a
single-ion anisotropy is defined using [`set_onsite_coupling!`](@ref). For
example:

```julia
# An easy axis anisotropy in the z-direction
S = spin_operators(sys, i)
set_onsite_coupling!(sys, -D*S[3]^3, i)

# The unique quartic single-ion anisotropy for a site with cubic point group
# symmetry
O = stevens_operators(sys, i)
set_onsite_coupling!(sys, O[4,0] + 5*O[4,4], i)

# An equivalent expression of this quartic anisotropy, up to a constant shift
set_onsite_coupling!(sys, 20*(S[1]^4 + S[2]^4 + S[3]^4), i)
```


## Renormalization procedure for `:dipole` mode

There are two allowed modes for a [`System`](@ref). The mode `:SUN` models each
spin as an SU(_N_) coherent state (i.e., as a set of $N$ complex amplitudes),
and is the most variationally accurate. The mode `:dipole` constrains the
SU(_N_) coherent-state dynamics to the space of pure dipoles. In either mode,
Sunny encourages specifying single-ion anisotropies as $N×N$ matrices. In
`:dipole` mode, Sunny will automatically renormalize the anisotropy operator to
achieve maximal consistency with `:SUN` mode. This procedure was derived in [D.
Dahlbom et al., [arXiv:2304.03874]](https://arxiv.org/abs/2304.03874). Here, we
summarize the final results.

The starting point is a quantum operator $\hat{\mathcal{H}}_{\mathrm{local}}$
giving the single-ion anisotropy for one site. It can be expanded in Stevens
operators,
```math
\hat{\mathcal H}_{\mathrm{local}} = \sum_{k, q} A_{k,q} \hat{\mathcal{O}}_{k,q}.
```

See the documentation of [`print_stevens_expansion`](@ref) for some explicit
examples of this expansion.

The traditional classical limit of a quantum spin Hamiltonian, which yields the
Landau-Lifshitz dynamics, can be derived by taking the formal $S \to\infty$
limit, such that each spin operator $\hat{\mathbf{S}}$ is replaced by its dipole
expectation value $\mathbf{s}$. Correspondingly, the Stevens operators
$\hat{\mathcal{O}}_{k,q}$ become polynomials $\mathcal{O}_{k,q}(\mathbf{s})$ in
the classical dipole. With this traditional approach, one would arrive at the
_bare_ expected energy,
```math
H_{\mathrm{bare}}(\mathbf{s}) = \sum_{k, q} A_{k,q} \mathcal{O}_{k,q}(\mathbf{s}).
```

In a real magnetic compound, however, $S$ may not be very large, and one can
achieve a better approximation by avoiding the $S \to\infty$ limit. The strategy
is to begin with the full dynamics of SU(_N_) coherent states, and then
constrain it to the space of dipoles $\mathbf{s}$. Doing so will again yield the
Landau-Lifshitz dynamics, but now involving the _renormalized_ expected energy,
```math
H_{\mathrm{renormalized}}(\mathbf{s}) = \sum_{k, q} c_k A_{k,q} \mathcal{O}_{k,q}(\mathbf{s}).
```
The $k$-dependent renormalization factors are
```math
\begin{align*}
c_2 &= 1-\frac{1}{2}S^{-1} \\
c_4 &= 1-3S^{-1}+\frac{11}{4}S^{-2}-\frac{3}{4}S^{-3} \\
c_6 &= 1-\frac{15}{2}S^{-1}+\frac{85}{4}S^{-2}-\frac{225}{8}S^{-3}+\frac{137}{8}S^{-4}-\frac{15}{4}S^{-5}.
\end{align*}
```

Sunny will use $H_{\mathrm{renormalized}}(\mathbf{s})$ in its classical dynamics
of dipoles. Because of this renormalization, **Sunny is more variationally
accurate than traditional codes like SpinW**.

## How and when to disable renormalization?

Although we generally recommend the above renormalization procedure, there are
circumstances where it is not desirable. Examples include reproducing a
model-system study, or describing a micromagnetic system for which the
$S\to\infty$ limit is quantitatively realized. To get symbolic operators in the
large-$S$ limit, use [`spin_matrices`](@ref) or [`stevens_matrices`](@ref) with
the argument `Inf`. Sunny will not perform any renormalization on anisotropy
operators constructed through these primitives.

Note that Sunny will _also_ renormalize scalar biquadratic exchange interactions
by default. Disable this renormalization by setting `large_S = true` in the call
to [`set_exchange!`](@ref).

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
\hat{\mathcal{O}}_{0,0} & = 1 \\
\\
\hat{\mathcal{O}}_{2,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm \hat{S}_{-}^{2})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{2,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm \hat{S}_{-})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{2,0} & =3\hat{S}_{z}^{2}-X\\
\\
\hat{\mathcal{O}}_{4,\pm4} & =\phi_{\pm}(\hat{S}_{+}^{4}\pm \hat{S}_{-}^{4})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm \hat{S}_{-}^{3})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm \hat{S}_{-}^{2})(7\hat{S}_{z}^{2}-(X+5))+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm \hat{S}_{-})(7\hat{S}_{z}^{3}-(3X+1)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{4,0} & =35\hat{S}_{z}^{4}-(30X-25)\hat{S}_{z}^{2}+(3X^{2}-6X)\\
\\
\hat{\mathcal{O}}_{6,\pm6} & =\phi_{\pm}(\hat{S}_{+}^{6}\pm \hat{S}_{-}^{6})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm5} & =\phi_{\pm}(\hat{S}_{+}^{5}\pm \hat{S}_{-}^{5})\hat{S}_{z}+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm4} & =\phi_{\pm}(\hat{S}_{+}^{4}\pm \hat{S}_{-}^{4})(11\hat{S}_{z}^{2}-X-38)+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm3} & =\phi_{\pm}(\hat{S}_{+}^{3}\pm \hat{S}_{-}^{3})(11\hat{S}_{z}^{3}-(3X+59)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm2} & =\phi_{\pm}(\hat{S}_{+}^{2}\pm \hat{S}_{-}^{2})(33\hat{S}_{z}^{4}-(18X+123)\hat{S}_{z}^{2}+X^{2}+10X+102)+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,\pm1} & =\phi_{\pm}(\hat{S}_{+}\pm \hat{S}_{-})(33\hat{S}_{z}^{5}-(30X-15)\hat{S}_{z}^{3}+(5X^{2}-10X+12)\hat{S}_{z})+\mathrm{h.c.}\\
\hat{\mathcal{O}}_{6,0} & =231\hat{S}_{z}^{6}-(315X-735)\hat{S}_{z}^{4}+(105X^{2}-525X+294)\hat{S}_{z}^{2}-5X^{3}+40X^{2}-60X
\end{align*}
```

Stevens operators $\hat{\mathcal{O}}_{k,q}$ for odd $k$ are disallowed from the
single-ion anisotropy under the assumption of time-reversal symmetry.
Computer-generated tables of Stevens operators with larger k are available from
C. Rudowicz and C. Y. Chung, J. Phys.: Condens. Matter 16, 5825 (2004).

For each $k$ value, the collection of operators $\{\hat{\mathcal{O}}_{k,q'}\}$
for $q' = -k, \dots, k$ is an irreducible representation of the group of
rotations O(3). In particular, a physical rotation will transform
$\hat{\mathcal{O}}_{k,q}$ into a linear combination of
$\hat{\mathcal{O}}_{k,q'}$ where $q'$ varies but $k$ remains fixed. 

In taking the large-$S$ limit, each dipole operator is replaced by its
expectation value $\mathbf{s} = \langle \hat{\mathbf{S}} \rangle$, and only
leading-order terms are retained. The operator $\hat{\mathcal{O}}_{k,q}$ becomes
a homogeneous polynomial $O_{k,q}(\mathbf{s})$ of order $k$ in the spin
components. One can see these polynomials by constructing
[`stevens_matrices`](@ref) with the argument `S = Inf`. Due to the normalization
constraint, each dipole can be expressed in polar angles, $(\theta, \phi)$. Then
the Stevens functions $O_{k,q}(\mathbf{s})$ correspond to the spherical harmonic
functions $Y_l^m(\theta, \phi)$ where $l=k$ and $m=q$, and modulo $k$ and
$q$-dependent rescaling factors.
