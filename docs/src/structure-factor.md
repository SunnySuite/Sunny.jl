# Structure Factor Conventions

## Dynamical correlations

Dynamical correlations are a fundamental observable in condensed matter systems,
and facilitate comparison between theory and experimental data.

Frequently, spin-spin correlations are of interest. More generally, one may
consider correlations between arbitrary operator fields $\hat{A}$ and $\hat{B}$.
In the Heisenberg picture, operators evolve in time. Consider the dynamical
correlation as an equilibrium expectation value,

```math
\begin{equation}
C(𝐫, t) = \int_V ⟨\hat{B}^†(𝐫₀, 0) \hat{A}(𝐫₀ + 𝐫, t)⟩ d𝐫₀,
\end{equation}
```

where the integral runs over some macroscopic volume $V → ∞$. By construction,
$C(𝐫, t)$ is an extensive quantity. Ignoring surface effects, the correlation
becomes an ordinary product in momentum-space,

```math
\begin{equation}
C(𝐪, t) = ⟨\hat{B}_𝐪^†(0) \hat{A}_𝐪(t)⟩.
\end{equation}
```

Our convention for the Fourier transform from position $𝐫$ to momentum $𝐪$ is,

```math
\begin{equation}
\hat{A}_𝐪 ≡ \int_V e^{+ i 𝐪⋅𝐫} \hat{A}(𝐫) d𝐫.
\end{equation}
```

For a Hermitian operator $\hat{A}^†(𝐫) = \hat{A}(𝐫)$, it follows that
$\hat{A}_𝐪^† ≡ (\hat{A}_𝐪)^† = \hat{A}_{-𝐪}$ in momentum space. 


## Lehmann representation in frequency space

Dynamical correlations are most conveniently calculated in frequency space. We
use the convention,

```math
\begin{equation}
C(𝐪, ω) ≡ \frac{1}{2π} \int_{-∞}^{∞} e^{-iωt} C(𝐪, t) dt.
\end{equation}
```

On the right-hand side, substitute the definition of Heisenberg time evolution,

```math
\begin{equation}
\hat{A}(t) ≡ e^{i t \hat{H}} \hat{A} e^{-i t \hat{H}},
\end{equation}
```

and the definition of the thermal average,

```math
\begin{equation}
⟨\hat{O}⟩ ≡ \frac{1}{\mathcal{Z}}\mathrm{tr}\, e^{-β \hat{H}} \hat{O}, \quad \mathcal{Z} ≡ \mathrm{tr}\, e^{-β \hat{H}}.
\end{equation}
```

Let $ϵ_μ$ and $|μ⟩$ denote the exact eigenvalues and eigenstates of the
full Hamiltonian $\hat{H}$. The eigenstates comprise a complete, orthonormal
basis. The operator trace can be evaluated as a sum over eigenbasis states
$|ν⟩$. Insert a resolution of the identity, $|μ⟩⟨μ| = 1$, with implicit
summation on the repeated $μ$ index. Then, collecting results and applying,

```math
\begin{equation}
\int_{-∞}^{∞} e^{-iωt} dt = 2πδ(ω),
\end{equation}
```

the result is the Lehmann representation of dynamical correlations,

```math
\begin{equation}
C(𝐪,ω) = \frac{1}{\mathcal{Z}} e^{-β ϵ_ν} δ(ϵ_μ - ϵ_ν - ω) ⟨ν|\hat{B}^†_𝐪|μ⟩⟨μ|\hat{A}_𝐪 |ν⟩,
\end{equation}
```

with implicit summation over eigenbasis indices $μ$ and $ν$. This representation
is the usual starting point for quantum calculations such as linear spin wave
theory and its various generalizations.

Using the Lehmann representation, it can be shown that positive and negative
frequencies are linked through a detailed balance condition,

```math
\begin{equation}
C_{⟨BA^†⟩}(𝐪,-ω) = e^{-β ω}  C_{⟨B^†A⟩}(𝐪, ω)^*.
\end{equation}
```

This subscript notation indicates that the left-hand side is a correlation of
Hermitian-conjugated operators. Typically $\hat{A}$ and $\hat{B}$ will be
Hermitian in real-space, and then detailed balance becomes $C(-𝐪,-ω) = e^{-β ω}
C(𝐪, ω)^*$.


## Discrete sums on the lattice

A chemical (crystallographic) unit cell is associated with three lattice vectors
$𝐚_{\{1,2,3\}}$. Site positions are

```math
\begin{equation}
𝐫_{𝐦,j} ≡ m_1 𝐚_1 + m_2 𝐚_2 + m_3 𝐚_3 + δ𝐫_j,
\end{equation}
```

for integers $𝐦 = \{m_1, m_2, m_3\}$. If the crystal is decorated, then $δ𝐫_j$
denotes the relative displacement of the Bravais sublattice $j$. 

Let $\hat{A}(𝐫)$ be decomposed into discrete contributions $\hat{A}_{𝐦,j}
δ(𝐫-𝐫_{𝐦,j})$ at each lattice point $𝐫_{𝐦,j}$. The Fourier transform
becomes a discrete sum,

```math
\begin{equation}
\hat{A}_𝐪 = \sum_j \sum_𝐦 e^{i 𝐪⋅𝐫_{𝐦,j}} \hat{A}_{𝐦,j} ≡ \sum_j \hat{A}_{𝐪,j}.
\end{equation}
```

The second equality above introduces $\hat{A}_{𝐪,j}$ as the Fourier transform
of $\hat{A}_{𝐦,j}$ for _single_ sublattice $j$. It can also be written,

```math
\begin{equation}
\hat{A}_{𝐪,j} = e^{i 𝐪⋅δ𝐫_j} \sum_𝐦 e^{i 2π \tilde{𝐪}⋅𝐦} \hat{A}_{𝐦,j},
\end{equation}
```

where $\tilde{𝐪}$ expresses momentum in dimensionless reciprocal lattice units
(RLU),

```math
\begin{equation}
𝐪 = \tilde{q}_1 𝐛_1 + \tilde{q}_2 𝐛_2 + \tilde{q}_3 𝐛_3,
\end{equation}
```

and $𝐛_{\{1,2,3\}}$ are the reciprocal lattice vectors. Equivalently,
$\tilde{q}_μ ≡ 𝐪 ⋅ 𝐚_μ / 2π$.

It will be convenient to introduce a dynamical correlation for the operators on
sublattices $i$ and $j$ only,

```math
\begin{equation}
C_{ij}(𝐪,t) ≡ ⟨\hat{B}^†_{𝐪,i}(0) \hat{A}_{𝐪,j}(t)⟩.
\end{equation}
```

By the linearity of expectation values,

```math
\begin{equation}
C(𝐪, t) = \sum_{ij} C_{ij}(𝐪,t).
\end{equation}
```


## Quantum sum rule

Integrating over all frequencies $ω$ yields the instant correlation at real-time
$t = 0$,

```math
\begin{equation}
\int_{-∞}^∞ C(𝐪,ω) dω = C(𝐪, t=0) = \sum_{ij} C_{ij}(𝐪, t=0).
\end{equation}
```

Here, we will investigate spin-spin correlations. For this, select
$\hat{B}_{𝐪,i} = \hat{𝐒}_{𝐪,i}$ and $A_{𝐪,j} = \hat{𝐒}_{𝐪,j}$,
such that the dynamical correlations become tensor valued,

```math
\begin{equation}
C_{ij}^{αβ}(𝐪, t=0) = ⟨\hat{S}_{𝐪,i}^{α†} \hat{S}_{𝐪,j}^{β}⟩.
\end{equation}
```

In the quantum spin-$s$ representation, the spin dipole on one site satisfies

```math
\begin{equation}
|\hat{𝐒}|^2 = \hat{S}^α \hat{S}^α = s(s+1),
\end{equation}
```

with implicit summation on the repeated $α$ index.

Suppose that each site of sublattice $j$ carries quantum spin of magnitude
$s_j$. Then there is a quantum sum rule of the form,

```math
\begin{equation}
\int_{\tilde{V}_\mathrm{BZ}} \frac{C_{jj}^{αα}(𝐪, t=0)}{N_\mathrm{cells}} d\tilde{𝐪} = s_j (s_j + 1),
\end{equation}
```

with summation on $α$, but not $j$, implied. The integral runs over the cubic
volume in reciprocal lattice units $\tilde{𝐪}$,

```math
\begin{equation}
\tilde{V}_\mathrm{BZ} ≡ [0,1]^3.
\end{equation}
```

This volume represents one unit cell on the reciprocal lattice, and has the
shape of a parallelepiped in physical momentum units $𝐪$. This volume is
equivalent to the first Brillouin zone because of the reciprocal-space
periodicity inherent to the Bravais sublattice. Note that the integral over
$\tilde{𝐪} ∈ \tilde{V}_\mathrm{BZ}$ could be converted to an integral over
physical momentum $𝐪$ by applying a Jacobian transformation factor, $d
\tilde{𝐪} = d𝐪 V_\mathrm{cell} / (2π)^3$, where $V_\mathrm{cell} = |𝐚_1 ⋅
(𝐚_2 × 𝐚_3)|$ is the volume of the chemical unit cell. The scaling factor

```math
\begin{equation}
N_\mathrm{cells} ≡ V / V_\mathrm{cell}
\end{equation}
```

denotes the number of chemical unit cells in the macroscopic volume $V$.

The derivation of the sum rule proceeds as follows. Substitute twice the
definition,

```math
\begin{equation}
\hat{S}^α_{𝐪,j} ≡ e^{i 𝐪⋅δ𝐫_j} \sum_𝐦 e^{i 2π \tilde{𝐪}⋅𝐦} \hat{S}^α_{𝐦,j}.
\end{equation}
```

Accounting for complex conjugation, the two phase factors $e^{i 𝐪⋅δ𝐫_j}$
cancel. The remaining $𝐪$-dependence can be integrated to yield a
Kronecker-$δ$,

```math
\begin{equation}
\int_{\tilde{V}_\mathrm{BZ}} e^{2πi \tilde{𝐪} ⋅ (𝐦 - 𝐦') } d\tilde{𝐪} = δ_{𝐦, 𝐦'}.
\end{equation}
```

Note that $⟨\hat{S}_{𝐦,j}^α \hat{S}_{𝐦,j}^α⟩ = s_j(s_j+1)$ is constant,
independent of the cell $𝐦$. This leaves a double sum over integers $𝐦$, which
evaluates to $\sum_{𝐦, 𝐦'} δ_{𝐦, 𝐦'} = N_\mathrm{cells}$. Combined, these
results verify the above-stated quantum sum rule for the sublattice $j$.

One can also derive a quantum sum rule on the full dynamical correlation $C^{α,
β}(𝐪, ω)$. Contributions from distinct sublattices $i ≠ j$ introduce a phase
factor $e^{- i 𝐪⋅(δ𝐫_i - δ𝐫_j)}$ that cancels when the momentum $𝐪$ is
averaged over a large number $N_\mathrm{BZ} → ∞$ of Brillouin zones. The final
result is a sum over contributions $C_{jj}(𝐪, t=0)$ for each sublattice $j$,

```math
\begin{equation}
\frac{1}{N_\mathrm{BZ}} \int_{N_\mathrm{BZ} × \tilde{V}_\mathrm{BZ}} \int_{-∞}^∞ \frac{C^{αα}(𝐪, ω)}{ N_\mathrm{cells}} dω d\tilde{𝐪} = \sum_j s_j (s_j + 1).
\end{equation}
```

## Neutron scattering cross section

The magnetic moment of a neutron is $\hat{\boldsymbol{μ}}_\mathrm{neutron} = - 2
γ μ_N \hat{𝐒}_\mathrm{neutron}$, where $γ = 1.913…$, $μ_N$ is the nuclear
magneton, and $\hat{𝐒}_\mathrm{neutron}$ is spin-1/2 angular momentum. Neutrons
interact with the magnetic moments of a material. These have the form
$\hat{\boldsymbol{μ}} = -μ_B g \hat{𝐒}$, where $μ_B$ is the Bohr magneton and
$\hat{𝐒}$ is the effective angular momentum. For a single electron, $g =
2.0023…$ is known to high precision. Within a crystal, however, the appropriate
$g_j$ for each sublattice $j$ may be any $3×3$ matrix consistent with point
group symmetries.

Each idealized magnetic moment $\hat{\boldsymbol{μ}}_{𝐦,j}$ is, in reality,
smoothly distributed around the site position $𝐫_{𝐦, j}$. This can be modeled
through convolution with a density function $f_j(𝐫)$. Fourier transform the
full magnetic density field $𝐌(𝐫)$ to obtain

```math
\begin{equation}
\hat{𝐌}_𝐪 ≡ \sum_j \hat{𝐌}_{𝐪,j},
\end{equation}
```

where,

```math
\begin{equation}
\hat{𝐌}_{𝐪,j} ≡ - μ_B e^{i 𝐪⋅δ𝐫_j} g_j \sum_𝐦 e^{i 2π \tilde{𝐪}⋅𝐦} \hat{𝐒}_{𝐦,j} f_j(𝐪).
\end{equation}
```

In Fourier space, $f_j(𝐪)$ is called the _magnetic form factor_. Frequently, it
will be approximated as an isotropic function of $q = |𝐪|$. Tabulated formula,
for various magnetic ions and charge states, are available in Sunny via the
[`FormFactor`](@ref) function. The idealized case $f_j(𝐪) = 1$ would describe
completely localized magnetic moments.

Neutron scattering intensities are given by the total differential
cross-section, $d^2 σ(𝐪, ω)/dωdΩ$, where $𝐪 = 𝐪_i - 𝐪_f$ is the momentum
transfer to the sample, $ω$ is the energy transfer to the sample, and $Ω$ is the
solid angle. Experimental intensity data will typically be provided in units of
$q_f / q_i$. Within the dipole approximation, the result for an unpolarized
neutron beam is,

```math
\begin{equation}
\frac{d^2 σ(𝐪, ω)}{dω dΩ} = \frac{q_f}{q_i} \left(\frac{γ r_0}{2}\right)^2 \sum_{α,β} \left(δ_{α,β} - \frac{q^α q^β}{q^2}\right) \frac{\mathcal{S}^{αβ}(𝐪, ω)}{μ_B^2}.
\end{equation}
```

The prefactor $q_f/q_i$ will be provided experimentally. Dimensions of area
arise from the characteristic scattering length, $γ r_0 / 2 ≈ 2.69×10^{-5}
\mathrm{Å}$, where $r_0$ is the classical electron radius.

The structure factor is of central importance to neutron scattering,

```math
\begin{equation}
\mathcal{S}^{αβ}(𝐪, ω) ≡ \frac{1}{2π} \int_{-∞}^{∞} e^{-iωt} ⟨\hat{M}_𝐪^{α†}(0) \hat{M}_𝐪^β(t)⟩ dt,
\end{equation}
```

and describes dynamical correlations of magnetic moments. It will differ
nontrivially from the spin-spin correlations if the $g_j$-tensor varies with
sublattice $j$.

## Conventions for the Sunny-calculated structure factor

Calculating the structure factor involves several steps, with various possible
settings. Sunny provides tools to facilitate this calculation and to extract
information from the results. For details, please see our [tutorials](@ref "1.
Spin wave simulations of CoRh₂O₄") as well as the complete [Library
Reference](@ref).

Through [`ssf_custom`](@ref) and related functions, Sunny will calculate the
spin structure factor as a 3×3 matrix in dimensionless units,

```math
\begin{equation}
𝒮^{αβ}(𝐪, ω) ≡ \frac{1}{N_\mathrm{cells} μ_B^2} \mathcal{S}^{αβ}(𝐪, ω),
\end{equation}
```

This is an intensive quantity because $N_\mathrm{cells}$, the number of chemical
cells in the macroscopic sample, is extensive. Note that the Sunny-calculated
intensity will depend on the chemical cell convention: intensity scales linearly
with chemical cell size.

Use [`ssf_perp`](@ref) to contract with $δ_{α,β} - q^α q^β/q^2$, i.e., to
project in the direction perpendicular to momentum transfer $𝐪$.

Set `apply_g = false` to calculate the correlation $C_{⟨𝐒𝐒⟩}(𝐪, ω) /
N_\mathrm{cells}$ between pure spin operators, rather than between magnetic
moments. In this special case that $g$ is a uniform scalar, this is equivalent
to $𝒮(𝐪, ω) / g^2$

## Calculations with spin wave theory

Calculating the dynamical structure factor with linear [`SpinWaveTheory`](@ref)
is relatively direct. In the traditional approach, quantum spin operators are
expressed with Holstein-Primakoff bosons, and dynamical correlations are
calculated to leading order in inverse powers of the quantum spin-$s$. For
systems constructed with `mode = :SUN`, Sunny automatically switches to a
multi-flavor boson variant of spin wave theory, which captures more single-ion
physics. Use [`SpinWaveTheorySpiral`](@ref) to study generalized spiral phases,
which allow for an incommensurate propagation wavevector. The module
[`SpinWaveTheoryLarge`](@ref) uses iterated matrix-vector multiplication to
achieve linear-scaling cost in the magnetic cell size. It can be useful for
studying systems with large magnetic cells, including systems with
long-wavelength structures, or systems with quenched chemical disorder.


## Calculations with classical spin dynamics

Finite temperature structure factor intensities can be estimated from the
dynamical correlations of classical spin dynamics (e.g. Landau-Lifshitz, or its
SU($N$) generalization). This is fundamentally a Monte Carlo approach, as the
trajectories must be initialized to a spin configuration that is sampled from
the finite-temperature thermal equilibrium. Samples are accumulated into a
[`SampledCorrelations`](@ref), from which intensity information may be
extracted.

Creating a `SampledCorrelations` requires specifying three keyword arguments.
These will determine the dynamics used to calculate samples and, consequently,
the $ω$ information that will be available. 

1. `energies`: A uniform range of resolved energies.
2. `dt`: The step size for dynamical time-integration. Larger may reduce
   simulation time, but the choice will be limited by the stability and accuracy
   requirements of the [`ImplicitMidpoint`](@ref) integration method. The
   function [`suggest_timestep`](@ref) can recommend a good value. The timestep
   may be adjusted downward so that the specified `energies` are sampled
   exactly.
3. `measure`: Specification of the pair correlations. This will frequently be
   reduced from the spin structure factor using one of [`ssf_trace`](@ref),
   [`ssf_perp`](@ref), or [`ssf_custom_bm`](@ref).

A sample may be added by calling [`add_sample!`](@ref). The input `sys` must be
a spin configuration in good thermal equilibrium, e.g., using the continuous
[`Langevin`](@ref) dynamics or using single spin flip trials with
[`LocalSampler`](@ref). The statistical quality of the $𝒮(𝐪,ω)$ can be
improved by repeatedly generating decorrelated spin configurations in `sys` and
calling `add_sample!` on each configuration.

The outline of typical use case might look like this:
```
# Make a `SampledCorrelations`
sc = SampledCorrelations(sys; dt=0.05, energies=range(0.0, 10.0, 100))

# Add samples
for _ in 1:nsamples
   decorrelate_system(sys) # Perform some type of Monte Carlo simulation
   add_sample!(sc, sys)    # Use spins to calculate trajectory and accumulate new sample of 𝒮(𝐪,ω)
end
```

### Extracting intensities sampled correlation data 

Like in spin wave theory, the basic function for extracting intensities
information from a `SampledCorrelations` is [`intensities`](@ref). It takes a
`SampledCorrelations`, a collection of ``𝐪``-vectors, a collection of
`energies`, and possible other options.

Since classical dynamics simulation take place on a finite lattice, the
fundamental intensities measurements are only available at a discrete grid of
wave vectors. In reciprocal lattice units, available grid points are $𝐪 =
[\frac{n_1}{L_1}, \frac{n_2}{L_2}, \frac{n_3}{L_3}]$, where $n_i$ runs from
$(\frac{-L_i}{2}+1)$ to $\frac{L_i}{2}$ and $L_i$ is the linear dimension of the
lattice used in the calculation. (An internal function
`Sunny.available_wave_vectors` provides access to this grid.) By default
[`intensities`](@ref) will adjust each wavevector $𝐪$ to the nearest available
grid point. 

Similarly, the resolution in `energies` is controlled the dynamical trajectory
length in real-time. Because the dynamical trajectory is not periodic in time,
some blurring between neighboring energy bins is unavoidable. Sunny's algorithm
for estimating the structure factor from real-time dynamics is [specified
here](https://github.com/SunnySuite/Sunny.jl/pull/246#issuecomment-2119294846).

The temperature parameter `kT` is required for `SampledCorrelations`
calculations, and will be used to perform classical-to-quantum rescaling of
intensities. If `kT = nothing`, then intensities will be provided according to
the classical Boltzmann distribution.

### The instantaneous structure factor

Use [`intensities_static`](@ref) to calculate $\mathcal{S}(𝐪)$, i.e.,
 correlations that are "instantaneous" in real-time. Mathematically,
$\mathcal{S}(𝐪)$ denotes an integral of the dynamical structure factor
$\mathcal{S}(𝐪, ω)$ over all energies $ω$. In [`SpinWaveTheory`](@ref), the
energy integral becomes a discrete sum over bands. In
[`SampledCorrelations`](@ref), a classical-to-quantum correction factor will be
applied within [`intensities`](@ref) prior to energy integration.

Sunny also supports a mechanism to calculate static correlations without any
spin dynamics. To collect such statistics, construct a
`SampledCorrelationsStatic` object. In this case, `intensities_static` will
return static correlations sampled from the classical Boltzmann distribution.
This dynamics-free approach is faster, but may miss important features that
derive from the quantum mechanical excitation spectrum.
