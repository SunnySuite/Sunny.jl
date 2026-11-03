# Structure Factor Conventions

## What Sunny calculates

The magnetic structure factor is the Fourier transform of the dynamical
correlations between magnetic moments,

```math
\begin{equation}
\mathcal{S}^{αβ}(𝐪, ω) ≡ \frac{1}{2π} \int_{-∞}^{∞} e^{-iωt} ⟨\hat{M}_𝐪^{α†}(0) \hat{M}_𝐪^β(t)⟩ dt.
\end{equation}
```

The Fourier-transformed magnetic moment density is

```math
\begin{equation}
\hat{𝐌}_𝐪 ≡ - μ_B \sum_j f_j(𝐪) g_j \hat{𝐒}_{𝐪,j}, \quad \hat{𝐒}_{𝐪,j} ≡ \sum_𝐦 e^{i 𝐪⋅𝐫_{𝐦,j}} \hat{𝐒}_{𝐦,j}.
\end{equation}
```

Here $\hat{𝐒}_{𝐦,j}$ is the spin dipole operator at position $𝐫_{𝐦,j}$ for
cell $𝐦$ and sublattice $j$ of the chemical cell, the $g_j$-tensor maps
effective spin to magnetic moment, $μ_B$ is the Bohr magneton, and $f_j(𝐪)$ is
the magnetic form factor.

Sunny reports the normalized structure factor in units of inverse energy,

```math
\begin{equation}
𝒮^{αβ}(𝐪, ω) ≡ \frac{1}{N_\mathrm{cells} μ_B^2} \mathcal{S}^{αβ}(𝐪, ω).
\end{equation}
```

It is normalized by $N_\mathrm{cells}$, the number of chemical cells in the
macroscopic sample. Observe that the scale of $𝒮^{αβ}$ depends on the chemical
cell convention. To avoid confusion, it is strongly recommended to build the
[`System`](@ref) from a chemical cell that follows ITA crystallographic
standards. For example, use a cubic cell to model a fcc or bcc crystal. The
chemical cell convention will be preserved throughout the calculation,
independent of subsequent calls to [`reshape_supercell`](@ref).

Each [`intensities`](@ref) calculator accepts a `measure` argument that
determines how $𝒮^{αβ}(𝐪, ω)$ is contracted:

- [`ssf_trace`](@ref) contracts with $δ^{αβ}$, yielding the trace $\sum_α
  𝒮^{αα}$.
- [`ssf_perp`](@ref) contracts with the polarization factor $δ^{αβ} - q^α
  q^β/q^2$.
- [`ssf_custom`](@ref) accepts a user-provided $𝐪$-dependent contraction
  function. As one option, it may return the full 3×3 matrix.
- [`ssf_custom_bm`](@ref) is similarly flexible, but expects $𝐪$ in the
  Blume-Maleev axis system.

The default form factors are $f_j(𝐪) = 1$. Pass `formfactors` to select the
[`FormFactor`](@ref) appropriate to each magnetic ion. Pass `apply_g=false` to
replace the magnetic moments $\hat{𝐌}$ by the spin operators $\hat{𝐒}$ as the
fundamental observables.

A temperature `kT` passed to [`intensities`](@ref) sets the thermal occupation
of excitations. An energy `kernel` models the instrument resolution function
$K(\omega, \omega')$. Sunny will use it to calculate broadened intensities,

```math
\begin{equation}
I(𝐪, ω) = \int K(ω, ω') 𝒮(𝐪, ω') dω'.
\end{equation}
```

### Unpolarized neutron scattering

Use [`ssf_perp`](@ref) to calculate the differential cross section for an
unpolarized neutron beam. Within the dipole approximation,

```math
\begin{equation}
\frac{d^2 σ(𝐪, ω)}{dω dΩ} = \frac{q_f}{q_i} \left(\frac{γ r_0}{2}\right)^2 N_\mathrm{cells} \sum_{α,β} \left(δ_{α,β} - \frac{q^α q^β}{q^2}\right) 𝒮^{αβ}(𝐪, ω).
\end{equation}
```

Here $𝐪$ and $ω$ denote momentum and energy transferred _to_ the sample. The
prefactor $q_f/q_i$ for $𝐪 = 𝐪_i - 𝐪_f$ is provided experimentally. The
characteristic magnetic scattering length is $γ r_0 / 2 ≈ 2.69×10^{-5}
\mathrm{Å}$, where $γ = 1.913…$ is the neutron gyromagnetic ratio and $r_0$ is
the classical electron radius. Contraction with the polarization factor can be
obtained with a  measure.

### Quantum sum rule

Use [`ssf_trace`](@ref) with `apply_g=false` to validate the quantum sum rule.
For pure spin correlations, the trace of the exact $𝒮^{αβ}(𝐪, ω)$, integrated
over all energies and averaged over momenta, satisfies

```math
\begin{equation}
\left\langle \int_{-∞}^∞ \sum_α 𝒮^{αα}(𝐪, ω) dω \right\rangle_𝐪 = \sum_j s_j (s_j + 1).
\end{equation}
```

The sum runs over sublattices $j$ of the chemical cell, each carrying quantum
spin $s_j$. The momentum average runs over many Brillouin zones, as detailed in
the [derivation](@ref "Derivation of the quantum sum rule"). Sunny's calculators
are only approximate, and may not recover this exact total.


## Dissipation and the retarded susceptibility

Generalizing the structure factor to arbitrary observables $\hat{A}$ and
$\hat{B}$, consider the dynamical correlation

```math
\begin{equation}
C(𝐪, ω) ≡ \frac{1}{2π} \int_{-∞}^{∞} e^{-iωt} ⟨\hat{B}_𝐪^†(0) \hat{A}_𝐪(t)⟩ dt.
\end{equation}
```

The choice $\hat{A} = \hat{M}^β$ and $\hat{B} = \hat{M}^α$ gives
$\mathcal{S}^{αβ}(𝐪, ω)$. Operators evolve in the Heisenberg picture,
$\hat{A}(t) ≡ e^{i t \hat{H}} \hat{A} e^{-i t \hat{H}}$, and $⟨\hat{O}⟩ ≡
\mathrm{tr}\, e^{-β \hat{H}} \hat{O} / \mathcal{Z}$ is the thermal average. Let
$ϵ_μ$ and $|μ⟩$ denote the exact eigenvalues and eigenstates of $\hat{H}$. In
the Lehmann representation,

```math
\begin{equation}
C(𝐪,ω) = \frac{1}{\mathcal{Z}} e^{-β ϵ_ν} δ(ϵ_μ - ϵ_ν - ω) ⟨ν|\hat{B}^†_𝐪|μ⟩⟨μ|\hat{A}_𝐪 |ν⟩,
\end{equation}
```

with implicit summation over $μ$ and $ν$. Each term is a transition $|ν⟩ →
|μ⟩$ that absorbs energy $ω$, weighted by the thermal population of the initial
state. This weight carries a universal dependence on temperature. On the
support of the delta function, $e^{-β ϵ_μ} = e^{-β ω} e^{-β ϵ_ν}$, so the
weight factors as $e^{-β ϵ_ν} = [1 + n_B(ω)] (e^{-β ϵ_ν} - e^{-β ϵ_μ})$, where
$n_B(ω) ≡ 1 / (e^{β ω} - 1)$ is the Bose function. The result is the
fluctuation-dissipation theorem,

```math
\begin{equation}
C(𝐪, ω) = \frac{1}{π} \left[1 + n_B(ω)\right] χ''(𝐪, ω),
\end{equation}
```

where

```math
\begin{equation}
χ''(𝐪, ω) ≡ \frac{π}{\mathcal{Z}} \left(e^{-β ϵ_ν} - e^{-β ϵ_μ}\right) δ(ϵ_μ - ϵ_ν - ω) ⟨ν|\hat{B}^†_𝐪|μ⟩⟨μ|\hat{A}_𝐪 |ν⟩.
\end{equation}
```

The population difference counts absorption minus stimulated emission, i.e., the
energy dissipated in the sample. It makes $χ''$ odd: for Hermitian $\hat{A}$ and
$\hat{B}$, $χ''(-𝐪, -ω) = -χ''(𝐪, ω)^*$. All asymmetry between positive and
negative $ω$ resides in the thermal factor, which satisfies $1 + n_B(-ω) =
-e^{-βω} [1 + n_B(ω)]$. This yields detailed balance, $C(-𝐪, -ω) = e^{-βω}
C(𝐪, ω)^*$. At zero temperature, the thermal factor becomes the step function
$θ(ω)$.

The dissipation is the anti-Hermitian part, $χ'' = (χ - χ^†)/2i$, of the retarded
susceptibility,

```math
\begin{equation}
χ(𝐪, z) = \frac{1}{\mathcal{Z}} \frac{e^{-β ϵ_ν} - e^{-β ϵ_μ}}{ϵ_μ - ϵ_ν - z} ⟨ν|\hat{B}^†_𝐪|μ⟩⟨μ|\hat{A}_𝐪 |ν⟩,
\end{equation}
```

evaluated at $z → ω + i0^+$. In the time domain, the population difference
arises from a commutator,

```math
\begin{equation}
χ(𝐪, z) = i \int_0^∞ e^{izt} ⟨[\hat{B}_𝐪^†(t), \hat{A}_𝐪(0)]⟩ dt.
\end{equation}
```

By the Kubo formula, $χ$ is the linear response of $⟨\hat{B}_𝐪^†⟩$ to a field
that couples to $\hat{A}_𝐪$. Causality makes it analytic in the upper
half-plane, which implies the Kramers-Kronig relation,

```math
\begin{equation}
χ(𝐪, z) = \frac{1}{π} \int_{-∞}^∞ \frac{χ''(𝐪, ω')}{ω' - z} dω'.
\end{equation}
```


## Calculator details

Each of Sunny's calculators estimates $χ''(𝐪, ω)$, multiplies by the thermal
factor $[1 + n_B(ω)]/π$, and convolves with the `kernel`. They differ in how
$χ''$ is obtained.

### Spin wave theory

In linear spin wave theory, quantum spin operators are expressed with
Holstein-Primakoff bosons, and the Hamiltonian is expanded to quadratic order,
i.e., to leading order in inverse powers of the quantum spin-$s$. For systems
constructed with `mode = :SUN`, Sunny automatically switches to a multi-flavor
boson variant of spin wave theory, which captures more single-ion physics.

The population difference $e^{-β ϵ_ν} - e^{-β ϵ_μ}$ in $χ''$ arises from a
commutator, and the commutator of two bosons is a c-number. In linear spin wave
theory, $χ$ is therefore independent of temperature, and the temperature
dependence of the intensity enters only through $1 + n_B(ω)$. The susceptibility
$χ(𝐪, z)$ is a finite sum of magnon poles. The calculators differ in how they
evaluate it.

- [`SpinWaveTheory`](@ref) diagonalizes the quadratic Hamiltonian at each $𝐪$,
  which gives every pole and its residue. [`intensities_bands`](@ref) returns
  these directly, and [`intensities`](@ref) broadens them.
- [`SpinWaveTheorySpiral`](@ref) studies generalized spiral phases, which allow
  for an incommensurate propagation wavevector $𝐤$. The magnons at $𝐪$ combine
  three branches, at $𝐪$ and $𝐪 ± 𝐤$ of the co-rotating frame.
- [`SpinWaveTheoryLarge`](@ref) uses iterated matrix-vector multiplication to
  achieve linear-scaling cost in the magnetic cell size. The Lanczos method
  produces $χ(𝐪, z)$ as a continued fraction, which is evaluated at complex
  energies without resolving individual poles. Iterations stop when the
  estimated error is below `tol` relative to the largest intensity at each $𝐪$.
  It can be useful for studying systems with large magnetic cells, including
  systems with long-wavelength structures, or systems with quenched chemical
  disorder.

Spin wave calculators also accept a regulator $η > 0$, which evaluates $χ$ at
$ω+iη$ rather than $ω + i0^+$. Finite $η$ defines a half width for intrinsic
Lorentzian broadening and damps intensity of soft modes with excitation energy
$ϵ ≲ η$.

### Classical spin dynamics

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

Like in spin wave theory, intensities are extracted with
[`intensities`](@ref). Since classical dynamics simulation take place on a
finite lattice, intensities are only available at a discrete grid of wave
vectors. In reciprocal lattice units, available grid points are $𝐪 =
[\frac{n_1}{L_1}, \frac{n_2}{L_2}, \frac{n_3}{L_3}]$, where $n_i$ runs from
$(\frac{-L_i}{2}+1)$ to $\frac{L_i}{2}$ and $L_i$ is the linear dimension of the
lattice used in the calculation. (An internal function
`Sunny.available_wave_vectors` provides access to this grid.) By default
[`intensities`](@ref) will adjust each wavevector $𝐪$ to the nearest available
grid point.

Similarly, the resolution in `energies` is controlled by the dynamical
trajectory length in real-time. Because the dynamical trajectory is not periodic
in time, some blurring between neighboring energy bins is unavoidable. Sunny's
algorithm for estimating the structure factor from real-time dynamics is
[specified
here](https://github.com/SunnySuite/Sunny.jl/pull/246#issuecomment-2119294846).

Classical trajectories sample the classical correlation $C_\mathrm{cl}$, for
which the fluctuation-dissipation theorem reads $C_\mathrm{cl}(𝐪, ω) = (kT/πω)
χ''(𝐪, ω)$. That is, the thermal factor $1 + n_B(ω)$ is replaced by its
high-temperature limit $kT/ω$. Given `kT`, [`intensities`](@ref) multiplies
$C_\mathrm{cl}$ by $βω [1 + n_B(ω)]$, which is equivalent to extracting $χ''$
and applying the quantum thermal factor. This classical-to-quantum rescaling
makes the result comparable with spin wave theory. If `kT = nothing`, the factor
is omitted, and intensities follow the classical Boltzmann distribution.

### Static correlations

Use [`intensities_static`](@ref) to calculate $𝒮(𝐪)$, i.e., correlations that
are "instantaneous" in real-time. Mathematically, $𝒮(𝐪)$ denotes an integral
of the dynamical structure factor $𝒮(𝐪, ω)$ over all energies $ω$. In
[`SpinWaveTheory`](@ref), the energy integral becomes a discrete sum over bands.
In [`SampledCorrelations`](@ref), a classical-to-quantum correction factor will
be applied within [`intensities`](@ref) prior to energy integration.

Sunny also supports a mechanism to calculate static correlations without any
spin dynamics. To collect such statistics, construct a
`SampledCorrelationsStatic` object. In this case, `intensities_static` will
return static correlations sampled from the classical Boltzmann distribution.
This dynamics-free approach is faster, but may miss important features that
derive from the quantum mechanical excitation spectrum.


## Derivation of the quantum sum rule

Integrating over all energies yields the equal-time correlation,

```math
\begin{equation}
\int_{-∞}^∞ \mathcal{S}^{αα}(𝐪, ω) dω = ⟨\hat{M}_𝐪^{α†} \hat{M}_𝐪^α⟩,
\end{equation}
```

with summation on $α$ implied. For pure spin correlations, $\hat{𝐌}_𝐪$ is
replaced by $\hat{𝐒}_𝐪 = \sum_{𝐦,j} e^{i 𝐪⋅𝐫_{𝐦,j}} \hat{𝐒}_{𝐦,j}$, and
the right-hand side becomes a double sum over sites,

```math
\begin{equation}
⟨\hat{𝐒}_𝐪^† ⋅ \hat{𝐒}_𝐪⟩ = \sum_{𝐦,j} \sum_{𝐦',j'} e^{i 𝐪⋅(𝐫_{𝐦',j'} - 𝐫_{𝐦,j})} ⟨\hat{𝐒}_{𝐦,j} ⋅ \hat{𝐒}_{𝐦',j'}⟩.
\end{equation}
```

Averaging over momenta $𝐪$ eliminates every term with $𝐫_{𝐦',j'} ≠
𝐫_{𝐦,j}$. Within a single Brillouin zone, the phases cancel between distinct
cells $𝐦 ≠ 𝐦'$. Between distinct sublattices $j ≠ j'$ of the same cell, the
average must extend over many Brillouin zones. Only the on-site terms survive,
and each contributes $⟨\hat{𝐒}^2⟩ = s_j(s_j+1)$. There are $N_\mathrm{cells}$
such terms for each sublattice $j$, and the normalization of $𝒮$ divides by
$N_\mathrm{cells}$. This recovers the sum rule stated
[above](@ref "Quantum sum rule").

The calculators recover this total only in part.

- Spin wave theory omits the elastic Bragg weight from its band intensities.
  For a polarized ferromagnet, the inelastic weight per spin is exactly $s$,
  and adding the elastic weight $s^2$ recovers $s(s+1)$. For a Néel
  antiferromagnet, the inelastic weight is $s + δ$, where $δ$ is the zero-point
  reduction of the ordered moment. The sum rule then holds only to leading
  order in $1/s$.
- Classical spin dynamics conserves the classical magnitude $|𝐒|^2 = s^2$,
  which is the total obtained with `kT = nothing`. With finite `kT`, the
  classical-to-quantum rescaling of [`intensities`](@ref) modifies the total.
  For static correlations, [`set_spin_rescaling_for_static_sum_rule!`](@ref)
  enforces the quantum value $s(s+1)$.
