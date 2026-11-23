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
cell convention, which is set by [`Crystal`](@ref). To avoid confusion, Sunny
will warn if the chemical cell deviates from ITA crystallographic standards. For
example, to model fcc or bcc, one should begin with a cubic cell. Any subsequent
calls to [`reshape_supercell`](@ref) (e.g., to build a primitive cell) will
preserve the original chemical cell convention.

Each intensities calculator accepts a `measure` argument that determines how
$𝒮^{αβ}(𝐪, ω)$ is contracted:

- [`ssf_trace`](@ref) contracts with $δ^{αβ}$, yielding the trace of $𝒮$.
- [`ssf_perp`](@ref) contracts with the projector $δ^{αβ} - q^α q^β/q^2$, as
  relevant to unpolarized neutron scattering.
- [`ssf_custom`](@ref) accepts a user-provided $𝐪$-dependent contraction
  function. It may return off-diagonal elements or the full 3×3 matrix itself.
- [`ssf_custom_bm`](@ref) is similarly flexible, but expects $𝐪$ in the
  Blume-Maleev axis system.

The default form factors are $f_j(𝐪) = 1$, which model purely localized
dipoles. Pass `formfactors` to select the [`FormFactor`](@ref) appropriate to
each magnetic ion. Pass `apply_g=false` to replace the magnetic moments
$\hat{𝐌}$ by the spin operators $\hat{𝐒}$ as the fundamental observables.

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
\frac{d^2 σ(𝐪, ω)}{dω dΩ} = \frac{q_f}{q_i} \left(\frac{γ r_0}{2}\right)^2 N_\mathrm{cells} \sum_{α,β} \left(δ^{αβ} - \frac{q^α q^β}{q^2}\right) 𝒮^{αβ}(𝐪, ω).
\end{equation}
```

Here $𝐪$ and $ω$ denote momentum and energy transferred _to_ the sample. The
prefactor $q_f/q_i$ for $𝐪 = 𝐪_i - 𝐪_f$ is provided experimentally. The
characteristic magnetic scattering length is $γ r_0 / 2 ≈ 2.69×10^{-5}
\mathrm{Å}$, where $γ = 1.913…$ is the neutron gyromagnetic ratio and $r_0$ is
the classical electron radius.

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

A general dynamical correlation is,

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

Each of Sunny's dynamical calculators estimates $χ''(𝐪, ω)$, multiplies by the
thermal factor $[1 + n_B(ω)]/π$, and convolves with the `kernel`. They differ in
how $χ''$ is obtained.

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

Spin wave calculators also accept a regulator $η > 0$, which evaluates $χ$ at $z
= ω + iη$ rather than $ω + i0^+$. By the Kramers-Kronig relation, the
anti-Hermitian part of $χ(z)$ is the convolution of $χ''$ with a Lorentzian of
half width $η$. Because $χ''$ is odd, the broadened poles at $±ϵ$ partially
cancel near $ω = 0$, which damps the intensity of soft modes with $ϵ ≲ η$.

### Classical spin dynamics

[`SampledCorrelations`](@ref) estimates $𝒮(𝐪, ω)$ at finite temperature from
classical spin dynamics (Landau-Lifshitz, or its SU($N$) generalization). Each
call to [`add_sample!`](@ref) collects statistics from a dynamical trajectory
initialized to the current spin configuration. That configuration should be
drawn from thermal equilibrium, e.g., via [`Langevin`](@ref) sampling.
Statistical error decreases with the number of decorrelated samples.

```julia
sc = SampledCorrelations(sys; measure, dt, energies=range(0, 10, 100))
for _ in 1:nsamples
    for _ in 1:ndecorr
        step!(sys, langevin)
    end
    add_sample!(sc, sys)
end
res = intensities(sc, qpts; energies=:available, kT)
```

The `energies` must be a uniform range starting at zero. Its spacing sets the
trajectory length, and its maximum sets the sampling interval. A good choice of
`dt` can be estimated from [`suggest_timestep`](@ref). Trajectories are finite,
so some blurring between neighboring energy bins is unavoidable. Details of the
estimator are [specified
here](https://github.com/SunnySuite/Sunny.jl/pull/246#issuecomment-2119294846).

The simulation volume is periodic, so intensities are available only on a
discrete grid of wavevectors. For a system of $L_1 × L_2 × L_3$ chemical cells,
these are $q_i = n_i / L_i$ in reciprocal lattice units. The
[`intensities`](@ref) function rounds each requested $𝐪$ to the nearest grid
point.

Classical trajectories sample the classical correlation, for which the
fluctuation-dissipation theorem reads $C_\mathrm{cl}(𝐪, ω) = (kT/πω) χ''(𝐪,
ω)$. The quantum thermal factor $1 + n_B(ω)$ is replaced by its high-temperature
limit $kT/ω$. Given `kT`, [`intensities`](@ref) multiplies by $βω [1 + n_B(ω)]$,
which restores the quantum factor and makes the result comparable with spin wave
theory. With `kT = nothing`, the classical correlations are reported unchanged.

### Static correlations

The static structure factor $𝒮(𝐪) = \int 𝒮(𝐪, ω) dω$ is the equal-time
correlation. Sunny's [`intensities_static`](@ref) obtains it in one of two ways.

Frequency-dependent calculators evaluate the $ω$ integral, with optional
restriction to energy `bounds`. For [`SpinWaveTheory`](@ref) and
[`SpinWaveTheorySpiral`](@ref), the integral becomes a sum over bands, and omits
the elastic Bragg weight. For [`SampledCorrelations`](@ref), it is a sum over
the energy grid, which retains the elastic line at $ω = 0$. With finite `kT`,
the sum is taken after the classical-to-quantum factor. With `kT = nothing`, it
gives the classical equal-time correlations.

The [`SampledCorrelationsStatic`](@ref) calculator skips the dynamics and
instead samples $𝒮(𝐪)$ directly from the _classical_ Boltzmann distribution.
This accelerates the calculation, but reflects classical rather than quantum
thermal statistics. It is most accurate at higher temperatures, where the
fluctuations become effectively classical. The [`SCGA`](@ref) calculator runs
even faster, and makes the additional approximation that spin fluctuations are
effectively Gaussian. SCGA can often be quite good when deep in the paramagnetic
phase. For both calculators, it is recommended to first call
[`set_spin_rescaling_for_static_sum_rule!`](@ref).

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
