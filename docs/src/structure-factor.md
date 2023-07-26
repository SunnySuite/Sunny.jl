# Structure Factor Calculations

## Overview
The dynamical structure factor is of fundamental importance for characterizing a
magnetic system, and facilitates quantitative comparison between theory and
experimental scattering data.

Consider, for example, a two-point dynamical spin correlation function,
$⟨s^α(𝐱+Δ𝐱, t+Δt) s^β(𝐱, t)⟩$. Here $s^α(𝐱, t)$ represents the time dynamics
of a spin dipole component $α$ at position $𝐱$, and brackets represent an
average over equilibrium initial conditions and over $(𝐱, t)$. The dynamical
structure factor is defined as the Fourier transform of this two-point
correlation in both space and time, up to an overall scaling factor. Using the
convolution theorem, the result is,

$$𝒮^{αβ}(𝐪, ω) = \frac{1}{V} ⟨ŝ^α(𝐪, ω)^\ast ŝ^β(𝐪, ω) ⟩,$$

with $V$ the system volume. We will restrict attention to lattice systems with
periodic boundaries.

Consider a crystal unit cell defined by three lattice vectors $𝐚_1, 𝐚_2,
𝐚_3$, and linear system sizes $L_1, L_2, L_3$ measured in unit cells. The
allowed momentum vectors take on discrete values $𝐪 = \sum_{α=1}^{3} m_α 𝐛_α /
L_α$, where $m_α$ are an integers and the reciprocal lattice vectors $𝐛_α$ are
defined to satisfy $𝐚_α ⋅ 𝐛_β = 2π δ_{α,β}$. For a Bravais lattice, $𝐪$ will
be periodic in the first Brillouin zone, i.e., under any shift $𝐪 → 𝐪 ± 𝐛_α$.
More generally, consider a non-Bravais lattice such that each unit cell may
contain multiple spins. By partitioning spins $s_j(𝐱,t)$ according to their
sublattice index $j$, the relevant momenta $𝐪$ remain discretized as above, but
now periodicity in the first Brillouin zone is lost. The structure factor may be
written as a phase-average over the displacements between sublattices
$𝐫_{j,k}$,

$$𝒮^{αβ}(𝐪, ω) = ∑_{j,k} e^{i 𝐫_{j,k} ⋅ 𝐪} 𝒮̃^{αβ}_{j,k}(𝐪, ω) ⟩,$$

From a theoretical perspective, the quantity

$$𝒮̃^{αβ}_{j,k}(𝐪, ω) = \frac{1}{V} ⟨ŝ_j^α(𝐪, ω)^\ast ŝ_k^β(𝐪, ω)⟩$$

is fundamental. For each sublattice $j$, the data $ŝ_j^α(𝐪, ω)$ can be
efficiently obtained by fast Fourier tranformation of a real space configuration
$s_j^α(𝐱, t)$. Internally, Sunny will calculate and store the discrete
$𝒮̃^{αβ}_{j,k}(𝐪, ω)$ correlation data, and use this to construct
$𝒮^{αβ}(𝐪,ω)$ intensities that can be compared with experiment.

Calculating this structure factor involves several steps, with various possible
settings. Sunny provides a number of tools to facilitate this calculation and to
extract information from the results. These tools are briefly outlined below.
Please see the Examples for a "real life" use case. Detailed function
information is available in the Library API.

## Estimating stucture factors with classical dynamics

The basic approach to estimating structure factor information using classical
dynamics relies on the generation of spin-spin correlation data from dynamical
trajectories. This is fundamentally a Monte Carlo calculation, as the
trajectories must be started from an initial spin configuration that is sampled
from thermal equilibrium. (Note that it is not possible to estimate a true T=0
dynamical structure factor using this approach, but the temperature may be very
small.) Samples are accumulated into a `SampledCorrelations`, from which
intensity information may be extracted. The user does not typically build their
own `SampledCorrelations`, but instead initializes one using either
`dynamical_correlations` or `instant_correlations`, as described below.

### Estimating a dynamical structure factor: ``𝒮(𝐪,ω)``

A `SampledCorrelations` for estimating the dynamical structure factor,
$𝒮^{αβ}(𝐪,ω)$, may be created by calling [`dynamical_correlations`](@ref). This
requires three keyword arguments. These will determine the dynamics used to
calculate samples and, consequently, the $ω$ information that will be available. 

1. `Δt`: Determines the step size used for simulating the dynamics. A smaller
   number will require proportionally more calculation time. While a smaller
   `Δt` will enable the resolution of higher energies, `Δt` is typically
   selected to ensure numerical stability rather than to maximize the largest
   $ω$ value. A safe choice is to use the smaller value of `Δt = 0.1/(J* S^2)`
   or `Δt = 0.1/(D * S)`, where `S` is magnetic moment of the largest local spin
   (as specified in [`SpinInfo`](@ref)), `J` is the parameter governing the
   largest bilinear interaction (e.g. exchange), and `D` is the parameter
   governing the largest single-site term of the Hamiltonian (e.g., anisotropy
   or Zeeman term).
2. `ωmax`: Sets the maximum resolved energy. Note that this is not independent
   of `Δt`. If `ωmax` too large, Sunny will throw an error and ask you to choose
   a smaller `Δt`. 
3. `nω`: Determines the number of energy bins to resolve. A larger number will
   require more calculation time.

A sample may be added by calling `add_sample!(sc, sys)`. The input `sys` must be
a spin configuration in good thermal equilibrium, e.g., using the continuous
[`Langevin`](@ref) dynamics or using single spin flip trials with
[`LocalSampler`](@ref). The statistical quality of the $𝒮^{αβ}(𝐪,ω)$ can be
improved by repeatedly generating decorrelated spin configurations in `sys` and
calling `add_sample!` on each configuration.

The outline of typical use case might look like this:
```
# Make a `SampledCorrelations`
sc = dynamical_correlations(sys; Δt=0.05, ωmax=10.0, nω=100) 

# Add samples
for _ in 1:nsamples
   decorrelate_system(sys) # Perform some type of Monte Carlo simulation
   add_sample!(sc, sys)    # Use spins to calculate trajectory and accumulate new sample of 𝒮(𝐪,ω)
end
```
The calculation may be configured in a number of ways; see the
[`dynamical_correlations`](@ref) documentation for a list of all keywords.


### Estimating an instantaneous ("static") structure factor: ``𝒮(𝐪)``

Sunny provides two methods for calculating instantaneous, or static, structure
factors: $𝒮^{αβ}(𝐪)$. The first involves calculating spatial spin-spin
correlations at single time slices. The second involves calculating a dynamic
structure factor first and integrating out the $ω$ information. The advantage of
the latter approach is that it enables application of an $ω$-dependent
classical-to-quantum rescaling of structure factor intensities, a method that
should be preferred whenever comparing results to experimental data or spin wave
calculations. A disadvantage of this approach is that it is computationally more
expensive. There are also many cases when it is not straightforward to calculate
a meaningful dynamics, as when working with Ising spins. In this section we will
discuss how to calculate instantaneous structure factors from static spin
configurations. Information about calculating instantaneous data from a
dynamical correlations can be found in the following section.

The basic usage for the instantaneous case is very similar to the dynamic case,
except one calls [`instant_correlations`](@ref) instead of
`dynamical_correlations` to configure a `SampledCorrelations`. Note that there
are no required keywords as there is no need to specify any dynamics.
`instant_correlations` will return a `SampledCorrelations` containing no data.
Samples may be added by calling `add_sample!(sc, sys)`, where `sc` is the
`SampledCorrelations`. When performing a finite-temperature calculation, it is
important to ensure that the spin configuration in the `sys` represents a good
equilibrium sample, as in the dynamical case. Note, however, that we recommend
calculating instantaneous correlations at finite temperature calculations by
using full dynamics (i.e., using `dynamical_correlations`) and then integrating
out the energy axis. An approach to doing this is described in the next section.

### Extracting information from correlation data 

The basic function for extracting information from a `SampledCorrelations`
at a particular wave vector, $𝐪$, is [`intensities_interpolated`](@ref). It takes a
`SampledCorrelations` and a list of wave vectors. For example,
`intensities_interpolated(sf, [[0.0, 0.5, 0.5]])` will calculate intensities for the
wavevector $𝐪 = (𝐛_2 + 𝐛_3)/2$. The keyword argument `formula` can be used to
specify an [`intensity_formula`](@ref) for greater control over the intensity calculation.
The default formula performs a contraction of $𝒮^{αβ}(𝐪,ω)$ that includes
polarization corrections. `intensities_interpolated`
returns a list of `nω` elements at each wavevector. The corresponding $ω$ values can be retrieved
by calling [`ωs`](@ref) on `sf`.

Since Sunny currently only calculates the structure factor on a finite lattice,
it is important to realize that exact information is only available at a
discrete set of wave vectors. Specifically, for each axis index $i$, we will get
information at $q_i = \frac{n}{L_i}$, where $n$ runs from $(\frac{-L_i}{2}+1)$
to $\frac{L_i}{2}$ and $L_i$ is the linear dimension of the lattice used for the
calculation. If you request a wave vector that does not fall into this set,
Sunny will automatically round to the nearest $𝐪$ that is available. If
`intensities_interpolated` is given the keyword argument
`interpolation=:linear`, Sunny will use trilinear interpolation to determine a
result at the requested wave vector. 

To retrieve the intensities at all wave vectors for which there is exact data,
first call the function [`all_exact_wave_vectors`](@ref) to generate a list of
`qs`. This takes an optional keyword argument `bzsize`, which must be given a
tuple of three integers specifying the number of Brillouin zones to calculate,
e.g., `bzsize=(2,2,2)`. The resulting list of wave vectors may then be passed to
`intensities_interpolated`.

Alternatively, [`intensities_binned`](@ref) can be used to place the exact data
into histogram bins for comparison with experiment.

The convenience function [`connected_path`](@ref) returns a list of wavevectors
sampled along a path that connects specified $𝐪$ points. This list can be used
as an input to `intensities`. Another convenience method,
[`spherical_shell`](@ref) will provide a list of wave vectors on a sphere of a
specified radius. This is useful for powder averaging. 

A number of arguments for [`intensity_formula`](@ref) are available which
modify the calculation of structure factor intensity. It is generally recommended
to provide a value of `kT` corresponding to the temperature of sampled configurations.
Given `kT`, Sunny will include an energy- and temperature-dependent classical-to-quantum 
rescaling of intensities in the formula.

To retrieve intensity data from a instantaneous structure factor, use
[`instant_intensities_interpolated`](@ref), which accepts similar arguments to
`intensities_interpolated`. This function may also be used to calculate
instantaneous information from a dynamical structure factor, i.e. from a
`SampledCorrelations` created with `dynamical_correlations`. Note that it is
important to supply a value to `kT` to reap the benefits of this approach over
simply calculating a static structure factor at the outset. 