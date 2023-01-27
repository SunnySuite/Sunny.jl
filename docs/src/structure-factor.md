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
"Real life" use cases can be found in our tutorials and detailed function
information is available in the Library API.

## Basic Usage

The basic data type for calculating, storing and retrieving structure factor
data is `StructureFactor`. Rather than creating a `StructureFactor` directly,
one should call either `DynamicStructureFactor`, for $𝒮^{αβ}(𝐪,ω)$, or
`StaticStructureFactor`, for $𝒮^{αβ}(𝐪)$. These functions will configure and
return an appropriate `StructureFactor`.

### Calculating a dynamical stucture factor

Calling `DynamicStructureFactor(sys; Δt, ωmax, nω)` will create a
`StructureFactor` for the user and calculate an initial sample. There are three
keywords that must be specified. These keywords will determine the dynamics
used to calculate the sample and, consequently, the $ω$ information that will be
available after the calculation has completed.

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

Structure factor data is calculated from classical dynamics using a Monte Carlo
approach. Each sample of $𝒮^{αβ}(𝐪,ω)$ is generated from a trajectory, the
initial condition for which must be a sample spin configuration from the
equilibrium distribution at a particular temperature. Therefore it is important
that the spin configuration in your `sys` represent such a sample prior to
calling `DynamicStructureFactor`. In other words, it is essential that `sys` be
properly thermalized before initiating the calculation. One approach to
accomplishing this is to use a [`LangevinSampler`](@ref).

Additional samples can be accumulated into a `StructureFactor` by calling
[`add_sample!(structure_factor, sys)`](@ref). Naturally, it is important that
the spin configuration in `sys` represent a new equilibrium sample before
calling `add_sample!`.

The outline of typical use case might look like this:
```
# Thermalize a `System`, `sys`, and set up a `LangevinSampler`, `sampler`
# prior to steps below. 

# Make a `StructureFactor` and calculate an initial sample
sf = DynamicStructureFactor(sys; Δt=0.05, ωmax=10.0, nω=100) 

# Add additional samples
for _ in 1:nsamples
   sample!(sys, sampler)   # Update spins in `sys` to generate a new initial condition
   add_sample!(sf, sys)    # Use spins to calculate and accumulate new sample of 𝒮(𝐪,ω)
end
```

The calculation may be configured in a number of ways; see the
[`DynamicStructureFactor`](@ref) documentation for a list of all keywords.


### Calculating a static structure factor

Sunny provides two methods for calculating static structure factors,
$𝒮^{αβ}(𝐪)$. The first involves calculating spin-spin correlations at single
instances of time. The second involves calculating a dynamic structure factor
first and integrating out the $ω$ information. The advantage of the latter
approach is that it enables application of an $ω$-dependent classical-to-quantum
rescaling of structure factor intensities, a method that should be preferred
whenever comparing results to experimental data or spin wave calculations. A
disadvantage of this approach is that it is computationally more expensive.
There are also many cases when it is not straightforward to calculate a
meaningful dynamics, as when working with Ising spins. In this section we will
discuss how to calculate static structure factors from static spin
configurations. Information about calculating static data from a dynamic
structure factor can be found in the following section.

The basic usage for the static case is very similar to the dynamic case, except
one calls `StaticStructureFactor(sys)` instead of `DynamicStructureFactor`. Note
that there are no required keywords as there is no need to specify any dynamics.
`StaticStructureFactor` will immediately calculate a sample of $𝒮(𝐪)$ using
the spin configuration contained in `sys`. It is therefore important that 
`sys` be properly thermalized before calling this function. Additional samples
may be added with `add_sample!(sf, sys)`, just as was done in the dynamic case.
As was true there, it is important to ensure that the spins in `sys` represents
a new equilibrium sample before calling `add_sample!`.

### Extracting information from structure factors

The basic function for extracting information from a dynamic `StructureFactor`
at a particular wave vector, $𝐪$, is [`intensities`](@ref). It takes a
`StructureFactor`, a list of wave vectors, and a contraction mode. For example,
`intensities(sf, [[0.0, 0.5, 0.5]], :trace)` will calculate intensities for the
wavevector $𝐪 = (𝐛_2 + 𝐛_3)/2$. The option `:trace` will contract spin
indices, returning $𝒮^{αα}(𝐪,ω)$. The option `:perp` will instead perform a
contraction that includes polarization corrections. The option `:full` will
return data for the full tensor $𝒮^{αβ}(𝐪,ω)$. `intensities` returns a list of
`nω` elements. The corresponding $ω$ values are given by `ωvals(sf)`, where `sf`
is the `StructureFactor`.

Since Sunny currently only calculates the structure factor on a finite lattice,
it is important to realize that exact information is only available at a
discrete set of wave vectors. Specifically, for each axis index $i$, we will get
information at $q_i = \frac{n}{L_i}$, where $n$ runs from $(\frac{-L_i}{2}+1)$
to $\frac{L_i}{2}$ and $L_i$ is the linear dimension of the lattice used for the
calculation. If you request a wave vector that does not fall into this set,
Sunny will automatically round to the nearest $𝐪$ that is available. If
`intensities` is given the keyword argument `interpolation=:linear`, Sunny will
use trilinear interpolation to determine the results at the requested wave
vector. 

To retrieve the intensities at all wave vectors for which there is exact data,
first call the function [`all_exact_wave_vectors`](@ref) to generate a list of
`qs`. This takes an optional keyword argument `bzsize`, which must be given a
tuple of three integers specifying the number of Brillouin zones to calculate,
e.g., `bzsize=(2,2,2)`. The resulting list of wave vectors may then be passed to
`intensities`.

The convenience function [`connected_path`](@ref) returns a list of wavevectors
sampled along a path that connects specified $𝐪$ points. This list can be used
as an input to `intensities`.

A number of keyword arguments are available which modify the calculation of
structure factor intensity. See the documentation of [`intensities`](@ref) for a
full list. It is generally recommended to provide a value to `kT` corresponding
to the temperature of sampled configurations. Given `kT`, Sunny will apply an
energy- and temperature-dependent classical-to-quantum rescaling of intensities. 

To retrieve intensity data from a static structure factor, use
[`static_intensities`](@ref), which shares keyword arguments with
[`intensities`](@ref). This function may also be used to calculate static
information from a dynamical structure factor. Note that it is important to
supply a value to `kT` to reap the benefits of this approach over simply
calculating a static structure factor at the outset. 