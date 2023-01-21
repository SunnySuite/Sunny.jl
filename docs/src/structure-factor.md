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
$𝒮̃^{αβ}_{j,k}(𝐪, ω)$ correlation data, and use this to construct $𝒮^{αβ}(𝐪,
ω)$ intensities that can be compared with experiment.

Calculating this structure factor involves several steps, with various possible
settings. Sunny provides a number of tools to facilitate this calculation and to
extract information from the results. These tools are briefly outlined below.
"Real life" use cases can be found in our tutorials and detailed function
information is available in the Library API.

## Basic Usage

### Calculating a dynamical stucture factor

The basic function for calculating dynamical structure factors is
[`calculate_structure_factor`](@ref). The steps for using it effectively are the
following:

1. Build a [`System`](@ref) and ensure that it is properly equilibrated at the
   temperature you wish to study. For example, if the `System` is in a ground
   state, one could use a [`LangevinHeunP`](@ref) integrator to thermalize it.
2. Set up a sampler that will generate decorrelated samples of spin
   configurations at the desired temperature, for example, by using a
   [`LangevinSampler`](@ref).
3. Call `calculate_structure_factor(sys, sampler; kwargs...)`, which will return
   a `StructureFactor`, containing all $𝒮̃^{αβ}_{jk}(𝐪, ω)$ data.

The calculation can be configured in a number of ways; see
[`calculate_structure_factor`](@ref) documentation for a list of all keywords.
In particular, note that an argument `nω` greater than one must be specified to
get a dynamical structure factor.

### Extracting information

The basic function for extracting information from a `StructureFactor` at a
particular wave vector, $𝐪$, is [`get_intensities`](@ref). It takes a
`StructureFactor`, a list of wave vectors, and a contraction mode. For example,
`get_intensities(sf, [[0.0, 0.5, 0.5]], :trace)` will calculate intensities for
the wavevector $𝐪 = (𝐛_2 + 𝐛_3)/2$. The option `:trace` will contract spin
indices, returning $𝒮^{αα}(𝐪,ω)$. The option `:perp` will instead perform a
contraction that includes polarization corrections. The option `:full` will
return data for the full tensor $𝒮^{αβ}(𝐪,ω)$. `get_intensities` returns a
list of `nω` elements. The corresponding $ω$ values are given by `ωvals(sf)`,
where `sf` is the `StructureFactor`.

The convenience function [`connected_path`](@ref) returns a list of wavevectors
sampled along a path that connects specified $𝐪$ points. This list can be used
as an input to `get_intensities`.

Since Sunny currently only calculates the structure factor on a finite lattice,
it is important to realize that exact information is only available at a
discrete set of wave vectors. Specifically, for each axis index $i$, we will get
information at $q_i = \frac{n}{L_i}$, where $n$ runs from
$(\frac{-L_i}{2}+1)$ to $\frac{L_i}{2}$ and $L_i$ is the linear dimension of
the lattice used for the calculation. If you request a wave vector that does not
fall in this set, Sunny will automatically round to the nearest $𝐪$ that is
available. If `get_intensities` is given the keyword argument
`interpolation=:linear`, Sunny will use trilinear interpolation to determine the
results at the requested wave vector. 

To retrieve the intensities at all wave vectors for which there is exact data,
one can use the function [`intensity_grid`](@ref). This takes an optional
keyword argument `bzsize`, which must be given a tuple of three integers
specifying the number of Brillouin zones to calculate, e.g., `bzsize=(2,2,2)`.

Many keyword arguments are available which modify the calculation of structure
factor intensity. See the documentation of [`get_intensities`](@ref) for a full
list. It is generally recommended to provide a value to `kT` corresponding to
the temperature of sampled configurations. Given `kT`, Sunny will apply a
classical-to-quantum rescaling of the energy intensities. 

### Static structure factors

A static structure will be calculated if the `nω` keyword of
`calculate_structure_factor` or `StructureFactor` is left at its default value
of 1. Static structure factors may also be calculated from a dynamical structure
factor simply by summing over all the energies (i.e., the $ω$-axis) provided
by `get_intensities`. We recommend calculating static structure factors in this
way in most cases (though it is of course much more expensive). The
static-from-dynamic approach makes it possible to apply the classical-to-quantum
intensity rescaling, which is energy dependent. Sunny provides the function
[`get_static_intensities`](@ref), which will perform the summation for you.