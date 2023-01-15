# Structure Factor Calculations

To be written...

This page gives information on the static and dynamical spin structure factors, how to use Sunny's high and low-level interfaces for computing it, and what is happening behind the scenes in these functions!

The central type implementing all of the computation behind the scenes is
[`StructureFactor`](@ref).

## Background

The structure factor is one representation in which to examine how spins are correlated within spin
configurations sampled from the thermal spin distribution defined by the system's Hamiltonian.
Specifically, we will write our spin degrees of freedom as ``S^α_j(𝐫, t)``, where
``𝐫 = n_a 𝐚 + n_b 𝐛 + n_c 𝐜`` is the coordinate of the unit cell, ``t`` is the time during some
evolved dynamics, ``j`` is the index into the basis/sublattice within the unit cell,
and ``α = \{x,y,z\}`` is the spin-component index.

Spin-spin correlations in real space and time can be characterized by the two-point correlation
function:

```math
C^{αβ}_{jk}(𝐫, t) = ⟨S^α_j(𝐫_0, t_0) S^β_k(𝐫_0 + 𝐫, t_0 + t)⟩_{𝐫_0, t_0}
```
where ``⟨⋅⟩_{𝐫_0, t_0}`` means we are taking a thermal average over different spin configurations,
as well as an average over reference positions ``𝐫_0`` in the lattice and times ``t_0`` in the
dynamics. Colloquially, this function is telling us "how much is the ``S^α`` component
of one spin on sublattice ``j`` correlated with the ``S^β`` component of another spin on
sublattice ``k`` which is displaced in position and time by ``(𝐫, t)``?".

The full _dynamic structure factor_ is the Fourier transform of the two-point correlation function.

```math
𝒮^{αβ}_{jk}(𝐪, ω) = \frac{1}{\sqrt{2π}} \sum_𝐫 \int dt e^{-i (𝐪 ⋅ 𝐫 + ωt)} C^{αβ}_{jk}(𝐫, t)
```

This is the quantity which the structure factor module computes. By explicitly performing the
spatial/time averages in our definition of the two-point correlation function, we can obtain
an alternate, more easily calculable form for the dynamic structure factor:

```math
𝒮^{αβ}_{jk}(𝐪, ω) = ⟨S^α_j(𝐪, ω) S^β_k(𝐪, ω)^∗⟩
```
where ``^∗`` refers to complex conjugation. This provides an easy direct route to calculating
the dynamic structure factor:

1. Obtain a bunch of thermal spin configurations ``S^α_j(𝐫)``
2. Using these as initial conditions, time evolve them all forward using Landau-Lifshitz
    dynamics to obtain ``S^α_j(𝐫, t)``.
3. Discrete Fourier transform them all to obtain ``S^α_j(𝐪, ω)``
4. Perform a complex-conjugated outer product to obtain a contribution
   ``S^α_j(𝐪, ω)S^β_k(𝐪, ω)^∗``
5. Average across all sampled spin configurations

Note that in typical presentations, the basis indices are not present as they are included
in the sum/integral over position. However, because spin simulations can resolve basis-dependent
correlations, we may as well keep them around for now. Neutron scattering experiments, however,
cannot resolve basis-dependent correlations, instead seeing only:

```math
𝒮^{αβ}(𝐪, ω) = \sum_{j,k=1}^{B} e^{-i 𝐪 ⋅ (𝐝_j - 𝐝_k)} 𝒮^{αβ}_{jk}(𝐪, ω)
```

where ``B`` is the number of basis sites within the unit cell and ``𝐝_j`` are the basis vectors.
Within this page, we will refer to going from the full structure factor to this reduced form as
performing the "phase-weighted sum" over basis sites.

The _static structure factor_ is the spatial Fourier transform of the equal-time correlation
function.

```math
𝒮^{αβ}_{jk}(𝐪) = \sum_𝐫 e^{-i𝐪 ⋅ 𝐫} C^{αβ}_{jk}(𝐫, 0)
               = \frac{1}{\sqrt{2π}} \int dω 𝒮^{αβ}_{jk}(𝐪, ω)
```

For both of these structure factors, neutron scattering experiments also do not resolve individual
spin components. Instead, the observed scattering intensity is proportional to the result
of applying the neutron dipole factor:

```math
𝒮(𝐪, ω) = ∑_{αβ} (δ_{αβ} - 𝐪̂_α 𝐪̂_β) 𝒮^{αβ}(𝐪, ω)
```

## High-level functions

Sunny exposes one main high-level function which performs the entirety of the steps (1)--(5)
outlined above for you: [`dynamic_structure_factor`](@ref). The documentation on that
function provides a relatively in-depth explanation of all of the arguments.

A helper function [`static_structure_factor`](@ref) also exists, which computes the
static structure factor simply by calling `dynamic_structure_factor` with `num_meas=1`.

## Manual incremental updates

If you are writing the lower-level simulation loop yourself, or have a stack of spin configurations on-hand that you want to compute the structure factor from, there is also an additional direct interface.

If you have a stack of snapshots on hand, then you can directly use them to 
construct a `StructureFactor`. A "stack of snapshots" can either be
represented as a `Vector{SpinSystem}` (all of which have the same underlying
lattice), or a `Vector{Array{Vec3, 4}}` along with a `Crystal` defining
the geometry.

```julia
sf = StructureFactor(spin_snaps; )
```

All of the Fourier transforms and computation
will be performed at construction time -- this may take considerable
memory and time! By default, this will produce the static structure factor.
To obtain the dynamic structure factor, set the `dyn_meas` keyword argument
to a non-zero value (the number of time-evolved snapshots to Fourier transform), along with proper settings for `dynΔt` and `meas_rate` (the
evolution timestep size, and how many timesteps are taken between snapshots).
See the documentation of [`StructureFactor`](@ref) for more details.

Alternatively, you can create a `StructureFactor` at the beginning of
your simulation by passing a single `SpinSystem`. The spin configuration
of this first system does not enter the averaged structure factor, as the
system is purely used to obtain information about the geometry.

Then, during the course of your simulation, call the `update!` function on
the `StructureFactor` along with your `SpinSystem` -- the current spin
configuration in the `SpinSystem` will be Fourier transformed and accumulated
into the structure factor. A bare-bones loop achieving this (assuming that you've already created a `system::SpinSystem` and a `sampler`) would look like:

```julia
sf = StructureFactor(system)
for _ in 1:therm_samples
   sample!(sampler)
   update!(sf, sys)
end
```

(In fact, the code for [`dynamic_structure_factor`](@ref) is not much more complex than this!) At the end of the loop, `sf` will hold the structure factor,
averaged across all of the thermal spin configurations it was provided.

## Calculating Structure Factors

This section details how the lower-level functions perform each step of computing the
structure factor.

To begin, we will assume you have on hand a large array storing a single spin trajectory. By
"single spin trajectory", we mean a single trajectory of a full system's worth of spins.
From an initial spin _configuration_, you can either just get the
static structure factor, or you will first need to perform Landau-Lifshitz dynamics using one of
the [Integrators](@ref) to construct a trajectory.

We will refer to this array as `spin_traj`, which can can be a `Array{SVector{3, Float64}}`
of size `size(spin_traj) == [B, D1, D2, D3, T]` (with `B` the number of basis sites,
`[D1, D2, D3]` the number of unit cells along each axis, and `T` the time axis).
Alternatively, this can be an `Array{ComplexF64}` of `size(spin_traj) == [3, B, D1, D2, D3, T]`
with the spins encoded into the real components. The former is more intuitive, but the
latter allows for in-place FFTs.

First, we need to perform a standard Fast Fourier Transform along the spatial and time axes.
This is be done with one of the following functions:

```_AT_docs
Sunny.fft_spin_traj
Sunny.fft_spin_traj!
```

As the documentation for the functions mentions, you will now have an array of `ComplexF64` of
size `[3, B, D1, D2, D3, T]`. (The spin component has been unfolded out into the first axis
regardless of the input format).

This could now be outer-producted with itself to form a contribution to the basis-resolved
structure factor. In particular, if `spin_traj_ft` is the name of your FFT'd spin trajectory,
the following function will perform this for you:

```julia
outerprod_conj(spin_traj_ft, (1, 2))
```

Which should result in a `ComplexF64` array of size `[3, 3, B, B, D1, D2, D3, T]`.
The documentation for this function can be seen below:

```_AT_docs
Sunny.outerprod_conj
Sunny.outerprod_conj!
```

Alternatively, if you only care about the post-basis-summation structure factor, you would
first want to instead perform the phase-weighted basis sum. This can be done manually, or
by using one of the following functions:

```_AT_docs
Sunny.phase_weight_basis
Sunny.phase_weight_basis!
```

As documented, this will return an array of `ComplexF64` of size `[3, Q1, ..., Qd, T]`,
where `Q1, ..., Qd` are the possibly expanded range of ``𝐪`` space requested through
`bz_size`.

As before, we can outer product this resulting array with itself to get a contribution to the
structure factor, now only required in the first axis as:

```julia
Sunny.outerprod_conj(spin_traj_ft, 1)
```
which should result in a `ComplexF64` array of size `[3, 3, D1, D2, D3, T]`

Repeat this entire process for all thermal spin trajectories you have at a given temperature,
average the result across all of them, and you have a dynamic structure factor! Note that if
you performed this entire process with an array containing a single spin configuration but
an extra "dummy" axis of length 1 (i.e. a size `[B, D1, D2, D3, 1]`), you would be left with
the _static_ structure factor!

There are additional functions which perform these accumulations while simulataneously
applying the neutron dipole form factor to reduce the spin components to a single
observable scalar. These are a bit of a mess currently, though.

```_AT_docs
Sunny.accum_dipole_factor!
Sunny.accum_dipole_factor_wbasis!
```
