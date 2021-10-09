# Internals

This page documents various components of how Sunny works internally, and how one
might be able to extend it with new functionalities. This page is very under-complete.
Documenting internals as I go.

## Handling `Interaction`s

Interactions exist at two levels:

1. The types that the user create and interface with (living in `Interactions.jl`)
2. The types that this gets converted behind the scenes upon creating a `SpinSystem`
    (which are scattered throughout the codebase).

To define a new type of `Interaction` which can appear in a [`Hamiltonian`](@ref),
one must provide both of these types (which may be the same!), a way to convert from
the first to the second, and a collection of functionalities on the second
needed for various simulation tasks. We provide the list of required steps below, and for
explicit examples see the core interactions defined across `src/Interactions.jl`
and `src/PairInteractions.jl`.

**(1)** Define a new struct which is a subtype of `Interaction`. This is the user-facing
         type, which should be the minimal specification needed to specify the
         interaction. As much as possible, instances of this type should be agnostic to
         the final crystal geometry they'll be placed on. We'll refer to this type
         here as `MyInt`.

**(2)** Create _another_ struct which will handle how to actually explicitly compute
         terms arising from your interaction on a specific lattice. We'll refer to this
         type here as `MyIntInternal`. This should expose a constructor
         `MyIntInternal(int::MyInt, crystal::Crystal, latsize)`.

**(3)** Provide the following methods:

- `Sunny.energy(sys::SpinSystem, int::MyIntInternal)` which computes the total term in the
     Hamiltonian given the state of the system.
- `Sunny._accum_field!(B::Array{Vec3}, spins::Array{Vec3}, int::MyIntInternal)` which
    accumulates the local "field" coming from this interaction into `B` given the state
    of the `spins`. Specifically, this function call should perform:
        ``𝐁ᵢ = 𝐁ᵢ - ∇_{𝐒ᵢ} ℋ_I``
    where ``ℋ_I`` is your new interaction term.


**(4)** Edit the `HamiltonianCPU` struct (in `Hamiltonian.jl`) to store a `Vector{MyIntInternal}`,
    and to call your defined `energy` and `_accum_field!` functions within
    the existing `energy(spins, ℋ::HamiltonianCPU)` and `field!(B, spins, ℋ)` functions.

**(5)** (Optional, for Metropolis sampling support) In `Metropolis.jl`,
edit `local_energy_change(sys, idx, newspin)` to compute the change in energy from your
interaction resulting from changing the spin of `sys` at `idx` to `newspin`.

The current design, while slightly unwieldy, is as it is to (1) separate Hamiltonian *definition*
from Hamiltonian *implementation* and (2) to avoid paying dispatch costs at
runtime. For example, it would be extremely clean and simple if `HamiltonianCPU` simply stored
a `Vector{<:IntInternal}`, and looped over this to call `energy` and `_accum_field!` functions.
However, we would have to pay at runtime to constantly look up in dispatch tables which version
of `energy` and `_accum_field!` we're calling at each loop iteration.

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

```@docs
fft_spin_traj
fft_spin_traj!
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

```@docs
outerprod_conj
outerprod_conj!
```

Alternatively, if you only care about the post-basis-summation structure factor, you would
first want to instead perform the phase-weighted basis sum. This can be done manually, or
by using one of the following functions:

```@docs
phase_weight_basis
phase_weight_basis!
```

As documented, this will return an array of `ComplexF64` of size `[3, Q1, ..., Qd, T]`,
where `Q1, ..., Qd` are the possibly expanded range of ``𝐪`` space requested through
`bz_size`.

As before, we can outer product this resulting array with itself to get a contribution to the
structure factor, now only required in the first axis as:

```julia
outerprod_conj(spin_traj_ft, 1)
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

```@docs
accum_dipole_factor!
accum_dipole_factor_wbasis!
```
