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
