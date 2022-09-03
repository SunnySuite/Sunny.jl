# Examples

Please take a look at our Jupyter notebook [tutorials].

Additional working examples are available in `examples/` as fully loadable files containing
executable functions. We will soon update this page with additional information
that will walk you through these examples. 

The high-level outline of performing a simulation is:

1. Create a [`Crystal`](@ref), either by providing explicit geometry information
    (Example 1), or by loading a `.cif` file (Example 2).
2. Using the `Crystal`, construct a collection of [Interactions](@ref).
3. Assemble a [`SpinSystem`](@ref) using the newly created `Crystal` and interactions, the size of the simulation box, and optionally information on the spin magnitude and ``g``-tensors of each site by passing a list of `SiteInfo`.
4. Construct a sampler, either a [`LangevinSampler`](@ref) (Example 1), or a 
    [`MetropolisSampler`](@ref) (Example 2).
5. Use the sampler directly to sample new states, or use it to perform [Structure factor calculations](@ref).

Sunny provides a full suite of symmetry analysis tools to ensure that the interactions
specified in step 2 are valid for the specified crystal.