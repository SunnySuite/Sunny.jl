# Parallelizing Calculations

Calculating structure factors with classical dynamics is computationally
expensive, and Sunny does not currently parallelize these calculations at a
fine-grained level. However, Julia provides facilities that allow users to run
multiple simulations in parallel with only a little extra effort. We will look
at two approaches to doing this:
[multithreading](https://docs.julialang.org/en/v1/manual/multi-threading/) and
Julia's
[`Distributed`](https://docs.julialang.org/en/v1/manual/distributed-computing/)
package. We'll present these approaches in a series of code snippets that can be
copied and pasted into your preferred Julia development environment.

## Review of the serial workflow

The serial approach to calculating a structure factor, covered in the [FeI₂
tutorial](@ref "4. Generalized spin dynamics of FeI₂ at finite *T*"), involves
thermalizing a spin [`System`](@ref) and then calling [`add_sample!`](@ref).
`add_sample!` uses the state of the `System` as an initial condition for the
calculation of a dynamical trajectory. The correlations of the trajectory are
calculated and accumulated into a running average of the ``\mathcal{S}(𝐪,ω)``.
This sequence is repeated to generate additional samples.

To illustrate, we'll set up a a simple model: a spin-1 antiferromagnet on a BCC
crystal. Constructing the `System` with a specific random number `seed` ensures
full reproducibility of the simulation.

```julia
using Sunny, GLMakie

function make_system(seed)
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    positions = [[0, 0, 0]/2, [1, 1, 1]/2]
    cryst = Crystal(latvecs, positions)
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(10, 10, 2), seed)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    return sys
end

sys = make_system(0)
```

A serial calculation of [`SampledCorrelations`](@ref) involving the
[`Langevin`](@ref) sampling method can now be performed as follows:

```julia
# Thermalize the system
dt = 0.05
integrator = Langevin(dt; damping=0.2, kT=0.5)
for _ in 1:5000
    step!(sys, integrator)
end

# Accumulator for S(q,ω) samples

energies = range(0.0, 10.0, 100)
measure = ssf_perp(sys)
sc = SampledCorrelations(sys; dt=0.1, energies, measure)

# Collect 10 samples
for _ in 1:10
    for _ in 1:1000
        step!(sys, integrator)
    end
    add_sample!(sc, sys)
end
```

This will take a second or two on a modern workstation, resulting in a single
`SampledCorrelations` that contains 10 samples.


## Multithreading approach
To use threads in Julia, you must launch your Julia environment appropriately.

- From the **command line**, launch Julia with `julia --threads=auto`. With this option, Julia will automatically use an optimal number of threads.
- **Jupyter notebook** users will need to to set up a multithreaded Julia kernel and restart into this kernel. The kernel can be created inside Julia with the command `IJulia.installkernel("Julia Multithreaded", env=Dict("JULIA_NUM_THREADS" => "auto"))`.
- **VSCode** users should open their settings and search for `Julia: Additional Args`. There will be link called `Edit in settings.json`. Click on this, add `"--threads=auto"` to the list `julia.additionalArgs`, and start a new REPL.

Before going further, make sure that `Threads.nthreads()` returns a number greater than 1.

We will use multithreading in a very simple way, essentially employing a
distributed memory approach to avoid conflicts around memory access. First
preallocate a number of systems and correlations. The integer `id` of each
system is used as its random number seed.

```julia
npar = Threads.nthreads()
systems = [make_system(id) for id in 1:npar]
scs = [SampledCorrelations(sys; dt=0.1, energies, measure) for _ in 1:npar]
```

!!! warning "Dealing with memory constraints"  
    If many threads are available and the system sizes are large, there may not
    be enough computer memory available to store the dynamical correlation data.
    In this case, consider reducing the number of parallel threads, `npar`.

When the `Threads.@threads` macro is applied before a `for` loop, the
iterations of the loop will execute in parallel using the available threads.
We will put the entire thermalization and sampling process inside the loop,
with each thread acting on a unique `System` and `SampledCorrelations`.

```julia
Threads.@threads for id in 1:npar
    integrator = Langevin(dt; damping=0.2, kT=0.5)
    for _ in 1:5000
        step!(systems[id], integrator)
    end
    for _ in 1:10
        for _ in 1:1000
            step!(systems[id], integrator)
        end
        add_sample!(scs[id], systems[id])
    end
end
```

You may find this takes a little bit longer than the serial example, but, at the
end of it, you will have generated `npar` correlation objects, each with 10
samples. We can merge these into a summary `SampledCorrelations` with `10*npar`
samples with [`merge_correlations`](@ref):

```julia
sc = merge_correlations(scs)
```

## Using `Distributed`
Julia also provides a distributed memory approach to parallelism through the
standard library package `Distributed`. This works by launching independent
Julia environments on different "processes." An advantage of this approach is
that it scales naturally to clusters since the processes are easily distributed
across many different compute nodes. A disadvantage, especially when working on
a single computer, is the increased memory overhead associated with launching
all these environments.

We begin by importing the package,

```julia
using Distributed
```

and launching some new processes. It is often sensible to create as many
processes as there are cores.

```julia
ncores = length(Sys.cpu_info())
addprocs(ncores)
```

You can think of each process as a separate computer running a fresh Julia
environment, so we'll need to import Sunny and define functions in each of these
environments. This is easily achieved with the `@everywhere` macro.
```julia
@everywhere using Sunny

@everywhere function make_system(seed)
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    positions = [[0, 0, 0]/2, [1, 1, 1]/2]
    cryst = Crystal(latvecs, positions)
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; seed)
    sys = resize_supercell(sys, (10, 10, 2))
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    return sys
end
```

A simple way to perform work on these processes is to use the parallel map
function, `pmap`. This will apply a function to each element of some iterable,
such as a list of numbers, and return a list of the results. It is a _parallel_
map because these function calls may occur at the same time on different Julia
processes. The `pmap` function takes care of distributing the work among the
different processes and retrieving the results.

In the example below, we give `pmap` a list of RNG seeds to iterate over, and
we define the function that will be applied to each of these seeds in a `do`
block. The contents of this block are essentially the same as what we put
inside our parallel `for` loop in the multithreading example. The main
difference is that the `System`s and `SampledCorrelations` are not created in
advance of the parallelization but are instead created inside each Julia
process. The `do` block returns a `SampledCorrelations`, and the output of all
the parallel computations are collected into list of `SampledCorrelations`
called `scs`.

```julia
scs = pmap(1:ncores) do id
    sys = make_system(id)
    sc = SampledCorrelations(sys; dt=0.1, energies, measure)
    integrator = Langevin(0.05; damping=0.2, kT=0.5)

    for _ in 1:5000
        step!(sys, integrator)
    end
    for _ in 1:10
        for _ in 1:1000 
            step!(sys, integrator)
        end
        add_sample!(sc, sys)
    end

    return sc
end
```

Finally, merge the results into a summary `SampledCorrelations`.

```julia
sc = merge_correlations(scs)
```
