# Parallelizing Classical Structure Factor Calculations

Calculating structure factors using classical dynamics is computationally
expensive, and Sunny does not currently parallelize these calculations at a
fine-grained level. However, Julia provides facilities that allow users to run
multiple simulations in parallel with only a little extra effort. We will look
at two approaches to doing this:
[multithreading](https://docs.julialang.org/en/v1/manual/multi-threading/) and
Julia's
[`Distributed`](https://docs.julialang.org/en/v1/manual/distributed-computing/)
package.
We'll present these approaches in a series of code snippets that can be
copied and pasted into your preferred Julia development environment.

## Review of the serial workflow

The serial approach to calculating a structure factor, covered in [Structure
Factors with Classsical Dynamics](@ref), involves thermalizing a spin `System`
and then calling [`add_sample!`](@ref). `add_sample!` uses the state of the
`System` as an initial condition for the calculation of a dynamical
trajectory. The correlations of the trajectory are calculated and accumulated
into a running average of the ``𝒮(𝐪,ω)``. This processes is repeated to
generate additional samples.

To illustrate, we'll set up a a simple model: a spin-1 antiferromagnet on an
FCC crystal. We'll write a function to generate the model to make things
easier later on.

```julia
using Sunny, GLMakie

function make_system(cryst; J, dims, seed=nothing)
    sys = System(cryst, dims, [SpinInfo(1, S=1, g=2)], :dipole; seed)
    set_exchange!(sys, J, Bond(1,1,[1,0,0]))
    return sys
end

cryst = Sunny.fcc_primitive_crystal()
sys = make_system(cryst; dims=(10,10,2), J=1.0)
```

We call [`dynamical_correlations`](@ref) to create a `SampledCorrelations`.

```julia
sc = dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0)
```

To thermalize and generate equilibrium samples, we'll need a
[`Langevin`](@ref) integrator. We'll set the temperature to $k_B T = 0.5$ meV.

```julia
kT = 0.5
Δt = 0.05
integrator = Langevin(Δt; kT, λ=0.1)
```

The serial calculation can then be performed as follows:

```julia
nsamples = 10

# Thermalize the system
for _ in 1:5000  # Sufficient number of steps to thermalize
    step!(sys, integrator)
end

for _ in 1:nsamples
    # Generate new sample by running Langevin dynamics
    for _ in 1:1000  # Sufficient number of steps to decorrelate
        step!(sys, integrator)
    end
    # Add the sample to the correlation data
    add_sample!(sc, sys)
end
```

This will take a second or two on a modern workstation, resulting in a single
`SampledCorrelations` that contains 10 samples.


## Multithreading approach
To use threads in Julia, you must launch your Julia environment appropriately.
From the command line, this can be achieved with `julia --threads=N`, where
`N` is the number of threads you'd like to use. If you don't know how many
threads you'd like, you can let Julia determine an appropriate number with
`julia --threads=auto`. If you are working in a Jupyter notebook, you will
need to to set up a multithreaded Julia kernel. To do this, open a Julia REPL
and type the following.
```
using IJulia
IJulia.installkernel("Julia Multithreaded", env=Dict("JULIA_NUM_THREADS" => "auto"))
```
After doing this, you'll have to restart Jupyter and make sure to select this
kernel when launching a notebook.

We will use multithreading in a very simple way, essentially employing a
distributed memory approach to avoid conflicts around memory access. First
we will preallocate a number of systems and correlations.

```julia
npar = Threads.nthreads()
systems = [make_system(cryst; J=1.0, dims=(10,10,2), seed=i) for i in 1:npar]
scs = [dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0) for _ in 1:npar]
```

!!! warning "Dealing with memory constraints"

    If you have many threads available and are working with a large system, you
    may not have enough memory to store all these systems and correlations. In
    that case, simply reduce `npar` to a small enough value that you can make the
    necessary allocations.

When the `Threads.@threads` macro is applied before a `for` loop, the
iterations of the loop will execute in parallel using the available threads.
We will put the entire thermalization and sampling process inside the loop,
with each thread acting on a unique `System` and `SampledCorrelations`.

```julia
Threads.@threads for id in 1:npar
    integrator = Langevin(Δt; kT, λ=0.1)
    for _ in 1:5000
        step!(systems[id], integrator)
    end
    for _ in 1:nsamples
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
samples with `merge_correlations`:

```julia
sc = merge_correlations(scs)
```

## Using `Distributed`
Julia also provides a distributed memory approach to parallelism through the
standard library package `Distributed`. This works by launching
independent Julia environments on different "processes". An advantage of this
approach is that it scales naturally to clusters -- the processes are easily
distributed across many different compute nodes. A disadvantage, especially when
working on a single computer, is the increased memory overhead associated with
launching many Julia environments.

We begin by importing the package.

```julia
using Distributed
```

We now need to launch the Julia processes. The number of cores you have on
your computer is often a good choice for the number of processes.

```julia
ncores = length(Sys.cpu_info())
addprocs(ncores)
```

You can think of each process as a separate computer running Julia. We now need
to import Sunny into all of these Julia environments. This can be achieved with
the `@everywhere` macro.

```julia
@everywhere using Sunny
```

We also need to define our functions on all of these processes. This again
can be achieved with `@everywhere`.

```julia
@everywhere function make_system(cryst; J, dims, seed=nothing)
    sys = System(cryst, dims, [SpinInfo(1, S=1, g=2)], :dipole; seed)
    set_exchange!(sys, J, Bond(1,1,[1,0,0]))
    return sys
end
```

A simple way to perform work on each of these processes is to use the parallel
map function, `pmap`. This will apply a function to each element of some
iterable, such as a list of numbers, and return a list of the results. It is a
_parallel_ map because these function calls may occur at the same time on
different Julia processes. The `pmap` function takes care of distributing the
work among the different processes and retrieving the results.

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
@time scs = pmap(1:ncores) do seed
    sys = make_system(Sunny.fcc_primitive_crystal(); J=1.0, dims=(10,10,2), seed)
    sc = dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0)
    integrator = Langevin(Δt; kT, λ=0.1)

    for _ in 1:5000      # Thermalize
        step!(sys, integrator)
    end
    for _ in 1:nsamples 
        for _ in 1:1000 
            step!(sys, integrator)
        end
        add_sample!(sc, sys)
    end

    return sc
end
```

Finally, we will merge the results into a summary `SampledCorrelations`.

```julia
sc = merge_correlations(scs)
```