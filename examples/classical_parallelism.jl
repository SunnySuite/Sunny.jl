# # Parallelizing Classical Structure Factor Calculations

using Sunny, GLMakie

# Calculating structure factors using classical dynamics is computationally
# expensive, and Sunny does not currently parallelize these calculations at a
# fine-grained level. However, Julia provides facilities that allow users to run
# multiple simulations in parallel with only a little extra effort. We will look
# at two approaches to doing this: [multithreading](https://docs.julialang.org/en/v1/manual/multi-threading/)
# and Julia's [`Distributed`](https://docs.julialang.org/en/v1/manual/distributed-computing/) package. 


# ## Review of the serial workflow

# The serial approach to calculating a structure factor, covered in [Structure
# Factors with Classsical Dynamics](@ref), involves generating equilibrium
# samples one at a time and then calling [`add_sample!`](@ref) on each of these.
# `add_sample!` uses the state of the `System` as an initial condition for the
# calculation of a dynamical trajectory. The correlations of the trajectory are
# calculated and accumulated into a running average.

# To illustrate, we'll set up a a simple model: a spin-1 antiferromagnet on an
# FCC crystal. We'll write a function to generate the model to make things
# easier later on.

function make_system(cryst; J, dims, seed=nothing)
    sys = System(cryst, dims, [SpinInfo(1, S=1, g=2)], :dipole; seed)
    set_exchange!(sys, J, Bond(1,1,[1,0,0]))
    return sys
end

cryst = Sunny.fcc_primitive_crystal()
sys = make_system(cryst; dims=(10,10,2), J=1.0)

# We call [`dynamical_correlations`](@ref) to create a `SampledCorrelations`.

sc = dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0)

# To thermalize and generate equilibrium samples, we'll need a
# [`Langevin`](@ref) integrator. We'll set the temperature to $k_B T = 0.5$ meV.

kT = 0.5
Δt = 0.05
integrator = Langevin(Δt; kT, λ=0.1);

# The serial calculation can then be performed as follows:

nsamples = 10

## Thermalize the system
for _ in 1:5000  # Sufficient number of steps to thermalize
    step!(sys, integrator)
end

for _ in 1:nsamples
    ## Generate new sample by running Langevin dynamics
    for _ in 1:1000  # Sufficient number of steps to decorrelate
        step!(sys, integrator)
    end
    ## Add the sample to the correlation data
    add_sample!(sc, sys)
end

# This takes a second or two on a modern workstation. Note that our
# `SampledCorrelations` now contains 10 samples. 

sc

#nb # Let's take a quick look at a path through reciprocal space. 

#nb qs = [[1,1,0], [1,1/2,0], [1,0,0]]
#nb path, xticks = reciprocal_space_path(cryst, qs, 40)
#nb formula = intensity_formula(sc, :perp; kT) 
#nb is = intensities_interpolated(sc, path, formula) 
#nb ωs = available_energies(sc) 
#nb heatmap(1:size(is, 1), ωs, is; colorrange=(0.0, 0.15), axis=(ylabel="Energy (meV)", xticks))

#nb # Note that the data is visibly noisy due to the small number of samples.

# ## Multithreading approach
# To use threads in Julia, you must launch your Julia environment appropriately.
# From the command line, this can be achieved with `julia --threads=N`, where
# `N` is the number of threads you'd like to use. If you don't know how many
# threads you'd like, you can let Julia determine an appropriate number with
# `julia --threads=auto`. If you are working in a Jupyter notebook, you will
# need to to set up a multithreaded Julia kernel. To do this, open a Julia REPL
# and type the following.
# ```
# using IJulia
# IJulia.installkernel("Julia Multithreaded", env=Dict(
#            "JULIA_NUM_THREADS" => "auto"
# ))
# ```
# After doing this, you'll have to restart Jupyter and make sure to select this
# kernel when launching a notebook.

# We will use multithreading in a very simple way, essentially employing a
# distributed memory approach to avoid conflicts around memory access. First
# we will preallocate a number of systems and correlations.  

npar = Threads.nthreads()
systems = [make_system(cryst; J=1.0, dims=(10,10,2), seed=i) for i in 1:npar]
scs = [dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0) for _ in 1:npar];


# !!! warning "Dealing with memory constraints"
#     
#   If you have many threads available and are working with a large system, you
#   may not have enough memory to store all these systems and correlations. In
#   that case, simply reduce `npar` to a small enough value that you can make the
#   necessary allocations.

# When the `Threads.@threads` macro is applied before a `for` loop, the
# iterations of the loop will execute in parallel using the available threads.
# We will put the entire thermalization and sampling process inside the loop,
# with each thread acting on a unique `System` and `SampledCorrelations`.

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

# You may find this takes a little bit longer than the serial example, but, at
# the end of it, you will have generated `npar` times more samples. Currently
# these samples are distributed across the separate `SampledCorrelations` in the
# `scs` array. We'll merge them into one summary `SampledCorrelations`. 

sc = merge_correlations(scs)

#nb # Finally we'll plot the results to compare with what we found before, finding
#nb # significantly improved results.
#nb 
#nb formula = intensity_formula(sc, :perp; kT)
#nb is = intensities_interpolated(sc, path, formula)
#nb ωs = available_energies(sc)
#nb heatmap(1:size(is, 1), ωs, is; colorrange=(0.0, 0.15), axis=(ylabel="Energy (meV)", xticks))

# ## Using `Distributed`
# Julia also provides a distributed memory approach to parallelism through the
# standard library package `Distributed` `Distributed` works by launching many
# different independent Julia environments. An advantage of this approach is
# that it scales naturally to clusters. A disadvantage, especially when working
# on a single computer, is the increased memory overhead associated with
# launching all these Julia environments. 
#
# We begin by importing the package.

using Distributed

# We now need to launch the Julia processes. The number of cores you have on
# your computer is often a good choice for the number of processes.

ncores = length(Sys.cpu_info())
addprocs(ncores);

# You can think of each of these processes as a separate computer running Julia.
# We now need to import Sunny into each Julia environment. This can be achieved
# with the `@everywhere` macro. 

@everywhere using Sunny

# We also need to define our functions on each of these processes. This again
# can be achieved with `@everywhere`.

@everywhere function make_system(cryst; J, dims, seed=nothing)
    sys = System(cryst, dims, [SpinInfo(1, S=1, g=2)], :dipole; seed)
    set_exchange!(sys, J, Bond(1,1,[1,0,0]))
    return sys
end;

# A simple way to perform work on each of these processes is to use the parallel
# map function, `pmap`. This will apply a function to each element of some
# iterable, such as a list of numbers, and return a list of the results. It is a
# _parallel_ map because these function calls may occur at the same time on
# different Julia processes. The `pmap` function takes care of distributing the
# work among the different processes and retrieving the results.

# In the example below, we give `pmap` a list of RNG seeds to iterate over, and
# we define the function that will be applied to each of these seeds in a `do`
# block. The contents of this block are essentially the same as what we put
# inside our parallel `for` loop in the multithreading example. The main
# difference is that the `System`s and `SampledCorrelations` are not created in
# advance of the parallelization but are instead created inside each Julia
# process. The `do` block returns a `SampledCorrelations`, and the output of all
# the parallel computations are collected into list of `SampledCorrelations`
# called `scs`.

@time scs = pmap(1:ncores) do seed
    ## Make our system and correlations
    sys = make_system(Sunny.fcc_primitive_crystal(); J=1.0, dims=(10,10,2), seed)
    sc = dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0)
    integrator = Langevin(Δt; kT, λ=0.1)

    ## Thermalize the system
    for _ in 1:5000
        step!(sys, integrator)
    end

    ## Collect samples
    for _ in 1:nsamples
        for _ in 1:1000
            step!(sys, integrator)
        end
        add_sample!(sc, sys) 
    end

    return sc
end;

# Finally, we will merge the results into a summary `SampledCorrelations`.

sc = merge_correlations(scs)

#nb # Because we've used the same seeds, the results look identical to those
#nb # achieved with multithreading.
#nb 
#nb formula = intensity_formula(sc, :perp; kT)
#nb is = intensities_interpolated(sc, path, formula)
#nb ωs = available_energies(sc)
#nb heatmap(1:size(is, 1), ωs, is; colorrange=(0.0, 0.15), axis=(ylabel="Energy (meV)", xticks))