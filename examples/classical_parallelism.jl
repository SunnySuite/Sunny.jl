# # Parallelizing Classical Structure Factor Calculations

using Sunny, GLMakie

# Calculating structure factors using classical dynamics can be computationally
# expensive, and Sunny does not currently parallelize these calculations at a
# fine-grained level. However, Julia provides facilities that allow users to run
# multiple simulations in parallel with only a little extra effort. We will look
# at two approaches to doing this, using multithreading and using Julia's
# `Distributed` package. 


# ## Review of the serial workflow

# The serial approach to calculating a structure factor, covered in [Structure
# Factors with Classsical Dynamics](@ref), involves generating thermodynamic
# samples one at a time and then calling [`add_sample!`](@ref) on each of these.
# `add_sample!` uses the state of the `System` as an initial condition for the
# generation of a dynamical trajectory. The correlations of the trajectory are
# calculated and accumulated into a running average.

# To illustrate, we'll set up a a simple model: a spin-1 antiferromagnet on an FCC
# crystal. We'll write a function to produce this system to make things easier
# later on.

function make_afm_system(cryst; dims=(10,10,2), J=1.0, seed=nothing)
    sys = System(cryst, dims, [SpinInfo(1, S=1, g=2)], :dipole; seed)
    set_exchange!(sys, J, Bond(1,1,[1,0,0]))
    return sys
end

cryst = Sunny.fcc_primitive_crystal()
sys = make_afm_system(cryst)

# And we will call [`dynamical_correlations`](@ref) to create an object store
# correlation data.

sc = dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0)

# To thermalize and generaite equilibrium samples, we'll need a
# [`Langevin`](@ref) integrator. We'll set the temperature to $k_B T = 0.5$ meV.

kT = 0.5
Δt = 0.05
integrator = Langevin(Δt; kT, λ=0.1)

# The serial calculation can then be performed as follows:
#!nb # We'll time this process to get a rough benchmark.

nsamples = 10

#!nb @time begin
## Thermalize the system
for _ in 1:5000  # Sufficient number of steps to thermalize
    step!(sys, integrator)
end

@time for _ in 1:nsamples
    ## Generate new sample by running Langevin dynamics
    for _ in 1:1000  #Sufficient number of steps to decorrelate
        step!(sys, integrator)
    end
    ## Add the sample to the correlation data
    add_sample!(sc, sys)
end
#!nb end

# This takes about 1.5 seconds on a modern workstation. Let's take a quick look
# at a path through reciprocal space.

qs = [[1,1,0], [1,1/2,0], [1,0,0]]
path, xticks = reciprocal_space_path(cryst, qs, 40)
formula = intensity_formula(sc, :perp; kT)
is = intensities_interpolated(sc, path, formula)

ωs = available_energies(sc)
heatmap(1:size(is, 1), ωs, is; colorrange=(0.0, 0.15), axis=(ylabel="Energy (meV)", xticks))

# Note that the data is visibly noisy due to the small number of samples
# collected. 

# ## Multithreading approach
# To use [threads](https://docs.julialang.org/en/v1/manual/multi-threading/) in
# Julia, it is important to launch your Julia environment appropriately. From
# the command line, this can be achieved with `julia --threads=N`, where `N` is
# the number of threads you'd like to use. If you don't know how many threads
# you'd like, you can let Julia determine an appropriate number with `julia
# --threads=auto`. If you are working in a Jupyter notebook, you will need to
# to set up a multithreaded Julia kernel.
#nb # To do this, open a Julia REPL and type the following.
#nb # ```
#nb # using IJulia
#nb # IJulia.installkernel("Julia Multithreaded", env=Dict(
#nb #            "JULIA_NUM_THREADS" => "auto"
#nb # ))
#nb # ```
#nb # This will create a new kernel using as many threads as you have available.
#nb # After doing this, you'll have to restart Jupyter and make sure to select this kernel when
#nb # launching a notebook.

# We will use multithreading in a very simple way, essentially employing a
# distributed memory approach to avoid issues around data races. First we'll
# preallocate a number of systems and correlations.  

npar = Threads.nthreads()
systems = [make_afm_system(cryst; seed=i) for i in 1:npar]
scs = [dynamical_correlations(sys; Δt=0.1, nω=100, ωmax=10.0) for _ in 1:npar]

# If you have many threads available and are working with a large system, you
# may not have enough memory to store all these systems and correlations. In
# that case, simply reduce `npar` to a small enough value that you can make the
# necessary allocations.

# When the `Threads.@threads` macro is applied before a `for` loop, the
# iterations of the loop will execute in parallel using the available threads.
# We will put the entire thermalization and sampling process inside the loop,
# with each thread acting on a unique `System` and `SampledCorrelations`.

@time Threads.@threads for id in 1:npar
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

# You should find this takes a little bit longer than the serial example, but,
# at the end of it, you will have generated `npar` more samples. Currently these
# samples are split across each element of the `scs` array. We'll merge them
# into one summary `SampledCorrelations`. 

sc = Sunny.merge_correlations(scs)

# Finally we'll plot the results to compare with what we found before, finding
# significantly improved results.

formula = intensity_formula(sc, :perp; kT)
is = intensities_interpolated(sc, path, formula)
ωs = available_energies(sc)
heatmap(1:size(is, 1), ωs, is; colorrange=(0.0, 0.15), axis=(ylabel="Energy (meV)", xticks))

# ## Using `Distributed`
# Julia also provides a distributed memory approach to parallelism through the
# standard library package
# [`Distributed`](https://docs.julialang.org/en/v1/manual/distributed-computing/).
# This approach works by launching many different independent Julia
# environments. An advantage of this approach is that it scales naturally to
# clusters. A disadvantage, especially when working on a single computer, is the
# increased memory overhead associated with launching all these Julia
# environments. 
#
# We begin by importing the package.

using Distributed

# We now need to launch the Julia processes. The number of cores you have on
# your computer is often a good choice for the number of processes.

addprocs(32)

# This will launch 32 Julia environments. We now need to import Sunny into
# each of these environments. This can be achieved with the `@everywhere` macro. 

@everywhere using Sunny

# To make a function available in each of these Julia environments,
# we again need to use `@everywhere`.

@everywhere function make_afm_system(cryst; dims=(10,10,2), J=1.0, seed=nothing)
    sys = System(cryst, dims, [SpinInfo(1, S=1, g=2)], :dipole; seed)
    set_exchange!(sys, J, Bond(1,1,[1,0,0]))
    return sys
end

# A simple way to perform work on each of these process is to use the parallel
# map function, `pmap`. We'll execute on each process the same computation we
# executed above in an individual thread. Here, however, we'll do our
# allocations within each process instead of allocating in advance. 

npar = 32

@time scs = pmap(1:npar) do seed
    ## Make our system and correlations
    sys = make_afm_system(Sunny.fcc_primitive_crystal(); seed)
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
end

# `pmap` will collect the return results into a list, which we will merge
# into a summary `SampledCorrelations`.

sc = Sunny.merge_correlations(scs)

# Because we've used the same seeds, the results look identical to those
# achieved with multithreading.

formula = intensity_formula(sc, :perp; kT)
is = intensities_interpolated(sc, path, formula)
ωs = available_energies(sc)
heatmap(1:size(is, 1), ωs, is; colorrange=(0.0, 0.15), axis=(ylabel="Energy (meV)", xticks))