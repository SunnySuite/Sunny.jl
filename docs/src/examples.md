# Examples

The examples stepped through here are available in `examples/` as full loadable files containing executable functions. Specifically, we work through here the simulations
performed in `examples/reproduce_testcases.jl`. All of the plotting code in here
currently depends on an installation of [Plots.jl](http://docs.juliaplots.org/latest/).

The high-level outline of performing a simulation is:

1. Create a [`Crystal`](@ref), either by providing explicit geometry information
    (Example 1), or by loading a `.cif` file (Example 2).
2. Using the `Crystal`, construct a [`Lattice`](@ref) which specifies the system
    size, and a collection of [Interactions](@ref) assembled into a [`Hamiltonian`](@ref).
3. Assemble a [`SpinSystem`](@ref) using the newly created `Lattice` and `Interaction`s.
4. Construct a sampler, either a [`LangevinSampler`](@ref) (Example 1), or a 
    [`MetropolisSampler`](@ref) (Example 2).
5. Use the sampler directly to sample new states, or use it to perform [Structure factor calculations](@ref).

Defining interactions in step (2) can be aided by our utilities for symmetry analysis, demonstrated at the bottom of this page.

In all examples, we will assume that `FastDipole` and `StaticArrays` have been loaded.

## Example 1: Diamond lattice with antiferromagnetic Heisenberg interactions

In this example, we will step through the basic steps needed to set up and run
a spin dynamics simulation, with finite-$T$ statistics obtained by using Langevin
dynamics. The full example is contained in the function `test_diamond_heisenberg_sf()`
within `examples/reproduce_testcases.jl`.

**(1)** We construct a diamond lattice by explicitly defining the lattice geometry. We will use the conventional cubic Bravais lattice with an 8-site basis. Our simulation box
will be ``8 \times 8 \times 8`` unit cells along each axis. Since we're not thinking about a specific system, we label all of the sites with the arbitrary species label `"A"`.

```julia
lat_vecs = SA[4.0 0.0 0.0;
              0.0 4.0 0.0;
              0.0 0.0 4.0]
basis_vecs = [
    SA[0.0, 0.0, 0.0],
    SA[0.0, 2.0, 2.0],
    SA[2.0, 0.0, 2.0],
    SA[2.0, 2.0, 0.0],
    SA[3.0, 3.0, 3.0],
    SA[3.0, 1.0, 1.0],
    SA[1.0, 3.0, 1.0],
    SA[1.0, 1.0, 3.0]
]
basis_labels = fill("A", 8)
latsize = SA[8, 8, 8]
lattice = Lattice{3}(lat_vecs, basis_vecs, basis_labels, latsize)
crystal = Crystal(lattice)
```

**(2)** In step 1, we ended up already creating our `Lattice`, so all that is left is
to define our Heisenberg interactions. We want to set up nearest-neighbor antiferromagnetic interactions with a strength of ``J = 28.28~\mathrm{K}``. One nearest-neighbor bond is the one connecting basis site 3 with basis site 6 within a single unit cell. (We can figure this out using our tools for symmetry analysis, at the bottom of this page).

```julia
J = 28.28           # Units of K
interactions = [
    Heisenberg(J, crystal, Bond{3}(3, 6, SA[0,0,0])),
]
ℋ = Hamiltonian{3}(interactions)
```

**(3)** Assembling a `SpinSystem` is straightforward. Then, we will randomize the system so that all spins are randomly placed on the unit sphere.

```julia
sys = SpinSystem(lattice, ℋ)
rand!(sys)
```

**(4)** We will simulate this system using Langevin dynamics, so we need to create a [`LangevinSampler`](@ref). Note that the units of integration time and temperature are relative to the units implicitly used when setting up the interactions.

```julia
Δt = 0.02 / J       # Units of 1/K
kT = 4.             # Units of K
α  = 0.1
kB = 8.61733e-5     # Units of eV/K
nsteps = 20000
sampler = LangevinSampler(sys, kT, α, Δt, nsteps)
```

At this point we can call `sample!(sampler)` to produce new samples of the system, which will be reflected in the state of `sys`. Instead, we will proceed to calculate
the finite-$T$ structure factor using our built-in routines.

**(5)** The full process of calculating a structure factor is handled
by [`structure_factor`](@ref). Internally, this function:

1. Thermalizes the system for a while
2. Samples a new thermal spin configuration
3. Performs constant-energy LL dynamics to obtain a Fourier-transformed
    dynamics trajectory. Use this trajectory to calculate a structure
    factor contribution ``S^{\alpha,\beta}(\boldsymbol{q}, \omega)``.
4. Repeat steps (2,3), averaging structure factors across samples.

See the documentation of [`structure_factor`](@ref) for details of how
this process is controlled by the function arguments, and how to properly
index into the resulting array.

In this example, we will just look at the diagonal elements of this
matrix along some cuts in reciprocal space. To improve statistics,
we average these elements across the ``x, y, z`` spin directions since
they are all symmetry equivalent in this model.

```julia
meas_rate = 10
S = structure_factor(
    sys, sampler; num_samples=5, dynΔt=Δt, meas_rate=meas_rate,
    num_freqs=1600, bz_size=(1,1,2), therm_samples=10, verbose=true
)

# Retain just the diagonal elements, which we will average across the
#  symmetry-equivalent directions.
avgS = zeros(Float64, axes(S)[3:end])
for α in 1:3
    @. avgS += real(S[α, α, :, :, :, :])
end

```

We then plot some cuts using a function `plot_many_cuts` defined within
the example script. (I.e. this code block will not successfully execute unless
you `include("examples/reproduce_testcases.jl)`). We omit this code here as it's just
a large amount of indexing and plotting code, but for details see the script.

```julia
# Calculate the maximum ω present in our FFT
# Need to scale by (S+1) with S=3/2 to match the reference,
#  and then convert to meV.
maxω = 1000 * 2π / ((meas_rate * Δt) / kB) / (5/2)
p = plot_many_cuts(avgS; maxω=maxω, chopω=5.0)
display(p)
```


## Example 2: FeI₂ with a complex collection of interactions

In this example, we work through performing a more complicated and realistic
simulation. While the number of interactions is much larger, the process is
no more complicated. We will also see how to perform sampling using Metropolis
Monte Carlo through the [`MetropolisSampler`](@ref) type. The full example is
contained in the function `test_FeI2_MC()` within
`examples/reproduce_testcases.jl`.

(To be filled out).

## Symmetry analysis

When defining pair interactions, we are always defining the interactions on
entire symmetry classes at the same time. To do this, we need to provide the
exchange matrix ``J`` on a specific reference `Bond`, which is then automatically
propagated to all symmetry-equivalent bonds. However, on any given bond, the
exchange matrix must live within a restricted space of ``3 \times 3`` matrices
that is confined by the symmetry properties of the underlying crystal.

To discover all symmetry classes of bonds up to a certain distance while simultaneously learning what the allowed form of the `J` matrix is, construct a `Crystal` then call the function [`print_bond_table`](@ref).

```
julia> lattice = FastDipole.diamond_conventional(1.0, (8, 8, 8))
julia> crystal = Crystal(lattice)
julia> print_bond_table(crystal, 4.0)

Bond{3}(1, 1, [0, 0, 0])
Distance 0, multiplicity 1
Connects [0, 0, 0] to [0, 0, 0]
Allowed coupling:  |A 0 0 |
                   |0 A 0 |
                   |0 0 A |

Bond{3}(3, 6, [0, 0, 0])
Distance 1.732, multiplicity 4
Connects [0.5, 0, 0.5] to [0.75, 0.25, 0.25]
Allowed coupling:  | A  B -B |
                   | B  A -B |
                   |-B -B  A |

Bond{3}(1, 2, [0, 0, 0])
Distance 2.828, multiplicity 12
Connects [0, 0, 0] to [0, 0.5, 0.5]
Allowed coupling:  | B  D  D |
                   |-D  C  A |
                   |-D  A  C |

Bond{3}(1, 6, [0, 0, 0])
Distance 3.317, multiplicity 12
Connects [0, 0, 0] to [0.75, 0.25, 0.25]
Allowed coupling:  |B C C |
                   |C D A |
                   |C A D |

Bond{3}(1, 1, [1, 0, 0])
Distance 4, multiplicity 6
Connects [0, 0, 0] to [1, 0, 0]
Allowed coupling:  |A 0 0 |
                   |0 B 0 |
                   |0 0 B |
```

Each block represents one symmetry equivalence class of bonds, along with a single
representative ("canonical") `Bond` for that class and the allowed exchange coupling
matrix on that canonical bond.

You can also query what the allowed exchange matrix is on a specific bond using [`allowed_J`](@ref).

```
julia> allowed_J(crystal, Bond{3}(1, 5, [1,-1,0]))

3×3 Matrix{String}:
 "D"  "A"  "B"
 "A"  "E"  "C"
 "B"  "C"  "F"
```