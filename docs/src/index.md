# Sunny.jl

This is the documentation for [Sunny.jl](https://github.com/MagSims/Sunny.jl). This package aims to be a general-purpose library for simulating and analyzing classical spin systems obtained as the ``\mathrm{SU(N)}`` classical limit of quantum spin Hamiltonians.

**Current Features**
- ``\mathrm{SU(2)}`` Landau-Lifshitz spin dynamics
- Define arbitrary spin Hamiltonians involving supported interactions:
    - External fields
    - Arbitrary generalized pair interactions
    - On-site anisotropies
    - Long-range dipole-dipole interactions
- Accelerated dipole-dipole interactions using Ewald summation and Fourier-space evaluation
- Finite-``T`` sampling using either Langevin dynamics or Metropolis Monte Carlo
- Structure factor calculations
- Automated symmetry analysis and bond equivalency class discovery
- Interactive plotting using [Makie.jl](https://github.com/JuliaPlots/Makie.jl)

**Planned Features**
- CPU parallelization of dynamics and Monte Carlo simulations
- Parallel tempering
- Generalized ``\mathrm{SU(N)}`` models
- GPU acceleration
- Accelerated local MC updates in the presence of long-range dipole interactions

To start running your first simulations, head over to [Getting Started](@ref) to install the package. Then, see [Examples](@ref) for a set of examples demonstrating various features.

Or, take a look at the full [Library](@ref) documentation.