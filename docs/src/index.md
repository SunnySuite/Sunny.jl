# Sunny.jl

[Sunny](https://github.com/SunnySuite/Sunny.jl/) is a package for simulating
classical spin systems, including the Landau-Lifshitz dynamics of spin dipoles
and its generalization to multipolar spin components. The latter is especially
powerful for modeling magnetic compounds with strong single-ion anisotropy
interactions.

Sunny provides the following features:

- Generalized spin dynamics using [SU(_N_) coherent states](https://arxiv.org/abs/2209.01265).
- Ability specify a crystal by a `.cif` file, or using its spacegroup symmetry.
- Symmetry analysis to classify allowed interaction terms, and to propagate them by symmetry.
- Single-ion anisotropy at arbitrary order, which can be specified using Stevens operators or as a polynomial of spin operators.
- Monte Carlo sampling of spin configurations in thermal equilibrium.
- Ewald summation for long-range dipole-dipole interactions, accelerated with the fast Fourier transform (FFT).
- Estimation of the $\mathcal{S}(\mathbf{q}, \omega)$ dynamical structure factor data, with options for various corrections (form factor, classical-to-quantum factors, ...)

Work in progress includes:

- Linear spin wave theory and its generalization to SU(_N_) coherent states.
- Interactive visualizations of the 3D crystals and structure factor data.
- MPI-distributed Monte Carlo sampling, including [parallel tempering](https://en.wikipedia.org/wiki/Parallel_tempering).
