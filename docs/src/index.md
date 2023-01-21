# Sunny.jl

[Sunny](https://github.com/SunnySuite/Sunny.jl/) is a package for simulating
classical spin systems, including the Landau-Lifshitz dynamics of spin dipoles
and its generalization to multipolar spin components. In the latter case, Sunny
resolves the local quantum structure of individual spins, making it particularly
suited for modeling magnetic compounds with strong local anisotropy.

Sunny provides the following features:

- Generalization of Landau-Lifshitz spin dynamics using the [formalism of SU(_N_) coherent states](https://arxiv.org/abs/2209.01265).
- Ability specify a crystal by a `.cif` file, or using its spacegroup symmetry.
- Symmetry analysis to classify allowed Hamiltonian terms, and to automatically populate all symmetry equivalent interactions.
- Single-ion anisotropy at arbitrary order, which can be specified using Stevens operators or as an arbitrary polynomial of spin operators.
- Sampling of spin configurations from the classical Boltzmann distribution at finite-$T$.
- Quasi-linear scaling dipole-dipole interactions via the fast Fourier transform (FFT) (Langevin mode only).
- Measurement of $\mathcal{S}(\mathbf{q}, \omega)$ structure factor data.
- Interactive visualizations of the 3D crystal structure and spin states (support for 3D structure factors is planned.)
- Distributed implementation of [parallel tempering](https://en.wikipedia.org/wiki/Parallel_tempering) to increase sampling efficiency.

A current limitation of Sunny is that it requires real-space dynamical
simulations to measure the structure factor, and this limits momentum-space
resolution. Support for linear spin wave theory and its SU(_N_) generalization
is **in progress**.
