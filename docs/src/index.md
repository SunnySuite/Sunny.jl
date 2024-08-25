# Overview

[Sunny](https://github.com/SunnySuite/Sunny.jl/) is a Julia package for modeling
atomic-scale magnetism. It provides powerful tools to study equilibrium and
non-equilibrium magnetic phenomena. In particular, it allows estimation of
dynamical structure factor intensities, $\mathcal{S}(𝐪,ω)$, to support
quantitative modeling of experimental scattering data.

**Features include**:

- Generalized spin dynamics using [SU(_N_) coherent
  states](https://arxiv.org/abs/2209.01265).
- Ability to specify a crystal from a `.cif` file or its spacegroup symmetry.
  Magnetic structures can be read from `.mcif` files.
- Interactive visualizations of the 3D crystals and magnetic ordering.
- Symmetry analysis to classify allowed interaction terms, and to propagate them
  by symmetry.
- Single-ion anisotropy at arbitrary order, which can be specified using Stevens
  operators or as a polynomial of spin operators.
- Monte Carlo sampling of spin configurations in thermal equilibrium, and
  optimization tools.
- Measurements of dynamical correlations. At low temperature, one can use linear
  spin wave theory and its multi-boson generalization. This generalizes to
  finite temperatures using the classical dynamics, which allows for strongly
  nonlinear effects.
- Long-range dipole-dipole interactions accelerated with the fast Fourier
  transform (FFT).
- Support for comparison with experimental data: form factor, dipole factor,
  temperature-dependent classical-to-quantum factors, etc.
