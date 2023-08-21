# Sunny.jl

[Sunny](https://github.com/SunnySuite/Sunny.jl/) is a Julia package for modeling atomic-scale magnetism using semi-classical spin dynamics. It provides powerful tools to estimate dynamical structure factor intensities, $\mathcal{S}(𝐪,ω)$, enabling quantitative comparison with experimental scattering data, e.g., neutrons or x-rays.

Features include:

- Generalized spin dynamics using [SU(_N_) coherent states](https://arxiv.org/abs/2209.01265).
- Ability specify a crystal by a `.cif` file, or using its spacegroup symmetry.
- Interactive visualizations of the 3D crystals and magnetic ordering.
- Symmetry analysis to classify allowed interaction terms, and to propagate them by symmetry.
- Single-ion anisotropy at arbitrary order, which can be specified using Stevens operators or as a polynomial of spin operators.
- Monte Carlo sampling of spin configurations in thermal equilibrium, and optimization tools.
- Measurements of dynamical correlation functions. For small supercells at low temperature, one can use linear spin wave theory and its multi-boson generalization. Alternatively, one can use the full semi-classical dynamics to study systems with large supercells (e.g., disordered systems), or anharmonic effects with thermal fluctuations.
- Long-range dipole-dipole interactions accelerated with the fast Fourier transform (FFT).
- Various correction factors to facilitate comparison with experimental data (form factor, dipole factor, temperature-dependent classical-to-quantum factors, intensity binning, etc.).
