<div align="center">
    <a href="https://github.com/SunnySuite/Sunny.jl/">
    <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.jpg" alt="Sunny.jl" width="350px">    
    </a>
</div>
<p>

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahub.com/docs/General/Sunny/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sunnysuite.github.io/Sunny.jl/dev)
[![Build Status](https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

Sunny is a Julia package for modeling atomic-scale magnetism using classical spin dynamics with quantum corrections. It provides powerful tools to estimate dynamical structure factor intensities, $\mathcal{S}(𝐪,ω)$, enabling quantitative comparison with experimental scattering data, e.g., neutrons or x-rays.

A unique feature of Sunny is its treatment of spins as [SU(_N_) coherent states](https://doi.org/10.48550/arXiv.2106.14125). Each quantum spin-_S_ state is a superposition of $N=2S+1$ levels, and evolves under unitary, SU(_N_) transformations. Through neglect of entanglement, the formalism allows to generalize the Landau-Lifshitz dynamics of spin dipoles to a dynamics of spin multipoles. The theory becomes especially useful for modeling materials with strong single-ion anisotropy effects (see our [FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/dev/examples/fei2_tutorial/)). In the future, the theory could also be used to model explicit spin-orbit coupling, or 'units' of locally entangled spins.

At low-temperatures, Sunny supports the usual linear spin wave theory (LSWT) for spin dipoles, and its ['multi-boson' generalization](https://doi.org/10.48550/arXiv.1307.7731). At finite temperatures, the full classical dynamics (with quantum correction factors) may be preferable to capture anharmonic effects. The [coupling of SU(_N_) spins to a thermal bath](https://doi.org/10.48550/arXiv.2209.01265) also makes possible the study of various non-equilibrium dynamics, e.g., thermal transport, pump-probe experiments, and spin-glass relaxation.

Sunny provides a number of tools to facilitate the specification and solution of spin Hamiltonians. This includes spacegroup symmetry analysis, powerful Monte Carlo sampling algorithms, and interactive 3D visualization. Efficient simulation is made possible by [several algorithmic developments](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).


## Try it out!

To see Sunny in action, a good starting point is our **[FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/dev/examples/fei2_tutorial/)**. This compound includes effective spin-1 moments with strong easy-axis anisotropy. The coupled dipole-quadrupole dynamics is efficiently described within the formalism of SU(3) coherent states, and is [crucial to explain neutron scattering data](https://doi.org/10.1038/s41567-020-01110-1).
<!-- 
In addition to the examples in the official [documentation](https://sunnysuite.github.io/Sunny.jl/dev/), a number of tutorials are available as Jupyter notebooks at the [SunnyTutorials](https://github.com/SunnySuite/SunnyTutorials/tree/main/Tutorials) repo.  -->

Try it yourself by downloading [Julia](https://julialang.org/), and installing Sunny within Julia's built-in package manager. New Julia users can refer to our **[Getting Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia)** guide. 

Sunny is evolving rapidly. See [Version History](https://sunnysuite.github.io/Sunny.jl/dev/versions/) for new features and breaking changes. To install a specific version of Sunny, say `v0.x`, use the command `add Sunny@0.x`.

## Related projects

Sunny is heavily inspired by the [SpinW](https://spinw.org/) code. In particular, Sunny's symmetry analysis, model specification, and LSWT features will be familiar to users of SpinW. Sunny differs from SpinW in its support for nonlinear classical spin dynamics. Such dynamics can be useful, e.g., to study thermal fluctuations, non-equilibrium dynamics, transport, or disorder in systems with large magnetic supercells.

Another LSWT code is [SpinWaveGenie](https://github.com/SpinWaveGenie/SpinWaveGenie) which is written in C++ and very fast. With future optimizations, Sunny aims to achieve comparable performance.

Conversely, there exist many powerful codes for studying classical spin dynamics, including [Spirit](https://github.com/spirit-code/spirit) and [Vampire](https://vampire.york.ac.uk/). Compared to these codes, Sunny puts more emphasis on capturing quantum effects of magnetic compounds.

A distinguishing feature of Sunny compared to previous codes is its support for simulating generalized spins via the theory of SU(_N_) coherent states.

## Join our community

If you're using Sunny, we'd like to interact with you. Please join our **[Slack community](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA)** and say hello. If you find an unexpected behavior in Sunny, you can also  [file an issue](https://github.com/SunnySuite/Sunny.jl/issues). If you write a paper using Sunny, please add it to our [Literature Wiki](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).

<div>
    <a href="https://www.lanl.gov">
    <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/lanl.png" alt="LANL" width="250px">
    </a>
    <a href="https://www.utk.edu">
    <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/utk.png" alt="U. Tennessee" width="250px">
    </a>
    <a href="https://www.gatech.edu/">
    <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/gatech.png" alt="Georgia Tech." width="250px">
    </a>
</div>

