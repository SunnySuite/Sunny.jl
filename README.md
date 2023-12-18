<div align="center">
    <a href="https://github.com/SunnySuite/Sunny.jl/">
        <picture>
            <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo-dark.svg">
            <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.svg" alt="Sunny.jl" width="350px">
        </picture>
    </a>
</div>
<p>

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sunnysuite.github.io/Sunny.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sunnysuite.github.io/Sunny.jl/dev)
[![Build Status](https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

Sunny is a Julia package for modeling atomic-scale magnetism using classical spin dynamics with quantum corrections. It provides powerful tools to estimate dynamical structure factor intensities, $\mathcal{S}(𝐪,ω)$, enabling quantitative comparison with experimental scattering data, e.g., neutrons or x-rays.

A unique feature of Sunny is its treatment of spins as [SU(_N_) coherent states](https://doi.org/10.48550/arXiv.2106.14125). Each quantum spin-_S_ state is a superposition of $N=2S+1$ levels, and evolves under unitary, SU(_N_) transformations. Through neglect of entanglement, the formalism allows to generalize the Landau-Lifshitz dynamics of spin dipoles to a dynamics of spin multipoles. The theory becomes especially useful for modeling materials with strong single-ion anisotropy effects (see our [FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/dev/examples/01_LSWT_SU3_FeI2.html)). In the future, the theory could also be used to model explicit spin-orbit coupling, or 'units' of locally entangled spins.

At low-temperatures, Sunny supports the usual linear spin wave theory for spin dipoles, and its ['multi-boson' generalization](https://doi.org/10.48550/arXiv.1307.7731). At finite temperatures, the full classical dynamics (with quantum correction factors) may be preferable to capture strongly nonlinear effects. The [coupling of SU(_N_) spins to a thermal bath](https://doi.org/10.48550/arXiv.2209.01265) also makes possible the study of various non-equilibrium dynamics, e.g., thermal transport, pump-probe experiments, and spin-glass relaxation.

Sunny provides a number of tools to facilitate the specification and solution of spin Hamiltonians. This includes spacegroup symmetry analysis, powerful Monte Carlo sampling algorithms, and interactive 3D visualization. Efficient simulation is made possible by [several algorithmic developments](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).


## Try it out!

A good starting point is our **[FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/dev/examples/01_LSWT_SU3_FeI2.html)**. This compound includes effective spin-1 moments with strong easy-axis anisotropy, and exemplifies the power of simulating SU(3) coherent states.
<!-- 
In addition to the examples in the official [documentation](https://sunnysuite.github.io/Sunny.jl/dev/), a number of tutorials are available as Jupyter notebooks at the [SunnyTutorials](https://github.com/SunnySuite/SunnyTutorials/tree/main/Tutorials) repo.  -->

Try it yourself by downloading [Julia](https://julialang.org/), and installing Sunny within Julia's built-in package manager. New Julia users can refer to our **[Getting Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia)** guide. 

Sunny is evolving rapidly. See [Version History](https://sunnysuite.github.io/Sunny.jl/dev/versions/) for new features and breaking changes. To install a specific version of Sunny, say `v0.x`, use the command `add Sunny@0.x`.

## Other spin wave codes

Sunny is inspired by SpinW, especially regarding symmetry analysis, model specification, and linear spin wave calculations. Relative to other spin wave codes, this table highlights Sunny's special features (as of 2023):

| | [McPhase](https://github.com/mducle/mcphase) | [SpinW](https://github.com/SpinW/spinw) | Sunny |
| -- | -- | -- | -- |
| Symmetry-guided modeling | ❌ | ✅ | ✅ |
| Interactive graphics | ❌ | ✅ | ✅ |
| General spin-multipole interactions | ✅ | ❌ | ✅ |
| [Interaction renormalization for dipoles](https://arxiv.org/abs/2304.03874) | ❌ | ❌ | ✅ |
| [Multi-flavor spin wave theory](https://arxiv.org/abs/1307.7731) | ✅ | ❌ | ✅ |
| [Classical SU(_N_) spin dynamics](https://arxiv.org/abs/2209.01265)</u> | ❌ | ❌ | ✅ |
| Ewald summation for dipole-dipole | ❌ | ❌ | 🟨⁽¹⁾ |
| Programming language | C++ | Matlab | [Julia](https://julialang.org/) |

_Fine print: (1) Dipole-dipole interactions currently supported in classical dyamics but not LSWT._

The classical SU(_N_) spin dynamics in Sunny generalizes the Landau-Lifshitz equation for $S > 1/2$ quantum spins. Linearizing and quantizing SU(_N_) dynamics yields a generalization of spin wave theory involving multi-flavor bosons.

Codes like [Spirit](https://github.com/spirit-code/spirit) and [Vampire](https://vampire.york.ac.uk/) focus less on capturing quantum effects, but may be good options for classical dynamics of pure dipoles.

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

