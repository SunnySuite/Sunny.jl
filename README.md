<div align="center">
    <a href="https://github.com/SunnySuite/Sunny.jl/">
    <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.jpg" alt="Sunny.jl" width="350px">    
    </a>
</div>
<p>

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sunnysuite.github.io/Sunny.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sunnysuite.github.io/Sunny.jl/dev)
[![Build Status](https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml?query=branch%3Amain)

A package for simulating classical spin systems, including the Landau-Lifshitz dynamics of spin dipoles and its generalization to multipolar spin components. In the latter case, Sunny resolves the local quantum structure of individual spins, making it particularly suited for modeling magnetic compounds with strong single-ion anisotropy, e.g., due to crystal field splitting.

Sunny additionally provides Monte Carlo algorithms for sampling from thermal equilibrium, as well as tools for measuring dynamical structure factors that can be compared with experimental neutron scattering data. Symmetry analyses are available to facilitate the design and specification of model Hamiltonians. Some interactive 3D visualization tools are available, with more planned.

Sunny is currently under heavy development with many **breaking changes**. See the [version history](https://sunnysuite.github.io/Sunny.jl/dev/versions/) for details.

## Examples

To see Sunny in action, we recommend browsing the [FeI2 case study](https://sunnysuite.github.io/Sunny.jl/dev/examples/fei2_tutorial/). This example uses SU(3) spin dynamics to measure the dynamical structure factor. The SU(3) formalism is essential to model the single-ion bound-state, which has been [observed experimentally](https://doi.org/10.1038/s41567-020-01110-1).

In addition to the examples in the official [documentation](https://sunnysuite.github.io/Sunny.jl/dev/), a number of tutorials are available as Jupyter notebooks at the [SunnyTutorials](https://github.com/SunnySuite/SunnyTutorials/tree/main/Tutorials) repo. 

## Technical description of SU(_N_) spin dynamics.

A quantum spin of magnitude _S_ has _N_ = 2 _S_ + 1 distinct angular momentum eigenstates. The traditional classical approximation effectively replaces each spin operator with a dipole expectation value, producing the Landau-Lifshitz spin-dipole dynamics. Alternatively, one can derive a dynamics that retains all multipolar expectation values (dipoles, quadrupoles, etc.) [[Zhang and Batista, Phys. Rev. B **104**, 104409 (2021)](https://arxiv.org/abs/2106.14125)]. This formalism can be understood as a dynamics of SU(_N_) coherent states, and coincides with the usual Laudau-Lifshitz dynamics when _N_=2.

Sunny uses recently developed algorithms to simulate this SU(_N_) spin dynamics in a highly efficient way:
* D. Dahlbom et al., _Geometric integration of classical spin dynamics via a mean-field Schrödinger equation_, Phys. Rev. B **106**, 054423 (2022) [[arXiv:2204.07563](https://arxiv.org/abs/2204.07563)].
* D. Dahlbom et al., _Langevin dynamics of generalized spins as SU(N) coherent states_, Phys. Rev. B **106**, 235154 (2022) [[arXiv:2209.01265](https://arxiv.org/abs/2209.01265)].

## Comparison with other tools

In addition to the generalized SU(_N_) spin dynamics, Sunny can also simulate spins in the dipole-only approximation. When running in this mode, the capabilities of Sunny are similar to [SpinW](https://spinw.org/), which was an inspiration for this project. Sunny supports general single-ion anisotropies, and currently has a focus on classical spin dynamics. The classical spin dynamics allows for simulations at finite temperatures, and (soon) the application of chemical inhomogeneities. Support for the SU(_N_) version of linear spin wave theory is coming soon; this approach will be faster, and enables arbitrary _k_-space resolution.

## Installation

Sunny is implemented in the [Julia programming language](https://julialang.org/). New Julia users may wish to read our [Getting Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia-and-Sunny) wiki page.

From the Julia prompt, one can install the latest Sunny using the built-in package manager:
```
julia> ]
pkg> add Sunny
```

New point-releases will include **breaking changes**. To install a specific version of Sunny, say `v0.x`, use the command `add Sunny@0.x`.

## Contact us

If you're using Sunny, we'd like to interact with you. Please join our [Slack User Community](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA) and say hello! If you find an unexpected behavior in Sunny, you can also [start a Github discussion](https://github.com/SunnySuite/Sunny.jl/discussions) or [file an issue](https://github.com/SunnySuite/Sunny.jl/issues). If you write a paper using Sunny, please add it to [our Wiki](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).

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

