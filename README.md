<div align="center">
    <img src="https://raw.githubusercontent.com/MagSims/Sunny.jl/master/assets/sunny_logo.jpg" width=50% alt="Sunny.jl">
</div>
<p>

<!--- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sunnysuite.github.io/Sunny.jl/stable) --->

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sunnysuite.github.io/Sunny.jl/dev)

A general-purpose library for simulating classical spin systems, including the Landau-Lifshitz dynamics of spin dipoles and its generalization to multipole spin moments. A local quantum spin _S_ is an $N = 2 S + 1$ level system. We call this local state an "SU(_N_) spin." An SU(2) spin is simply a spin dipole. By generalizing to SU(_N_) spins, we capture all quantum entanglement between the local _N_ levels, which is crucial to modeling magnetic compounds with strong local anisotropy.

Sunny provides Monte Carlo algorithms for sampling SU(_N_) spin configurations from thermal equilibrium, as well as tools for measuring dynamical structure factors that can be compared with experimental neutron scattering data. Sunny provides symmetry analyses to facilitate the design and specification of model Hamiltonians, and interactive tools to visualize 3D crystal structures and (coming soon) structure factor data.

SU(_N_) spin dynamics was originally proposed in the paper:

* H. Zhang and C. D. Batista, _Classical spin dynamics based on SU(N) coherent states_, Phys. Rev. B 104, 104409 (2021) [[arXiv:2106.14125](https://arxiv.org/abs/2106.14125)].

Our preprint reviews this theory, and describes a numerical method for symplectic integration of SU(_N_) dynamics:

* D. Dahlbom et al., _Geometric integration of classical spin dynamics via a mean-field Schrödinger equation_ [[arXiv:2204.07563](https://arxiv.org/abs/2204.07563)].

## Getting started with Julia

New Julia users should begin with our [Getting Started](GettingStarted.md) guide.

## Installation

Sunny is evolving rapidly, and users are recommended to install the package by
tracking the main branch:
```
julia> ]
pkg> add Sunny#main
```
Please keep up-to-date by periodically running the Julia "`pkg> update`" command, and report any bugs that you encounter.

<!-- Users who wish to contribute to Sunny source-code development should instead use the `dev` command:
```
julia> ]
pkg> dev Sunny
```

This will `git clone` the source code to the directory `~/.julia/dev/Sunny`. You can make changes to these files,
and they will be picked up by Julia.  The package manager will not touch
any package installed by `dev`, so you will be responsible
for keeping Sunny up to date, e.g., using the command `git pull` from Sunny package directory. -->

Check that Sunny is working properly by running the unit tests, "`pkg> test Sunny`".

<!-- 
For plotting, you may also wish to install
```
pkg> add Plots
pkg> add GLMakie
```

At the time of this writing, GLMakie has some rough edges, especially on Mac platforms. Run `test GLMakie` to make sure it is working properly. -->

A good way to interact with Sunny is through the Jupyter notebook interface. This support can be installed through the [IJulia](https://github.com/JuliaLang/IJulia.jl) package.

## Documentation

[Full documentation available here](https://sunnysuite.github.io/Sunny.jl/dev).
