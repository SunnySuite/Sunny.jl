<div align="center">
    <img src="https://raw.githubusercontent.com/MagSims/Sunny.jl/main/assets/sunny_logo.jpg" width=50% alt="Sunny.jl">
</div>
<p>

<!--- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sunnysuite.github.io/Sunny.jl/stable) --->

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sunnysuite.github.io/Sunny.jl/dev)

A general-purpose library for simulating classical spin systems, including the Landau-Lifshitz dynamics of spin dipoles and its generalization to multipole spin moments.
<!-- 
## Example notebooks

[**Coming soon**]  To get a feeling for what Sunny can do, we recommend to start by browsing through our Jupyter notebook examples. -->

## What it does

Sunny simulates the classical interactions between locally quantum spins. A quantum spin of magnitude _S_ is an $N = 2 S + 1$ level system. We call this local state an "SU(_N_) spin." Moving from SU(2) spin dipoles to the more general framework of SU(_N_) allows to capture all quantum entanglement between the local _N_ levels, which is crucial to modeling magnetic compounds with strong local anisotropy.

Sunny provides Monte Carlo algorithms for sampling SU(_N_) spin configurations from thermal equilibrium, as well as tools for measuring dynamical structure factors that can be compared with experimental neutron scattering data. Sunny provides symmetry analyses to facilitate the design and specification of model Hamiltonians, and interactive tools to visualize 3D crystal structures and (coming soon) structure factor data.

SU(_N_) spin dynamics was originally proposed in:

* H. Zhang and C. D. Batista, _Classical spin dynamics based on SU(N) coherent states_, Phys. Rev. B 104, 104409 (2021) [[arXiv:2106.14125](https://arxiv.org/abs/2106.14125)].

Sunny implements geometrically constrained numerical integration of this SU(_N_) dynamics:

* D. Dahlbom et al., _Geometric integration of classical spin dynamics via a mean-field Schrödinger equation_ [[arXiv:2204.07563](https://arxiv.org/abs/2204.07563)].

Sunny additionally implements many new algorithms for efficient equilibrium sampling of SU(_N_) spins.

_Dipole mode_. As a special case, Sunny can be restricted to the dipole-only approximation of spin. In this mode, the capabilities of Sunny are similar to [SpinW](https://spinw.org/). A key difference is that Sunny does not (currently) employ linear spin wave theory. Advantages are: (1) Applicability to finite temperature measurements and (2) Support for single-ion anisotropies beyond quadratic order.   A disadvantage is that structure factor measurements $\mathcal S(q,\omega)$ have momentum-space ($q$) resolution that is limited by the size of magnetic super cell.

## Installation

Sunny is implemented in the [Julia programming language](https://julialang.org/). New Julia users may wish to start with our [Getting Started](GettingStarted.md) guide.

From the Julia prompt, one can install Sunny using the built-in package manager. We currently recommend tracking the main branch:
```
julia> ]
pkg> add Sunny#main
```

Check that Sunny is working properly by running the unit tests: `pkg> test Sunny`. Please keep up-to-date by periodically running the Julia update command: `pkg> update`.

A good way to interact with Sunny is through the Jupyter notebook interface. This support can be installed through the [IJulia](https://github.com/JuliaLang/IJulia.jl) package.

## API Reference

[Full documentation available here](https://sunnysuite.github.io/Sunny.jl/dev).

## Contact us

If you discover bugs, or find Sunny useful, please contact us at kbarros@gmail.com and david.dahlbom@gmail.com.

<!-- Users who wish to contribute to Sunny source-code development should instead use the `dev` command:
```
julia> ]
pkg> dev Sunny
```

This will `git clone` the source code to the directory `~/.julia/dev/Sunny`. You can make changes to these files,
and they will be picked up by Julia.  The package manager will not touch
any package installed by `dev`, so you will be responsible
for keeping Sunny up to date, e.g., using the command `git pull` from Sunny package directory. -->


<!-- 
For plotting, you may also wish to install
```
pkg> add Plots
pkg> add GLMakie
```

At the time of this writing, GLMakie has some rough edges, especially on Mac platforms. Run `test GLMakie` to make sure it is working properly. -->

