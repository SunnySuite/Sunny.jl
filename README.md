<div align="center">
    <a href="https://github.com/SunnySuite/Sunny.jl/">
        <picture>
            <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo-dark.svg">
            <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.svg" alt="Sunny.jl" width="350px">
        </picture>
    </a>
    <br><br>
    <table>
    <tr>
        <td>
            <b>Documentation:</b>&nbsp;&nbsp;
            <a href="https://sunnysuite.github.io/Sunny.jl/stable"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Docs-stable"></a>&nbsp;&nbsp;
            <a href="https://sunnysuite.github.io/Sunny.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Docs-dev"></a>
        </td>
        <td>
            <b>Build status:</b>&nbsp;&nbsp;
            <a href="https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="CI"></a>
        </td>
    </tr>
    </table>
</div>

## Overview

Sunny is a Julia package for modeling atomic-scale magnetism using classical spin dynamics with quantum corrections. It provides powerful tools to estimate dynamical structure factor intensities, $\mathcal{S}(𝐪,ω)$, enabling quantitative comparison with experimental scattering data, e.g., neutrons or x-rays.

A unique feature of Sunny is its treatment of spins as [SU(_N_) coherent states](https://doi.org/10.48550/arXiv.2106.14125). Each quantum spin-_S_ state is a superposition of $N=2S+1$ levels, and evolves under unitary, SU(_N_) transformations. When entanglement is neglected, the formalism allows to generalize the Landau-Lifshitz dynamics of spin dipoles to a dynamics of spin multipoles. The theory becomes especially useful for modeling materials with strong single-ion anisotropy effects (see our [FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/dev/examples/01_LSWT_SU3_FeI2.html)). In the future, Sunny will also support explicit spin-orbit coupling, and 'units' of locally entangled spins.

At low-temperatures, Sunny supports the usual linear spin wave theory for spin dipoles, and its ['multi-boson' generalization](https://doi.org/10.48550/arXiv.1307.7731). At finite temperatures, the full classical dynamics (with quantum correction factors) may be preferable to capture strongly nonlinear fluctuations. The [coupling of SU(_N_) spins to a thermal bath](https://doi.org/10.48550/arXiv.2209.01265) also makes possible the study of various non-equilibrium dynamics, e.g., thermal transport, pump-probe experiments, and spin-glass relaxation.

Sunny provides a number of tools to facilitate the specification and solution of spin Hamiltonians. This includes spacegroup symmetry analysis, powerful Monte Carlo sampling algorithms, and interactive 3D visualization. Efficient simulation is made possible by several algorithmic developments. See our [Sunny publications](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature) for more information.


## Try it out!

New Julia users should begin with our **[Getting Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia)** guide. A showcase for Sunny is the **[FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/dev/examples/01_LSWT_SU3_FeI2.html)**. This compound has a strong easy-axis anisotropy which gives rise to a single-ion bound state, and serves to exemplify the power of simulating SU(_N_) coherent states. Regarding the traditional linear spin wave theory of dipoles, see our [adaptations of the SpinW tutorials](https://sunnysuite.github.io/Sunny.jl/dev/examples/spinw/SW08_Kagome_AFM.html).

Sunny is evolving rapidly. [Version History](https://sunnysuite.github.io/Sunny.jl/dev/versions/) lists the new features and breaking changes. To install a specific version of Sunny, say `v0.x`, use the command `add Sunny@0.x`.

## Other spin wave codes

Sunny is inspired by SpinW, especially regarding symmetry analysis, model specification, and linear spin wave calculations. Relative to other spin wave codes, this table highlights Sunny's special features (as of 2024):

| | [McPhase](https://github.com/mducle/mcphase) | [SpinW](https://github.com/SpinW/spinw) | Sunny |
| -- | -- | -- | -- |
| Symmetry-guided modeling | ❌ | ✅ | ✅ |
| Interactive graphics | ❌ | ✅ | ✅ |
| General spin-multipole interactions | ✅ | ❌ | ✅ |
| [Interaction renormalization for dipoles](https://arxiv.org/abs/2304.03874) | ❌ | ❌ | ✅ |
| [Multi-flavor spin wave theory](https://arxiv.org/abs/1307.7731) | ✅ | ❌ | ✅ |
| [Classical SU(_N_) spin dynamics](https://arxiv.org/abs/2209.01265)</u> | ❌ | ❌ | ✅ |
| Ewald summation for dipole-dipole | ❌ | ❌ | ✅ |
| Programming language | C++ | Matlab | [Julia](https://julialang.org/) |

The classical SU(_N_) spin dynamics in Sunny generalizes the Landau-Lifshitz equation for $S > 1/2$ quantum spins. Linearizing and quantizing SU(_N_) dynamics yields a generalization of spin wave theory involving multi-flavor bosons.

Codes like [Spirit](https://github.com/spirit-code/spirit) and [Vampire](https://vampire.york.ac.uk/) focus less on capturing quantum effects, but may be good options for classical dynamics of pure dipoles.

## Join our community

We want to interact with you! Please **[join our Slack community](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA)** and say hello. If you encounter a problem with Sunny, please ask on the Slack `#helpdesk` channel or  [file a Github issue](https://github.com/SunnySuite/Sunny.jl/issues). If you use Sunny in a paper, please let us know and add it to our [Literature Wiki](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).

<br>
<div>
    <a href="https://www.lanl.gov">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="assets/lanl-dark.svg">
        <img src="assets/lanl-light.svg" alt="LANL" height="38px">
    </picture>
    </a> &nbsp;&nbsp;
    <a href="https://www.utk.edu">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="assets/utk-dark.svg">
        <img src="assets/utk-light.svg" alt="UTK" height="38px">
    </picture>
    </a> &nbsp;&nbsp;
    <a href="https://www.gatech.edu">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="assets/gatech-dark.svg">
        <img src="assets/gatech-light.svg" alt="GATech" height="38px">
    </picture>
    </a> &nbsp;
    <a href="https://www.ornl.gov/">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="assets/ornl-dark.svg">
        <img src="assets/ornl-light.svg" alt="ORNL" height="38px">
    </picture>
    </a> &nbsp;&nbsp;
    <a href="https://www.lsu.edu/">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="assets/lsu-dark.svg">
        <img src="assets/lsu-light.svg" alt="LSU" height="38px">
    </picture>
    </a>
</div>

