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

A unique feature of Sunny is its treatment of spins as [SU(_N_) coherent states](https://arxiv.org/abs/2106.14125). This theory generalizes Landau-Lifshitz spin dynamics to a [nonlinear Schrödinger equation](https://arxiv.org/abs/2204.07563), which retains $N=2s+1$ levels for each quantum spin-_s_ state. This generalization is important for models with strong single-ion anisotropy (see our **[FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/dev/examples/03_LSWT_SU3_FeI2.html)**) and for [localized "units" of strongly entangled spins](https://arxiv.org/abs/2405.16315).  Efficient simulation is enabled by several algorithmic developments as listed on our [publications page](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).

At low-temperatures, Sunny supports the usual linear spin wave theory for spin dipoles, and its ['multi-boson' generalization](https://arxiv.org/abs/1307.7731). At finite temperatures, Sunny can calculate the dynamical structure factor using classical spin dynamics with quantum corrections. Langevin coupling to a thermal bath additionally makes possible the study of various non-equilibrium dynamics, e.g., thermal transport, pump-probe experiments, and spin-glass relaxation.

Sunny provides a number of tools to facilitate the specification and solution of spin Hamiltonians. This includes spacegroup symmetry analysis, powerful Monte Carlo sampling algorithms, and interactive 3D visualization.

## Try it out!

[Download Julia](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia) and try the **[Sunny Tutorials](https://sunnysuite.github.io/Sunny.jl/dev/examples/01_LSWT_CoRh2O4)**. For traditional linear spin wave theory, see also the **[SpinW ports](https://sunnysuite.github.io/Sunny.jl/dev/examples/spinw/SW01_FM_Heseinberg_chain.html)**.

Sunny is evolving rapidly. [Version History](https://sunnysuite.github.io/Sunny.jl/dev/versions/) lists the new features and breaking changes. To install a specific version of Sunny, say `v0.x`, use the command `add Sunny@0.x`.

## Other spin wave codes

Sunny is inspired by SpinW, especially regarding symmetry analysis, model specification, and linear spin wave calculations. Relative to other spin wave codes, this table highlights Sunny's special features (as of 2024):

| | [McPhase](https://github.com/mducle/mcphase) | [SpinW](https://github.com/SpinW/spinw) | Sunny |
| -- | -- | -- | -- |
| Symmetry-guided modeling | ❌ | ✅ | ✅ |
| Interactive graphics | ❌ | ✅ | ✅ |
| [Incommensurate spiral order](https://arxiv.org/abs/1402.6069) | ❌ | ✅ | ✅ |
| Arbitrary couplings of quantum operators | ✅ | ❌ | ✅ |
| [Interaction renormalization corrections](https://arxiv.org/abs/2304.03874) | ❌ | ❌ | ✅ |
| [Multi-flavor spin wave theory](https://arxiv.org/abs/1307.7731) | ✅ | ❌ | ✅ |
| [Classical SU(_N_) spin dynamics](https://arxiv.org/abs/2209.01265) | ❌ | ❌ | ✅ |
| [Ewald summation for dipole-dipole](https://sunnysuite.github.io/Sunny.jl/dev/examples/07_Dipole_Dipole.html) | ❌ | ❌ | ✅ |
| Programming language | C++ | Matlab | [Julia](https://julialang.org/) |

Codes like [Spirit](https://github.com/spirit-code/spirit) and [Vampire](https://vampire.york.ac.uk/) focus less on capturing quantum effects, but may be better options for the dynamics of classical dipoles, e.g., in the micromagnetics context.

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

