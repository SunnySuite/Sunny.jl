<div align="center">
    <a href="https://github.com/SunnySuite/Sunny.jl/">
        <picture>
            <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo-dark.svg">
            <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.svg" alt="Sunny.jl" width="350px">
        </picture>
    </a>
    <br>
    <a href="https://sunnysuite.github.io/Sunny.jl/stable"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Docs-stable"></a>
    <a href="https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="CI"></a>
</div>

## Overview

Sunny is a Julia package for modeling atomic-scale magnetism with quantum effects. Through spin dynamics simulations, it enables direct comparison with experimental scattering data, e.g., neutrons or x-rays. Ease of use is a priority, with tools for symmetry-guided modeling and interactive visualization.

At low-temperatures, Sunny supports the usual linear spin wave theory of spin dipoles with generalization to multi-flavor bosons. At finite temperatures, Sunny supports the classical Landau-Lifshitz spin dynamics with generalization to SU(_N_) coherent states. Such generalizations are useful for modeling the coupled dynamics of higher order spin multipoles (see, e.g., the [FeI₂ tutorial](https://sunnysuite.github.io/Sunny.jl/stable/examples/03_LSWT_SU3_FeI2.html)), and for capturing localized "units" of strongly entangled spins. Through dynamical coupling to a thermal bath, Sunny makes possible the study of non-equilibrium dynamics, e.g., thermal transport, pump-probe experiments, and spin-glass relaxation. Many of these features build on our team's [theoretical research](https://sunnysuite.github.io/Sunny.jl/stable/why.html#Advanced-theory-made-accessible).

## Try it out!

Start by browsing the **[Tutorials](https://sunnysuite.github.io/Sunny.jl/stable/examples/01_LSWT_CoRh2O4)**. For traditional linear spin wave theory, see also the [SpinW ports](https://sunnysuite.github.io/Sunny.jl/stable/examples/spinw/SW01_FM_Heseinberg_chain.html).

See [Getting Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia) for installation instructions. [Version History](https://sunnysuite.github.io/Sunny.jl/dev/versions) lists new features and breaking changes.

## Related packages

Sunny is inspired by SpinW, especially regarding symmetry analysis and traditional spin wave theory. Relative to other spin wave codes, this table highlights Sunny's special features (as of 2024):

| | [McPhase](https://github.com/mducle/mcphase) | [SpinW](https://github.com/SpinW/spinw) | Sunny |
| -- | -- | -- | -- |
| Symmetry-guided modeling | ❌ | ✅ | ✅ |
| Interactive graphics | ❌ | ✅ | ✅ |
| [Incommensurate spiral order](https://sunnysuite.github.io/Sunny.jl/stable/examples/spinw/SW15_Ba3NbFe3Si2O14.html) | ❌ | ✅ | ✅ |
| [Interaction renormalization](https://sunnysuite.github.io/Sunny.jl/stable/renormalization.html) | ❌ | ❌ | ✅ |
| [Multi-flavor spin wave theory](https://sunnysuite.github.io/Sunny.jl/stable/examples/03_LSWT_SU3_FeI2.html) | ✅ | ❌ | ✅ |
| [Arbitrary spin-multipole couplings](https://sunnysuite.github.io/Sunny.jl/stable/library.html#Sunny.set_pair_coupling!) | ✅ | ❌ | ✅ |
| [Classical SU(_N_) spin dynamics](https://sunnysuite.github.io/Sunny.jl/stable/examples/04_GSD_FeI2.html) | ❌ | ❌ | ✅ |
| [Linear-scaling spin wave theory](https://sunnysuite.github.io/Sunny.jl/stable/examples/09_Disorder_KPM.html) | ❌ | ❌ | ✅ |
| [Fast long-range dipole interactions](https://sunnysuite.github.io/Sunny.jl/stable/examples/07_Dipole_Dipole.html) | ❌ | ❌ | ✅ |
| Programming language | C++ | Matlab | [Julia](https://julialang.org/) |

Codes like [Spirit](https://github.com/spirit-code/spirit) and [Vampire](https://vampire.york.ac.uk/) focus less on capturing quantum effects, but might be preferred for large-scale classical spin dynamics, e.g., for micromagnetics.

## Join our community

We want to interact with you! Please [join our Slack community](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA) and say hello. If you encounter a problem, please ask on the Slack `#helpdesk` channel. If you use Sunny in a paper, please add it to our [Literature Wiki](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).

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

