<a href="https://github.com/SunnySuite/Sunny.jl/">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo-dark.svg">
        <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.svg" alt="Sunny.jl" width="350px">
    </picture>
</a>

| **Documentation**         | **Build Status**      | **Citation**            |
| :------------------------ | :-------------------- | :---------------------- |
| [![][docs-img]][docs-url] | [![][ci-img]][ci-url] | [![][doi-img]][doi-url] |

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://sunnysuite.github.io/Sunny.jl/stable
[ci-img]: https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/SunnySuite/Sunny.jl/actions/workflows/CI.yml?query=branch%3Amain
[doi-img]: https://img.shields.io/badge/DOI-10.48550-blue
[doi-url]: https://doi.org/10.21105/joss.08138

## Overview

Sunny is a Julia package for modeling magnetic materials. It emphasizes _symmetry-aware_ Hamiltonians, careful treatment of _quantum-spin_ degrees of freedom, and a toolkit for _quantitative comparison with scattering data_, e.g., neutrons or x-rays. Ease of use is a priority: Sunny is extensively documented, supports interactive visualization, and has model fitting capabilities.

## Try it out

Start with the [Tutorials](https://sunnysuite.github.io/Sunny.jl/stable/examples/01_LSWT_CoRh2O4). For traditional spin wave theory, jump to the [SpinW ports](https://sunnysuite.github.io/Sunny.jl/stable/examples/spinw/SW01_FM_Heseinberg_chain.html).

See [Getting Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia) for installation instructions. See [Release Notes](https://sunnysuite.github.io/Sunny.jl/dev/versions) for new features and breaking changes.

## Key features

Sunny supports many standard tools for modeling spin systems and introduces several unique ones.

- **Symmetry-guided modeling**, including enumeration of [symmetry-allowed couplings](https://sunnysuite.github.io/Sunny.jl/dev/examples/03_LSWT_SU3_FeI2.html#Symmetry-analysis) and propagation of interactions by symmetry equivalence.
- **General spin couplings**. [Arbitrary single-ion anisotropy](https://sunnysuite.github.io/Sunny.jl/dev/library.html#Sunny.set_onsite_coupling!) may be specified via Stevens operator expansion or as a general spin polynomial. [Arbitrary multipolar coupling](https://sunnysuite.github.io/Sunny.jl/dev/library.html#Sunny.set_pair_coupling!) between site-pairs is also supported.
- **Spin-wave theory** for calculating quantum spin excitations. This includes the usual dipole formalism ([CoRh₂O₄ SWT](https://sunnysuite.github.io/Sunny.jl/stable/examples/01_LSWT_CoRh2O4.html)) and its generalization to spin multipoles via multi-flavor bosons ([FeI₂ SU(3) SWT](https://sunnysuite.github.io/Sunny.jl/stable/examples/03_LSWT_SU3_FeI2.html)). Sunny also supports special solvers for incommensurate spiral order ([Ba₃NbFe₃Si₂O₁₄ spiral](https://sunnysuite.github.io/Sunny.jl/stable/examples/spinw/SW15_Ba3NbFe3Si2O14.html)) and for scalability to very large disordered magnetic cells ([KPM demo](https://sunnysuite.github.io/Sunny.jl/stable/examples/09_Disorder_KPM.html)).
- **Finite-temperature dynamics and sampling**. This includes the Landau-Lifshitz spin dynamics with Langevin coupling to a thermal bath ([CoRh₂O₄ dynamics](https://sunnysuite.github.io/Sunny.jl/stable/examples/02_LLD_CoRh2O4.html)) and its generalization to spin multipoles via SU(_N_) coherent states ([FeI₂ dynamics](https://sunnysuite.github.io/Sunny.jl/stable/examples/04_GSD_FeI2.html)). Monte Carlo methods such as parallel tempering accelerate the sampling of highly frustrated magnets ([Advanced MC demos](https://github.com/SunnySuite/Sunny.jl/tree/main/examples/extra/Advanced_MC)).
- **Self-consistent Gaussian approximation** [(SCGA)](https://sunnysuite.github.io/Sunny.jl/stable/library.html#Sunny.SCGA) for efficient paramagnetic-phase observables, e.g. susceptibility and diffuse scattering intensity.
- **Long-range dipole-dipole interactions** with proper Ewald summation and customizable demagnetization tensor ([Pyrochlore Ewald demo](https://sunnysuite.github.io/Sunny.jl/stable/examples/07_Dipole_Dipole.html)).
- **Model fitting** to experimental data including magnon bands ([LuVO₃ fitting](https://sunnysuite.github.io/Sunny.jl/dev/examples/spinw/SW35_LuVO3-fitting.html)) and diffuse scattering intensities ([MgCr₂O₄ fitting](https://sunnysuite.github.io/Sunny.jl/dev/examples/10_SCGA_fitting.html)).

Sunny is also a platform to connect [theoretical research](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature#methods) with a wide range of [magnetism applications](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature#applications).

Related packages include [SpinW](https://github.com/SpinW/spinw) (symmetry-guided spin wave theory), [McPhase](https://github.com/mducle/mcphase) (multi-flavor bosons for accurate treatment of spin multipoles), and [Spinteract](https://doi.org/10.48550/arXiv.2210.09016) (SCGA solvers and fitting). [Vampire](https://vampire.york.ac.uk/) might be a good choice for larger-scale micromagnetic simulation.

## Join our community

We want to interact with you! [Join our Slack workspace](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA) and say hello. If you encounter a problem, ask on the Slack `#helpdesk` channel. If you find Sunny useful, please cite the main [JOSS paper](https://doi.org/10.21105/joss.08138) and any relevant [methodology papers](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature#methods). Finally, share your work with others on the [Sunny Literature Wiki](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).

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


