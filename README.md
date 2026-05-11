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

Sunny is a Julia package for simulating magnetism. It emphasizes symmetry-aware spin Hamiltonians, efficient methods for capturing quantum effects, and comparison with experimental data, e.g., neutron or X-ray scattering. Sunny is extensively documented, supports interactive visualization, and facilitates model fitting.

## Try it out

Start with the [Tutorials](https://sunnysuite.github.io/Sunny.jl/stable/examples/01_LSWT_CoRh2O4). For traditional spin wave theory, jump to the [SpinW ports](https://sunnysuite.github.io/Sunny.jl/stable/examples/spinw/SW04_Frustrated_square.html).

See also the [Getting Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia) guide and [Release Notes](https://sunnysuite.github.io/Sunny.jl/stable/versions).

## Key features

- **Symmetry-guided modeling**. Load CIF or mCIF files, enumerate [symmetry-allowed couplings](https://sunnysuite.github.io/Sunny.jl/stable/examples/03_LSWT_SU3_FeI2.html#Symmetry-analysis), and propagate interactions by symmetry equivalence.
- **General spin couplings**. Specify [single-ion anisotropy](https://sunnysuite.github.io/Sunny.jl/stable/library.html#Sunny.set_onsite_coupling!) in Stevens operators or spin polynomials. Arbitrary [multipolar coupling](https://sunnysuite.github.io/Sunny.jl/stable/library.html#Sunny.set_pair_coupling!) on a bond is also supported.
- **Spin-wave theory** for quantum spin excitations, including the usual dipole formalism ([CoRh₂O₄ bands](https://sunnysuite.github.io/Sunny.jl/stable/examples/01_LSWT_CoRh2O4.html)) and its generalization to multipoles via multi-flavor bosons ([FeI₂ bands](https://sunnysuite.github.io/Sunny.jl/stable/examples/03_LSWT_SU3_FeI2.html)). Special calculators are available for incommensurate spiral order ([Ba₃NbFe₃Si₂O₁₄ spiral](https://sunnysuite.github.io/Sunny.jl/stable/examples/spinw/SW15_Ba3NbFe3Si2O14.html)) and for linear-scaling on very large systems ([KPM solver](https://sunnysuite.github.io/Sunny.jl/stable/library.html#Sunny.SpinWaveTheoryKPM)).
- **Finite-temperature spin dynamics** of dipoles via the Landau-Lifshitz equation with thermostat ([CoRh₂O₄ dynamics](https://sunnysuite.github.io/Sunny.jl/stable/examples/02_LLD_CoRh2O4.html)). This generalizes to a spin multipole dynamics via the theory of SU(_N_) coherent states ([FeI₂ dynamics](https://sunnysuite.github.io/Sunny.jl/stable/examples/04_GSD_FeI2.html)). Methods like parallel tempering accelerate equilibrium sampling ([Advanced MC]((https://github.com/SunnySuite/Sunny.jl/tree/main/examples/extra/Advanced_MC))).
- **Self-consistent Gaussian approximation** [(SCGA)](https://sunnysuite.github.io/Sunny.jl/stable/library.html#Sunny.SCGA) for static paramagnetic-phase observables, e.g. susceptibility and diffuse scattering.
- **Long-range dipole-dipole interactions** with proper Ewald summation and customizable demagnetization tensor ([Pyrochlore Ewald](https://sunnysuite.github.io/Sunny.jl/stable/examples/07_Dipole_Dipole.html)).
- **Chemical disorder** in large, inhomogeneous systems ([Disordered intensities](https://sunnysuite.github.io/Sunny.jl/stable/examples/09_Disorder_KPM.html)).
- **Entangled units** to capture the local entanglement of strongly-coupled sites ([Ba₃Mn₂O₈ intensities](https://sunnysuite.github.io/Sunny.jl/stable/examples/contributed/entangled_units.html)).
- **Model fitting** to experimental data, e.g., magnon bands ([LuVO₃ fitting](https://sunnysuite.github.io/Sunny.jl/stable/examples/spinw/SW35_LuVO3-fitting.html)), inelastic powder ([β-Na₂PrO₃ fitting](https://sunnysuite.github.io/Sunny.jl/stable/examples/11_Powder_fitting.html)), and diffuse scattering ([MgCr₂O₄ fitting](https://sunnysuite.github.io/Sunny.jl/stable/examples/10_SCGA_fitting.html)).

Sunny is also a platform for **methods development** that bridges theory and experiment ([team publications](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature#methods)).

Related packages include [SpinW](https://github.com/SpinW/spinw) (symmetry-guided spin wave theory), [McPhase](https://github.com/mducle/mcphase) (accurate treatment of spin multipoles), [Spinteract](https://doi.org/10.48550/arXiv.2210.09016) (SCGA solvers and fitting), and [Vampire](https://vampire.york.ac.uk/) (large-scale spin dynamics).

## Join our community

Say hello on [our Slack workspace](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA). Support is available on the `#helpdesk` channel.

## Citation

Please cite the main [JOSS paper](https://doi.org/10.21105/joss.08138) and any relevant [methodology papers](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature#methods) to help sustain Sunny development. 

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


