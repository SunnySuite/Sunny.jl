---
title: 'Sunny.jl: A Julia Package for Spin Dynamics[^1]'
tags:
  - Julia 
  - magnetism
  - spin dynamics
  - condensed matter physics 
  - neutron scattering 
authors:
  - name: David Dahlbom 
    corresponding: true
    orcid: 0000-0002-0221-5086 
    affiliation: "1, 2"

  - name: Hao Zhang 
    orcid: 0000-0002-9799-9118 
    affiliation: "3"

  - name: Cole Miles 
    orcid: 0000-0002-5581-4226 
    affiliation: "4"

  - name: Sam Quinn 
    orcid: 0000-0002-5043-3934 
    affiliation: "5, 6"

  - name: Alin Niraula
    affiliation: "7"

  - name: Bhushan Thipe
    affiliation: "7"

  - name: Matthew Wilson 
    affiliation: "8"

  - name: Sakib Matin
    affiliation: "3"

  - name: Het Mankad 
    orcid: 0000-0003-3890-2379 
    affiliation: "9"

  - name: Steven Hahn 
    orcid: 0000-0003-3890-2379 
    affiliation: "9"

  - name: Daniel Pajerowski 
    orcid: 0000-0003-3890-2379 
    affiliation: "1"

  - name: Steve Johnston 
    orcid: 0000-0002-2343-0113 
    affiliation: "2"

  - name: Zhentao Wang
    orcid: 0000-0001-7442-2933
    affiliation: "10" 

  - name: Harry Lane
    orcid: 
    affiliation: "11"

  - name: Ying Wai Li
    orcid: 0000-0003-0124-8262
    affiliation: "12"

  - name: Xiaojian Bai 
    orcid: 0000-0002-3974-626X 
    affiliation: "7"

  - name: Martin Mourigal
    orcid: 0000-0003-2772-8440 
    affiliation: "5"

  - name: Cristian D. Batista 
    orcid: 0000-0003-1667-3667 
    affiliation: "2"

  - name: Kipton Barros 
    corresponding: true
    orcid: 0000-0002-1333-5972 
    affiliation: "3"


affiliations:
 - name: Neutron Scattering Division, Oak Ridge National Laboratory 
   index: 1
 - name: Department of Physics and Astronomy, University of Tennessee 
   index: 2
 - name: Theoretical Division and CNLS, Los Alamos National Laboratory 
   index: 3
 - name: Kodiak Robotics 
   index: 4
 - name: School of Physics, Georgia Institute of Technology
   index: 5
 - name: Department of Physics and Astronomy, Univeriy of California, Los Angeles
   index: 6
 - name: Department of Physics and Astronomy, Louisiana State University
   index: 7
 - name: X-Computational Physics Division, Los Alamos National Laboratory
   index: 8
 - name: Computer Science and Mathematics Division, Oak Ridge National Laboratory 
   index: 9
 - name: School of Physics, Zhejiang University
   index: 10 
 - name: Department of Physics and Astronomy, University of Manchester 
   index: 11 
 - name: Computer, Computational, and Statistical Sciences Division, Los Alamos National Laboratory
   index: 12
date: 12 December 2024
bibliography: paper.bib
---

# Summary

Sunny is a Julia package designed to serve the needs of the quantum magnetism
community. It supports the specification of a very broad class of spin models
and a diverse suite of numerical solvers. These include powerful methods for
simulating spin dynamics both in and out of equilibrium. Uniquely, it features a
broad generalization of classical and semiclassical approaches to SU(_N_)
coherent states, which is useful for studying systems exhibiting strong
spin-orbit coupling or local entanglement effects. Sunny also offers a
well-developed framework for calculating the dynamical spin structure factor,
enabling direct comparison with scattering experiments. Ease of use is a
priority, with tools for symmetry-guided modeling and interactive visualization. 

[^1]: This manuscript has been authored by UT-Battelle, LLC, under contract
DE-AC05-00OR22725 with the US Department of Energy (DOE). The US government
retains and the publisher, by accepting the article for publication,
acknowledges that the US government retains a nonexclusive, paid-up,
irrevocable, worldwide license to publish or reproduce the published form of
this manuscript, or allow others to do so, for US government purposes. DOE will
provide public access to these results of federally sponsored research in
accordance with the DOE Public Access Plan
(https://www.energy.gov/doe-public-access-plan).


# Statement of need

Progress in quantum magnetism depends on the development of accurate models of
magnetic materials. Scattering techniques, such inelastic neutron scattering
(INS) and resonant inelastic X-ray scattering (RIXS), are among the most
informative methods available for probing the dynamics of quantum magnets,
yielding the dynamical spin structure factor $\mathcal{S}(\mathbf{q},\omega)$ as
experimental output. To evaluate the validity of a hypothetical model, it is
necessary to calculate $\mathcal{S}(\mathbf{q}, \omega)$ theoretically. This is
generally an intractable problem that must be treated numerically or with
various approximation schemes. The difficulty of this step represents a
bottleneck in the development of accurate models and impedes the advancement of
our understanding of quantum materials.

The Sunny project is a collaborative effort among theorists, experimentalists
and computational scientists aimed at developing theoretical and numerical
methodologies for modeling realistic quantum magnets. The central product of
this effort is the Sunny software package, which makes recent theoretical
advances available in a form readily accessible to students and researchers.
Distinguishing features of Sunny include:

- Symmetry analysis tools that facilitate model specification, visualization and
  data retrieval.
- A suite of optimizers, Monte Carlo samplers, and spin dynamics solvers that
  can all be applied to the same system specification. 
- Implementation of the SU(_N_) coherent state formalism for classical and
  semiclassical calculations.
- An interface tailored toward the needs of scattering scientists.
- Code written entirely in Julia, a language that can achieve speeds comparable
  to C++ or Fortran while maintaining the usability of a scripting language. 
- A well documented codebase, an extensive collection of correctness tests, and
  a website featuring many tutorials.

There are a number of existing codes that can calculate
$\mathcal{S}(\mathbf{q},\omega)$ using linear spin wave theory (LSWT), some of
which have served as inspiration to the Sunny project [@rotter:2004;
@SpinWaveGenie; @petit:2016; @li:2024]. The symmetry analysis tools of SpinW in
particular have served as a model [@toth:2015]. There are also codes that
perform classical spin simulations using Landau-Lifshitz (LL) dynamics
[@muller:2019; @evans:2014]. Sunny is unique in offering both approaches and
generalizing them through a formalism based on SU(_N_) coherent states
[@muniz:2014; @zhang_batista:2021]. Sunny additionally permits completely
general single-ion anisotropies and coupling of multipolar moments; provides an
efficient implementation of long-range dipole-dipole interactions; automates the
application of a number of quantum renormalizations [@dahlbom:2023]; and offers
iterative solvers for efficient LSWT on large magnetic cells [@lane:2024].

The value of collecting all these tools together in a modern, easy-to-use
package is evidenced by the large number of publications that have already made
use of Sunny, a partial list of which is maintained on the GitHub wiki [@Sunny].
We note a number of experimental studies that have relied on Sunny for analysis
[@bai:2021; @lee:2023; @do:2023; @bai:2023; @kim:2023; @paddison:2024;
@scheie:2023; @park:2023a; @nagl:2024; @na:2024; @park:2024a; @park:2024b;
@park:2024c]; as well as theoretical and methodological works [@zhang:2023;
@zhang:2024; @dahlbom:2024a; @dahlbom:2024b]. Additional papers documenting the
theoretical and algorithmic advances that have enabled the development of Sunny
are discussed below. 


# Feature Overview

## Symmetry analysis

By unifying and extending existing open source frameworks for the symmetry
analysis of crystals -- including Spglib [@togo:2024], Brillouin.jl
[@Brillouin], and CrystalInfoFramework.jl [@CrystalInfoFramework] -- Sunny
facilitates the process of determining the complete set of interactions allowed
by spacegroup symmetries. Similarly, any interaction specified on a site or bond
will be automatically propagated to all symmetry-equivalent sites and bonds, as
required by the spacegroup symmetries. Finally, the symmetry information enables
convenient specification of paths and slices through reciprocal space, aiding
visualization and comparison to experimental data. All these tools can be
applied just as easily to a user-specified crystal or to a crystal loaded from
an industry-standard CIF file [@hall:1991]. 

![a) Ground state of $\mathrm{FeI}_{2}$, found using Sunny's `minimize_energy!` function and visualized with `plot_spins`. b) The crystal of $\mathrm{FeI}_2$ visualized with the `view_crystal` function. Hovering the cursor over a bond reveals the exchange interaction, if already assigned, or a general expression for all symmetry-allowed interactions. \label{fig:symmetry}](figs/FigInteractions.png)


## Visualization 

Both the symmetry analysis and data retrieval features of Sunny include 3D
visualization tools built on the Makie package [@danisch:2021]. These can be
used to plot spin configurations, investigate the symmetries of a crystal
\autoref{fig:symmetry}, generate animations of dynamic behavior, and plot the
predicted results of scattering experiments \autoref{fig:Sqw}. 


## SU(_N_) Formalism and System Modes

Traditional classical and semiclassical approaches to spin dynamics are based on
the assignment of a classical dipole to each lattice site. Recent theoretical
work has generalized this picture, replacing dipoles with richer objects, namely
SU(_N_) coherent states. Within this formalism, quantum spin-$s$ is faithfully
represented as a linear combination of the $N = 2s + $ possible levels.
Capturing such local quantum effects is particularly important for describing
systems characterized by strong onsite anisotropies or local entanglement
effects. The SU(_N_) generalization applies equally to LSWT calculations
[@muniz:2014] and classical spin dynamics [@zhang_batista:2021]. Users can
access this formalism simply by setting the "mode" of a spin system to `:SUN`.
Sunny also offers a `:dipole` mode, which is similar to the traditional
classical approach but includes quantum renormalizations of biquadratic and
single-ion anisotropy terms [@dahlbom:2023]. Finally, there is a mode that
implements the traditional approach without any additional corrections,
`:dipole_uncorrected`. Most Sunny features are supported in all modes. 

![_Left_: Scattering intensities of $\mathrm{FeI}_2$ as measured on the SEQUOIA instrument at the Spallation Neutron Source, Oak Ridge National Laboratory [@bai:2021]. _Right_: Predicted scattering intensities calculated with Sunny's SU(_N_) linear spin wave solver. The figure was generated with Sunny's data retrieval and plotting functions. \label{fig:Sqw}](figs/FigSqw.png)

## Optimization and Monte Carlo Tools 

Identifying a classical ground state is often the first step when calculating
the scattering response of a magnet. Sunny provides several tools for finding
such states, including gradient optimizers built on the Optim.jl package
[@mogensen:2018]. These optimizers work on supercells as well as spiral
orderings. Sophisticated Monte Carlo tools are also provided, which can be used
both to anneal into ground states and to estimate finite temperature statistics.
In particular, the classical dynamics can be run with Langevin coupling to a
thermal bath [@dahlbom:2022b], and samplers are provided that implement the
Wang-Landau [@wang:2001] and parallel tempering algorithms [@swendsen:1986].


## Linear Spin Wave Theory (LSWT)

Sunny has extensive support for LSWT calculations, including for systems with
arbitrarily complex single-ion anisotropies, general bilinear and biquadratic
interactions, and long-range dipole-dipole interactions. Like SpinW, Sunny
provides efficient LSWT calculations for systems that exhibit incommensurate
spiral orderings [@toth:2015]. Additionally, Sunny provides tools to efficiently
calculate $\mathcal{S}(\mathbf{q},\omega)$ on very large magnetic cells using
iterative matrix-vector multiplications [@lane:2024]. The simulation of large
supercells is essential to study systems with chemical disorder or complex
magnetic orderings.


## Classical Dynamics

The efficiency of LSWT calculations makes it the preferred tool when studying
magnets near zero temperature. At elevated temperatures or in out-of-equilibrium
conditions, however, classical dynamics becomes a valuable technique. Sunny
supports both traditional Landau-Lifshitz dynamics and its generalization to
SU(_N_) coherent states [@zhang_batista:2021]. Dissipationless trajectories are
calculated using a symplectic integration scheme [@dahlbom:2022a], and a
generalization of the stochastic Landau-Lifshitz-Gilbert equations to SU(_N_)
coherent states [@dahlbom:2022b] enables the simulation of dynamics coupled to a
thermal bath. This is particularly valuable for simulating, e.g., thermal
transport, pump-probe experiments, and spin-glass relaxation.


## Sunny as a Platform for Future Developments
To make these existing features more widely available, work at ORNL is underway
to integrate Sunny into the Calvera platform for neutron data analysis
[@watson:2022]. Sunny itself can serve as a platform for new solvers and
analysis techniques, building on its mature model specification and data
retrieval features. Current efforts are directed at supporting: the
self-consistent Gaussian approximation for diffuse scattering, enabling
functionality inspired by [@paddison:2024]; non-perturbative corrections to LSWT
for the modeling of continua and bound states, which can be probed in INS and
terahertz spectroscopy experiments [@bai:2023; @legros:2021]; and observables
relevant to RIXS experiments. 



# Acknowledgements

We thank Mosé Giordano and Simon Danisch for valuable discussions.
This work was supported by the U.S. Department of Energy, Office of Science,
Office of Basic Energy Sciences, under Award Numbers DE-SC0022311,
DE-SC-0018660, and DE-SC0025426. Support was also provided by 
the LANL LDRD program. C.D.B. acknowledges partial support from the National
Science Foundation Materials Research Science and Engineering Center program
through the UT Knoxville Center for Advanced Materials and Manufacturing
(DMR-2309083).

# References