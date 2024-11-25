
# Why Choose Sunny?

## Powerful and easy to use

Feature highlights include:

- Ability to [specify a crystal](@ref Crystal) from a `.cif` file, from its
  spacegroup number and Wyckoffs, or from a full chemical cell with
  automatically inferred symmetry operations. Magnetic structures [can be
  read](@ref set_dipoles_from_mcif!) from `.mcif` files.
- Interactive visualization of [3D crystals](@ref view_crystal) and [magnetic
  structures](@ref plot_spins).
- Symmetry analysis to determine [allowed anisotropies and interaction
  terms](@ref print_symmetry_table), and to propagate them by symmetry
  equivalence.
- Single-ion anisotropy [can be specified](@ref set_onsite_coupling!) using
  either [Stevens operators](@ref stevens_matrices) or [spin polynomials](@ref
  spin_matrices). Arbitrary coupling between spin multipoles is [also
  supported](@ref set_pair_coupling!). Classical-to-quantum [renormalization
  factors](@ref "Interaction Renormalization") enhance fidelity.
- [Fast optimization](@ref minimize_energy!) of magnetic structures using
  supercells or the [propagation vector formalism](@ref
  minimize_spiral_energy!).
- Statistical sampling of spins in thermal equilibrium using [Langevin
  dynamics](@ref Langevin) or [local Monte Carlo updates](@ref LocalSampler).
  Advanced Monte Carlo methods, such as [parallel
  tempering](https://github.com/SunnySuite/Sunny.jl/tree/main/examples/extra/Advanced_MC),
  for simulations of classical spin liquids and frustrated magnetism.
- Classical dynamics of spin dipoles and its generalization to SU(_N_) coherent
  states. One can [sample dynamical correlations](@ref SampledCorrelations) at
  finite temperature. The [CP² skyrmion example](@ref "6. Dynamical quench into
  CP² skyrmion liquid") illustrates a highly non-equilibrium quench process that
  depends crucially on spin quadrupole degrees of freedom.
- Generalized [linear spin wave theory](@ref SpinWaveTheory) (LSWT) for
  low-temperature spin dynamics. Special support is provided for efficient
  calculations on [incommensurate spiral phases](@ref SpinWaveTheorySpiral) and
  on [very large magnetic cells](@ref SpinWaveTheoryKPM). The [FeI₂
  example](@ref "3. Multi-flavor spin wave simulations of FeI₂") showcases LSWT
  with multi-flavor bosons. The [disordered system example](@ref "9. Disordered
  system with KPM") demonstrates acceleration for large system sizes.
- [Dipole-dipole interactions](@ref enable_dipole_dipole!) with full Ewald
  summation, as illustrated in the [pyrochlore LSWT example](@ref "7. Long-range
  dipole interactions"). Dipole-dipole interactions in classical dynamics are
  accelerated with the fast Fourier transform (FFT).
- Tools for comparing to experimental data: [form factors](@ref FormFactor),
  [custom spin contractions](@ref ssf_custom_bm), averaging over [powder](@ref
  powder_average) and [domain orientations](@ref domain_average), etc.
- Programmatic interface in the [Julia language](https://julialang.org/) for
  full flexibility and performance.

But still evolving:

- Sunny does not yet have GPU acceleration of classical spin dynamics. An
  alternative here might be [Vampire](https://vampire.york.ac.uk/).

## Advanced theory made accessible

Sunny is also a platform for [disseminating foundational
advances](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature#methods)
in quantum magnetism and computational methods. The theory of [SU(_N_) coherent
states](https://doi.org/10.1103/PhysRevB.104.104409) offers a group theoretic
framework to formulate alternative classical limits of a microscopic quantum
model. [New algorithms](https://doi.org/10.1103/PhysRevB.106.235154) enable
highly efficient simulation of spin-multipoles and beyond. For reasons not fully
understood, such classical limits can be [remarkably accurate at finite
temperatures](https://doi.org/10.1103/PhysRevB.109.014427) when used in
conjunction with appropriate renormalization schemes. The SU(_N_) picture also
suggests new geometric interpretations of quantum spin. This leads to a [deeper
understanding](https://arxiv.org/abs/2304.03874) of existing renormalization
schemes for traditional spin wave theory, and suggests a whole landscape in
which to search for [novel topological
states](https://doi.org/10.1038/s41467-023-39232-8). Ongoing Sunny research aims
to incorporate more quantum entanglement into the classical picture. Local units
of strongly entangled spins will soon be supported and show great promise for
cases like [dimerized ladders](https://doi.org/10.1103/PhysRevB.110.104403).
Longer term, Sunny also aims to include perturbative corrections beyond linear
spin wave theory, as well as a non-perturbative treatment of quantum bound
states.
