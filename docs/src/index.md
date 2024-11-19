# Welcome

[Sunny](https://github.com/SunnySuite/Sunny.jl/) provides powerful tools to
simulate equilibrium and non-equilibrium magnetic phenomena from microscopic
models. It also facilitates calculation of the dynamical spin structure factor
$\mathcal{S}(𝐪,ω)$ for direct comparison to experimental scattering data.

To get a feel for Sunny, start by browsing the [CoRh₂O₄ example](@ref "1. Spin
wave simulations of CoRh₂O₄"). This and subsequent examples demonstrate a range
of Sunny features. See also the [SpinW ports](@ref "SW01 - FM Heisenberg
chain"), which focus on spin wave theory.

## Why choose Sunny?

Sunny is powerful and easy to use. Features include:

- Ability to [specify a crystal](@ref Crystal) from a `.cif` file, from its
  spacegroup number and representative Wyckoff positions, or from automatically
  inferred symmetry operations. Magnetic structures [can be read](@ref
  set_dipoles_from_mcif!) from `.mcif` files.
- Interactive visualization of [3D crystals](@ref view_crystal) and [magnetic
  structures](@ref plot_spins).
- Symmetry analysis to determine [allowed anisotropies and interaction
  terms](@ref print_symmetry_table), and to propagate them by symmetry
  equivalence.
- Single-ion anisotropy [can be specified](@ref set_onsite_coupling!) using
  either [Stevens operators](@ref stevens_matrices) or [spin polynomials](@ref
  spin_matrices). Arbitrary coupling between spin multipoles is [also
  supported](@ref set_pair_coupling!). Classical-to-quantum [renormalization
  factors](@ref "Interaction Renormalization") are included to enhance fidelity.
- [Fast optimization](@ref minimize_energy!) of magnetic structures using
  supercells or the [propagation vector formalism](@ref
  minimize_spiral_energy!).
- Statistical sampling of spins in thermal equilibrium using [Langevin
  dynamics](@ref Langevin) or [local Monte Carlo updates](@ref LocalSampler).
  Advanced Monte Carlo methods such as [parallel
  tempering](https://github.com/SunnySuite/Sunny.jl/tree/main/examples/extra/Advanced_MC)
  are effective for simulations of classical spin liquids and frustrated
  magnetism.
- [Non-equilibrium spin dynamics](@ref "6. Dynamical quench into CP² skyrmion
  liquid") with generalization to arbitrary spin multipoles via the theory of
  [SU(_N_) coherent states](https://arxiv.org/abs/2209.01265).
- Dynamical correlation measurement via [linear spin wave theory](@ref
  SpinWaveTheory) and its multi-boson generalization. Special support is
  provided for calculations on [incommensurate spiral phases](@ref
  SpinWaveTheorySpiral) and on [large, disordered magnetic cells](@ref
  SpinWaveTheoryKPM). At finite temperatures, one can use classical dynamics to
  [sample dynamical correlations](@ref SampledCorrelations) with strong
  nonlinearities.
- Long-range [dipole-dipole interactions](@ref enable_dipole_dipole!)
  accelerated with the fast Fourier transform (FFT).
- Tools for comparing to experimental data: [form factors](@ref FormFactor),
  [custom spin contractions](@ref ssf_custom_bm), averaging over [powder](@ref
  powder_average) and [domain orientations](@ref domain_average), etc.

## Join our community

We want to interact with you! Please [join our Slack
community](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA)
and say hello. If you encounter a problem with Sunny, please ask on the Slack
`#helpdesk` channel. If you use Sunny in a paper, please add it to our
[Literature Wiki](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).
