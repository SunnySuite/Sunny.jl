# Welcome

[Sunny](https://github.com/SunnySuite/Sunny.jl/) provides powerful tools to
study equilibrium and non-equilibrium magnetic phenomena. For example, it
facilitates calculation of dynamical structure factor intensities
$\mathcal{S}(𝐪,ω)$ that can be directly compared to experimental scattering
data.

## Look around

A good place to start is the [CoRh₂O₄ example](@ref "1. Spin wave simulations of
CoRh₂O₄"). This and subsequent tutorials demonstrate a range of Sunny features.
One can also browse the [SpinW ports](@ref "SW01 - FM Heisenberg chain"), which
focus on spin wave theory.

## Join our community

We want to interact with you! Please [join our Slack
community](https://join.slack.com/t/sunny-users/shared_invite/zt-1otxwwko6-LzPtp7Fazkjx2XEqfgKqtA)
and say hello. If you encounter a problem with Sunny, please ask on the Slack
`#helpdesk` channel. If you use Sunny in a paper, please add it to our
[Literature Wiki](https://github.com/SunnySuite/Sunny.jl/wiki/Sunny-literature).


## Why choose Sunny?

Sunny is powerful and easy to use. Features include:

- Ability to [specify a crystal](@ref Crystal) from a `.cif` file, from its
  spacegroup and Wyckoffs, or with symmetries inferred from the chemical cell.
  Magnetic structures [can be read](@ref set_dipoles_from_mcif!) from `.mcif`
  files.
- Interactive visualization of [3D crystals](@ref view_crystal) and [magnetic
  ordering](@ref plot_spins).
- Symmetry analysis to [classify allowed interaction terms](@ref
  print_symmetry_table), and to propagate them by symmetry.
- Single-ion anisotropy at arbitrary order, which [can be specified](@ref
  set_onsite_coupling!) using [Stevens operators](@ref stevens_matrices) or
  [spin polynomials](@ref spin_matrices). Arbitrary coupling between spin
  multipoles is [also allowed](@ref set_pair_coupling!).
- Statistical sampling of spins in thermal equilibrium using [Langevin
  dynamics](@ref Langevin) and [local Monte Carlo updates](@ref LocalSampler).
- [Fast optimization](@ref minimize_energy!) of the magnetic ground state.
- Support for non-equilibrium spin dynamics, with [generalization to spin
  multipoles](@ref "6. Dynamical quench into CP² skyrmion liquid") via the
  theory of [SU(_N_) coherent states](https://arxiv.org/abs/2209.01265).
- Dynamical correlation measurement via [linear spin wave theory](@ref
  SpinWaveTheory) and its multi-boson generalization. Special support is
  provided for calculations on [incommensurate spiral phases](@ref
  SpinWaveTheorySpiral) and on [large, disordered magnetic cells](@ref
  SpinWaveTheoryKPM). At finite temperatures, one can use classical dynamics to
  [sample dynamical correlations](@ref SampledCorrelations) with strong
  nonlinearities.
- Long-range [dipole-dipole interactions](@ref enable_dipole_dipole!)
  accelerated with the fast Fourier transform (FFT).
- Conveniences for comparing to experimental data: [form factors](@ref
  FormFactor), [custom spin contractions](@ref ssf_custom_bm),
  classical-to-quantum [renormalization factors](@ref "Interaction
  Renormalization"), etc.
