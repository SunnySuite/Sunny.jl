# # 1. Spin wave simulations of CoRh₂O₄
#
# This tutorial introduces Sunny through its features for performing
# conventional spin wave theory calculations. For concreteness, we consider the
# crystal CoRh₂O₄ and reproduce the calculations of [Ge et al., Phys. Rev. B 96,
# 064413](https://doi.org/10.1103/PhysRevB.96.064413).

# ### Get Julia and Sunny
# 
# Sunny is implemented in Julia, which allows for interactive development (like
# Python or Matlab) while also providing high numerical efficiency (like C++ or
# Fortran). New Julia users should begin with our **[Getting
# Started](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia)**
# guide. Sunny requires Julia 1.10 or later.
#
# From the Julia prompt, load both `Sunny` and `GLMakie`. The latter is needed
# for graphics.

using Sunny, GLMakie

# If these packages are not yet installed, Julia will offer to install them for
# you. If executing this script gives an error, you may need to `update` Sunny
# and GLMakie from the [built-in package
# manager](https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia#the-built-in-julia-package-manager).

# ### Units

# The [`Units`](@ref) object defines physical constants for conversions. Select
# meV as the default energy scale and angstrom as the default length scale.

units = Units(:meV, :angstrom);

# ### Crystal cell
#
# A crystallographic cell may be loaded from a `.cif` file, or can be specified
# from atom positions and types.
#
# Start by defining the shape of the conventional chemical cell. CoRh₂O₄ has
# cubic spacegroup 227 (Fd-3m). Its lattice constants are 8.5 Å, and the cell
# angles are 90°. With this information, [`lattice_vectors`](@ref) constructs a
# 3×3 matrix `latvecs`. Columns of `latvecs` define the lattice vectors ``(𝐚_1,
# 𝐚_2, 𝐚_3)`` in the global Cartesian coordinate system. Conversely, columns
# of `inv(latvecs)` define the global Cartesian axes ``(\hat{x}, \hat{y},
# \hat{z})`` in components of the lattice vectors.

a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)

# Construct the [`Crystal`](@ref) cell from the spacegroup number 227 and one
# representative atom of each occupied Wyckoff. In the standard setting of
# spacegroup 227, position `[0, 0, 0]` belongs to Wyckoff 8a, which is the
# diamond cubic crystal.

positions = [[0, 0, 0]]
cryst = Crystal(latvecs, positions, 227; types=["Co"], setting="1")

# [`view_crystal`](@ref) launches an interface for interactive inspection and
# symmetry analysis.

view_crystal(cryst)

# ### Spin system

# A [`System`](@ref) will define the spin model. This requires
# [`SpinInfo`](@ref) information for one representative atom per
# symmetry-distinct site. The cobalt atoms have quantum spin ``S = 3/2``. The
# ``g``-factor defines the magnetic moment ``μ = g 𝐒`` in units of the Bohr
# magneton. The option `:dipole` indicates a traditional model type, for which
# quantum spin is modeled as a dipole expectation value.

latsize = (1, 1, 1)
S = 3/2
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole)

# Previous work demonstrated that inelastic neutron scattering data for CoRh₂O₄
# is well described with a single antiferromagnetic nearest neighbor exchange,
# `J = 0.63` meV. Use [`set_exchange!`](@ref) with the bond that connects atom 1
# to atom 3, and has zero displacement between chemical cells. Applying the
# symmetries of spacegroup 227, Sunny will propagate this interaction to the
# other nearest-neighbor bonds. Calling [`view_crystal`](@ref) with `sys` now
# shows the antiferromagnetic Heisenberg interactions as blue polkadot spheres.

J = +0.63 # (meV)
set_exchange!(sys, J, Bond(1, 3, [0, 0, 0]))
view_crystal(sys)

# ### Optimizing spins

# To search for the ground state, call [`randomize_spins!`](@ref) and
# [`minimize_energy!`](@ref) in sequence. For this problem, optimization
# converges rapidly to the expected Néel order. See this with
# [`plot_spins`](@ref), where spins are colored according to their global
# ``z``-component.

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[s[3] for s in sys.dipoles])

# The diamond lattice is bipartite, allowing each spin to perfectly anti-align
# with its 4 nearest-neighbors. Each of these 4 bonds contribute ``-JS^2`` to
# the total energy. Two sites participate in each bond, so the energy per site
# is ``-2JS^2``. Check this by calling [`energy_per_site`](@ref).

@assert energy_per_site(sys) ≈ -2J*S^2

# ### Reshaping the magnetic cell

# The same Néel order can also be described with a magnetic cell that consists
# of the 2 cobalt atoms in the primitive cell. Columns of the 3×3 `shape` matrix
# below are the primitive lattice vectors in units of the conventional, cubic
# lattice vectors ``(𝐚_1, 𝐚_2, 𝐚_3)``. Use [`reshape_supercell`](@ref) to
# construct a system with this shape, and verify that the energy per site is
# unchanged.

shape = [0 1 1;
         1 0 1;
         1 1 0] / 2
sys_prim = reshape_supercell(sys, shape)
@assert energy_per_site(sys_prim) ≈ -2J*S^2

# Plotting the spins of `sys_prim` shows the primitive cell as a gray wireframe
# inside the conventional cubic cell.

plot_spins(sys_prim; color=[s[3] for s in sys_prim.dipoles])

# ### Spin wave theory

# With this primitive cell, we will perform a [`SpinWaveTheory`](@ref)
# calculation of the structure factor ``\mathcal{S}(𝐪,ω)``. The measurement
# [`ssf_perp`](@ref) indicates projection of the spin structure factor
# perpendicular to the direction of momentum transfer. This measurement is
# appropriate for unpolarized neutron scattering.

swt = SpinWaveTheory(sys_prim; measure=ssf_perp(sys_prim))

# Define a [`q_space_path`](@ref) that connects high-symmetry points in
# reciprocal space. The ``𝐪``-points are given in reciprocal lattice units
# (RLU) for the _original_ cubic cell. For example, `[1/2, 1/2, 0]` denotes the
# sum of the first two reciprocal lattice vectors, ``𝐛_1/2 + 𝐛_2/2``. A total
# of 500 ``𝐪``-points will be sampled along the path.

qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 500)

# Select [`lorentzian`](@ref) broadening with a full-width at half-maximum
# (FWHM) of 0.8 meV. The isotropic [`FormFactor`](@ref) for Co²⁺ dampens
# intensities at large ``𝐪``.

kernel = lorentzian(fwhm=0.8)
formfactors = [FormFactor("Co2")];

# Calculate the single-crystal scattering [`intensities`](@ref)` along the path,
# for 300 energy points between 0 and 6 meV. Use [`plot_intensities`](@ref) to
# visualize the result.

energies = range(0, 6, 300)
res = intensities(swt, path; energies, kernel, formfactors)
plot_intensities(res; units)

# To directly compare with the available experimental data, perform a
# [`powder_average`](@ref) over all possible crystal orientations. Consider 200
# ``𝐪`` magnitudes ranging from 0 to 3 inverse angstroms. Each magnitude
# defines spherical shell in reciprocal space, to be sampled with `2000`
# ``𝐪``-points. The calculation completes in just a couple seconds because the
# magnetic cell size is small.

radii = range(0, 3, 200) # (1/Å)
res = powder_average(cryst, radii, 2000) do qs
    intensities(swt, qs; energies, kernel, formfactors)
end
plot_intensities(res; units, saturation=1.0)

# This result can be compared to experimental neutron scattering data
# from Fig. 5 of [Ge et al.](https://doi.org/10.1103/PhysRevB.96.064413)
# ```@raw html
# <img width="95%" src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/docs/src/assets/CoRh2O4_intensity.jpg">
# ```


# ### What's next?
#
# * For more spin wave calculations of this traditional type, one can browse the
#   [SpinW tutorials ported to Sunny](@ref "SW01 - FM Heisenberg chain").
# * Spin wave theory neglects thermal fluctuations of the magnetic order. Our
#   [next tutorial](@ref "2. Landau-Lifshitz dynamics of CoRh₂O₄ at finite *T*")
#   demonstrates how to sample spins in thermal equilibrium, and measure
#   dynamical correlations from the classical spin dynamics.
# * Sunny also offers features that go beyond the dipole approximation of a
#   quantum spin via the theory of SU(_N_) coherent states. This can be
#   especially useful for systems with strong single-ion anisotropy, as
#   demonstrated in our [tutorial on FeI₂](@ref "3. Multi-flavor spin wave
#   simulations of FeI₂").