# # 2. Spin wave simulations of CoRh₂O₄
#
# This tutorial illustrates the conventional spin wave theory of dipoles. We
# consider a simple model of the diamond-cubic crystal CoRh₂O₄, with parameters
# extracted from [Ge et al., Phys. Rev. B 96,
# 064413](https://doi.org/10.1103/PhysRevB.96.064413).

using Sunny, GLMakie

# Construct a diamond [`Crystal`](@ref) in the conventional (non-primitive)
# cubic unit cell. Sunny will populate all eight symmetry-equivalent sites when
# given the international spacegroup number 227 ("Fd-3m") and the appropriate
# setting. For this spacegroup, there are two conventional translations of the
# unit cell, and it is necessary to disambiguate through the `setting` keyword
# argument. (On your own: what happens if `setting` is omitted?)

units = Units(:meV)
a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], 227, setting="1")

# In a running Julia environment, the crystal can be viewed interactively using
# [`view_crystal`](@ref).

view_crystal(cryst)

# Construct a [`System`](@ref) with quantum spin ``S=3/2`` constrained to the
# space of dipoles. Including an antiferromagnetic nearest neighbor interaction
# `J` will favor Néel order. To optimize this magnetic structure, it is
# sufficient to employ a magnetic lattice consisting of a single crystal unit
# cell, `latsize=(1,1,1)`. Passing an explicit random number `seed` will ensure
# repeatable results.

latsize = (1, 1, 1)
S = 3/2
J = 0.63 # (meV)
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=0)
set_exchange!(sys, J, Bond(1, 3, [0,0,0]))

# In the ground state, each spin is exactly anti-aligned with its 4
# nearest-neighbors. Because every bond contributes an energy of ``-JS^2``, the
# energy per site is ``-2JS^2``. In this calculation, a factor of 1/2 avoids
# double-counting the bonds. Due to lack of frustration, direct energy
# minimization is successful in finding the ground state.

randomize_spins!(sys)
minimize_energy!(sys)

@assert energy_per_site(sys) ≈ -2J*S^2

# Plotting the spins confirms the expected Néel order. Note that the overall,
# global rotation of dipoles is arbitrary.

s0 = sys.dipoles[1,1,1,1]
plot_spins(sys; color=[s'*s0 for s in sys.dipoles])

# For numerical efficiency, it is helpful to work with the smallest possible
# magnetic supercell; in this case, it is the primitive cell. The columns of the
# 3×3 `shape` matrix define the lattice vectors of the primitive cell as
# multiples of the conventional, cubic lattice vectors. After transforming the
# system with [`reshape_supercell`](@ref), the energy per site remains the same.
shape = [0 1 1;
         1 0 1;
         1 1 0] / 2
sys_prim = reshape_supercell(sys, shape)
@assert energy_per_site(sys_prim) ≈ -2J*S^2
plot_spins(sys_prim; color=[s'*s0 for s in sys_prim.dipoles])

# Now estimate ``𝒮(𝐪,ω)`` with [`SpinWaveTheory`](@ref).

swt = SpinWaveTheory(sys_prim; corrspec=ssf_perp(sys_prim))

# For the "single crystal" result, we use [`q_space_path`](@ref) to construct a
# path that connects high-symmetry points in reciprocal space.

qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 400)

# Select [`lorentzian2`](@ref) broadening with a full-width at half-maximum
# (FWHM) of 0.8 meV. Use [`ssf_perp`](@ref) to calculate unpolarized scattering
# intensities. The isotropic [`FormFactor`](@ref) for Cobalt(2+) dampens
# intensities at large ``𝐪``.

kernel = Sunny.lorentzian2(fwhm=0.8)
formfactors = [FormFactor("Co2")]
energies = range(0, 6, 300)
res = intensities(swt, path; energies, kernel, formfactors)
plot_intensities(res; units)

# We use [`powder_average`](@ref) to average intensities over all possible
# crystal orientations. Perform this calculation for 200 momentum magnitudes,
# ranging from 0 to 3 inverse angstroms. Each ``𝐪``-magnitude defines a
# spherical shell in reciprocal space. Sample it with `2000` wavevectors of
# quasi-uniform distribution.

radii = range(0, 3, 200) # (1/Å)
res = powder_average(cryst, radii, 2000) do qs
    intensities(swt, qs; energies, kernel, formfactors)
end
plot_intensities(res; units)

# This result can be compared to experimental neutron scattering data
# from Fig. 5 of [Ge et al.](https://doi.org/10.1103/PhysRevB.96.064413)
# ```@raw html
# <img width="95%" src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/docs/src/assets/CoRh2O4_intensity.jpg">
# ```
