# # SW18 - Distorted kagome
#
# This is a Sunny port of [SpinW Tutorial
# 18](https://spinw.org/tutorials/18tutorial), originally authored by Goran
# Nilsen and Sandor Toth. This tutorial illustrates spin wave calculations for
# KCu₃As₂O₇(OD)₃. The Cu ions are arranged in a distorted kagome lattice, and
# exhibit an incommensurate helical magnetic order, as described in [G. J.
# Nilsen, et al., Phys. Rev. B **89**, 140412
# (2014)](https://doi.org/10.1103/PhysRevB.89.140412). The model follows [Toth
# and Lake, J. Phys.: Condens. Matter **27**, 166002
# (2015)](https://arxiv.org/abs/1402.6069).


using Sunny, GLMakie

# Build the distorted kagome crystal, with spacegroup 12 ("C 1 2/m 1" setting).

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(10.2, 5.94, 7.81, 90, 117.7, 90)
positions = [[0, 0, 0], [1/4, 1/4, 0]]
types = ["Cu1", "Cu2"]
cryst = Crystal(latvecs, positions, "C 1 2/m 1"; types)
view_crystal(cryst)

# Define the interactions.

moments = [1 => Moment(s=1/2, g=2), 3 => Moment(s=1/2, g=2)]
sys = System(cryst, moments, :dipole)
J   = -2
Jp  = -1
Jab = 0.75
Ja  = -J/.66 - Jab
Jip = 0.01
set_exchange!(sys, J, Bond(1, 3, [0, 0, 0]))
set_exchange!(sys, Jp, Bond(3, 5, [0, 0, 0]))
set_exchange!(sys, Ja, Bond(3, 4, [0, 0, 0]))
set_exchange!(sys, Jab, Bond(1, 2, [0, 0, 0]))
set_exchange!(sys, Jip, Bond(3, 4, [0, 0, 1]))

# Use [`minimize_spiral_energy!`](@ref) to optimize the generalized spiral
# order. This determines the propagation wavevector `k`, and fits the spin
# values within the unit cell. One must provide a fixed `axis` perpendicular to
# the polarization plane. For this system, all interactions are rotationally
# invariant, and the `axis` vector is arbitrary. In other cases, a good `axis`
# will frequently be determined from symmetry considerations.

axis = [0, 0, 1]
randomize_spins!(sys)
k = minimize_spiral_energy!(sys, axis; k_guess=randn(3))
plot_spins(sys; ndims=2)

# If successful, the optimization process will find one two propagation
# wavevectors, `±k_ref`, with opposite chiralities. In this system, the
# [`spiral_energy_per_site`](@ref) is independent of chirality.

k_ref = [0.785902495, 0.0, 0.107048756]
k_ref_alt = [1, 0, 1] - k_ref
@assert isapprox(k, k_ref; atol=1e-6) || isapprox(k, k_ref_alt; atol=1e-6)
@assert spiral_energy_per_site(sys; k, axis) ≈ -0.78338383838

# Check the energy with a real-space calculation using a large magnetic cell.
# First, we must determine a lattice size for which k becomes approximately
# commensurate. 

suggest_magnetic_supercell([k_ref]; tol=1e-3)

# Resize the system as suggested, and perform a real-space calculation. Working
# with a commensurate wavevector increases the energy slightly. The precise
# value might vary from run-to-run due to trapping in a local energy minimum.

new_shape = [14 0 1; 0 1 0; 0 0 2]
sys2 = reshape_supercell(sys, new_shape)
randomize_spins!(sys2)
minimize_energy!(sys2)
energy_per_site(sys2)

# Return to the original system (with a single chemical cell) and construct
# [`SpinWaveTheorySpiral`](@ref) for calculations on the incommensurate spiral
# phase.

measure = ssf_perp(sys; apply_g=false)
swt = SpinWaveTheorySpiral(sys; measure, k, axis)

# Plot intensities for a path through ``𝐪``-space.

qs = [[0,0,0], [1,0,0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)

# Plot the powder-averaged intensities

radii = range(0, 2, 100) # (1/Å)
energies = range(0, 6, 200)
kernel = gaussian(fwhm=0.05)
res = powder_average(cryst, radii, 400) do qs
    intensities(swt, qs; energies, kernel)
end
plot_intensities(res; units)
