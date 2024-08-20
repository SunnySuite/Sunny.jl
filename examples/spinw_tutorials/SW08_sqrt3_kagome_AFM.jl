# # SW08 - √3×√3 kagome antiferromagnet
#
# This is a Sunny port of [SpinW Tutorial
# 8](https://spinw.org/tutorials/08tutorial), originally authored by Bjorn Fak
# and Sandor Toth. It calculates the linear spin wave theory spectrum for the
# ``\sqrt{3} \times \sqrt{3}`` order of a Kagome antiferromagnet.

# Load Sunny and the GLMakie plotting package

using Sunny, GLMakie

# Define the chemical cell of a kagome lattice with spacegroup 147 (P-3).

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(6, 6, 40, 90, 90, 120)
cryst = Crystal(latvecs, [[1/2, 0, 0]], 147)
view_crystal(cryst; dims=2)

# Construct a spin system with nearest neighbor antiferromagnetic exchange.

sys = System(cryst, (3, 3, 1), [SpinInfo(1; S=1, g=2)], :dipole)
J = 1.0
set_exchange!(sys, J, Bond(2, 3, [0, 0, 0]))

# Initialize to an energy minimizing magnetic structure, for which
# nearest-neighbor spins are at 120° angles.

k = -[1/3, 1/3, 0]
axis = [0,0,1]
set_spiral_order_on_sublattice!(sys, 1; k, axis, S0=[cos(0),sin(0),0])
set_spiral_order_on_sublattice!(sys, 2; k, axis, S0=[cos(0),sin(0),0])
set_spiral_order_on_sublattice!(sys, 3; k, axis, S0=[cos(2π/3),sin(2π/3),0])
plot_spins(sys; dims=2)

# Check energy. Each site participates in 4 bonds with energy ``J\cos(2π/3)``.
# Factor of 1/2 avoids double counting.

@assert energy_per_site(sys) ≈ (4/2)*J*cos(2π/3)

# Calculate and plot intensities for a path through ``𝐪``-space. Note the very
# intense flat band at zero energy transfer.

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
q_points = [[-1/2,0,0], [0,0,0], [1/2,1/2,0]]
path = q_space_path(cryst, q_points, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)

# Calculate and plot the powder averaged spectrum. (This calculation is disabled
# from the automatic builds because it takes about two minutes to run.)

# ```julia
# # radii = range(0, 2.5, 200)
# # energies = range(0, 3, 200)
# # kernel = gaussian(fwhm=0.05)
# # res = powder_average(cryst, radii, 1000) do qs
# #     intensities(swt, qs; energies, kernel)
# # end
# # plot_intensities(res; units, saturation=0.8)
# ```
