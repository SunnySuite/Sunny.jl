# # SW04 - Frustrated square lattice
#
# This is a Sunny port of [SpinW Tutorial
# 4](https://spinw.org/tutorials/04tutorial), originally authored by Bjorn Fak
# and Sandor Toth. It calculates the spin wave spectrum of the frustrated square
# lattice.

# Load Sunny and the GLMakie plotting package

using Sunny, GLMakie

# To model the 2D square lattice, create an elongated tetragonal cell with one
# atom.

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(3.0, 3.0, 6.0, 90, 90, 90) 
cryst = Crystal(latvecs, [[0, 0, 0]])
view_crystal(cryst; ndims=2)

# Construct a spin system with competing nearest-neighbor (AFM) and
# next-nearest-neighbor (FM) interactions. The Néel magnetic order requires a
# supercell of 2×2 chemical cells.

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2, 2, 1))
J1 = 1.0
J2 = -0.1
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J2, Bond(1, 1, [1, 1, 0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)

# Calculate and plot intensities for a path through ``𝐪``-space.

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)
