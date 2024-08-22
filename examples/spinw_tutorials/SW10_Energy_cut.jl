# # SW10 - Energy cut on square lattice
#
# This is a Sunny port of [SpinW Tutorial
# 10](https://spinw.org/tutorials/10tutorial), originally authored by Sandor
# Toth. It calculates the spin wave spectrum on a constant energy cut of the
# frustrated square lattice.

# Load Sunny and the GLMakie plotting package

using Sunny, GLMakie

# Define the chemical cell for the 2D square lattice.

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(1.0, 1.0, 3.0, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])

# Construct a spin system with nearest-neighbor antiferomagnetic interactions of
# 1.0 meV. Energy minimization yields the expected Néel order.

sys = System(cryst, (2,2,1), [SpinInfo(1, S=1, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; dims=2)

# Define a 2D slice through ``𝐪``-space with [`q_space_grid`](@ref). 

nQ = 201
grid = q_space_grid(cryst, [1, 0, 0], range(0, 2, nQ), [0, 1, 0], range(0, 2, nQ))

# Calculate and plot a constant energy cut at the precise value of 3.75 meV.
# Apply a line broadening with a full-width half-max of 0.2 meV to approximately
# capture intensities between 3.5 and 4.0 meV.

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, grid; energies=[3.75], kernel=gaussian(fwhm=0.2))
plot_intensities(res; units)

# Sunny does not include Horace integration. For now, one can write custom code
# to manually integrate intensity bands between 3.5 and 4 meV.

res = intensities_bands(swt, grid)
summed_intensity = zeros(size(res.disp)[2:3])
for ci in CartesianIndices(res.disp)
    if 3.5 < res.disp[ci] < 4.0
        summed_intensity[ci[2], ci[3]] += res.data[ci]
    end
end

# To make the plot look nice, for now we use a private function internal to
# Sunny.

plot_intensities(Sunny.InstantIntensities(cryst, grid, summed_intensity))
