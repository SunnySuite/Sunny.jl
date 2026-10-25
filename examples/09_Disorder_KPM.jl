# # 9. Disordered systems
#
# This example uses [`SpinWaveTheoryLarge`](@ref) to efficiently calculate the
# neutron scattering spectrum of a disordered triangular antiferromagnet. The
# model is inspired by YbMgGaO4, as studied in [Paddison et al, Nat. Phys.,
# **13**, 117–122 (2017)](https://doi.org/10.1038/nphys3971) and [Zhu et al,
# Phys. Rev. Lett. **119**, 157201
# (2017)](https://doi.org/10.1103/PhysRevLett.119.157201). Disordered occupancy
# of non-magnetic Mg/Ga sites can be modeled as a stochastic distribution of
# exchange constants and ``g``-factors. Including this disorder introduces
# broadening of the spin wave spectrum.

using Sunny, GLMakie, LinearAlgebra

# Set up minimal triangular lattice system. Include antiferromagnetic exchange
# interactions between nearest neighbor bonds, with a slight easy-plane
# anisotropy, as appropriate for YbMgGaO4. An optional `seed` for random number
# generation makes the calculations exactly reproducible on a fixed Julia
# version.

latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
cryst = Crystal(latvecs, [[0, 0, 0]])
sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
J = diagm([1, 1, 0.95])
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))

# Energy minimization yields the magnetic ground state with 120° angles between
# spins in triangular plaquettes, lying in the easy plane.

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles], ndims=2)

# Select a ``𝐪``-space path for the spin wave calculations.

qs = [[0, 0, 0], [1/3, 1/3, 0], [1/2, 0, 0], [0, 0, 0]]
labels = ["Γ", "K", "M", "Γ"]
path = q_space_path(cryst, qs, 150; labels)

# Perform a traditional spin wave calculation. The spectrum shows sharp modes
# associated with coherent excitations about the K-point ordering wavevector,
# ``𝐪 = [1/3, 1/3, 0]``. The Gaussian `kernel` models instrumental resolution.
# A small intrinsic broadening `η` gives each mode a Lorentzian line shape, and
# damps the divergent intensity of the Goldstone mode at the K point.

η = 0.05
kernel = gaussian(fwhm=0.4)
energies = range(0.0, 3.0, 150)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, path; energies, η, kernel)
plot_intensities(res)

# Use [`repeat_periodically`](@ref) to enlarge the system by a factor of 10 in
# each dimension, for a total of 30×30 chemical cells. Use
# [`to_inhomogeneous`](@ref) to disable symmetry constraints and allow for the
# addition of disordered interactions.

sys_inhom = to_inhomogeneous(repeat_periodically(sys, (10, 10, 1)))

# Iterate over all [`symmetry_equivalent_bonds`](@ref) for nearest neighbor
# sites in the inhomogeneous system. Scale the AFM exchange by a noise factor
# that has a standard deviation of 1/3. The newly minimized energy configuration
# allows for long wavelength modulations on top of the original 120° order. The
# easy-plane anisotropy makes the minimization fast and the result a stable
# energy minimum.

for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1, 1, [1, 0, 0]))
    noise = randn(sys_inhom.rng)/3
    set_exchange_at!(sys_inhom, (1.0 + noise) * J, site1, site2; offset)
end

minimize_energy!(sys_inhom)
plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], ndims=2)

# Traditional spin wave theory with exact diagonalization becomes impractical
# for large system sizes. For such systems, use instead
# [`SpinWaveTheoryLarge`](@ref), which employs a highly scalable algorithm.
#
# Observe that disorder in the nearest-neighbor exchange serves to broaden the
# discrete excitation bands into a continuum.

swt = SpinWaveTheoryLarge(sys_inhom; measure=ssf_perp(sys_inhom))
res = intensities(swt, path; energies, η, kernel, threaded=true)
plot_intensities(res)

# Now apply a magnetic field of magnitude 7.5 (energy units) along the global
# ``ẑ`` axis. This field fully polarizes the spins and opens a gap. The spin
# wave spectrum shows a sharp mode at the Γ-point (zone center) that broadens
# into a continuum along the K and M points (zone boundary).

set_field!(sys_inhom, [0, 0, 7.5])
randomize_spins!(sys_inhom)
minimize_energy!(sys_inhom)

energies = range(0.0, 9.0, 150)
swt = SpinWaveTheoryLarge(sys_inhom; measure=ssf_perp(sys_inhom))
res = intensities(swt, path; energies, η, kernel, threaded=true)
plot_intensities(res)

# Add disorder to the ``z``-component of each magnetic moment ``g``-tensor. This
# further broadens intensities, now across the entire path. Some intensity
# modulation within the continuum is also apparent. This modulation is a
# finite-size effect and would be mitigated by further enlarging the system.

for site in eachsite(sys_inhom)
    noise = randn(sys_inhom.rng)/8
    sys_inhom.gs[site] = [1 0 0; 0 1 0; 0 0 1+noise]
end
randomize_spins!(sys_inhom)
minimize_energy!(sys_inhom)

swt = SpinWaveTheoryLarge(sys_inhom; measure=ssf_perp(sys_inhom))
res = intensities(swt, path; energies, η, kernel, threaded=true)
plot_intensities(res)

# For reference, the equivalent non-disordered system shows a single coherent
# mode.

set_field!(sys, [0, 0, 7.5])
randomize_spins!(sys)
minimize_energy!(sys)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, path; energies, η, kernel)
plot_intensities(res)
