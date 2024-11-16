# # 9. Disordered system with KPM
#
# This example uses the [kernel polynomial method](@ref SpinWaveTheoryKPM) to
# efficiently calculate the neutron scattering spectrum of a disordered
# triangular antiferromagnet. The model is inspired by YbMgGaO4, as studied in
# [Paddison et al, Nature Phys., **13**, 117–122
# (2017)](https://doi.org/10.1038/nphys3971) and [Zhu et al, Phys. Rev. Lett.
# **119**, 157201 (2017)](https://doi.org/10.1103/PhysRevLett.119.157201).
# Disordered occupancy of non-magnetic Mg/Ga sites can be modeled as a
# stochastic distribution of exchange constants and ``g``-factors. Including
# this disorder introduces broadening of the spin wave spectrum.

using Sunny, GLMakie

# Set up minimal triangular lattice system. Include antiferromagnetic exchange
# interactions between nearest neighbor bonds. Energy minimization yields the
# magnetic ground state with 120° angles between spins in triangular plaquettes.

latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
cryst = Crystal(latvecs, [[0, 0, 0]])
sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1))
set_exchange!(sys, +1.0, Bond(1, 1, [1,0,0]))

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles], ndims=2)

# Select a ``𝐪``-space path for the spin wave calculations.

qs = [[0, 0, 0], [1/3, 1/3, 0], [1/2, 0, 0], [0, 0, 0]]
labels = ["Γ", "K", "M", "Γ"]
path = q_space_path(cryst, qs, 150; labels)

# Perform a traditional spin wave calculation. The spectrum shows sharp modes
# associated with coherent excitations about the K-point ordering wavevector,
# ``𝐪 = [1/3, 1/3, 0]``.

kernel = lorentzian(fwhm=0.4)
energies = range(0.0, 3.0, 150)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

# Use [`repeat_periodically`](@ref) to enlarge the system by a factor of 10 in
# each dimension. Use [`to_inhomogeneous`](@ref) to disable symmetry
# constraints, and allow for the addition of disordered interactions.

sys_inhom = to_inhomogeneous(repeat_periodically(sys, (10, 10, 1)))

# Use [`symmetry_equivalent_bonds`](@ref) to iterate over all nearest neighbor
# bonds of the inhomogeneous system. Modify each AFM exchange with a noise term
# that has variance of 1/3. The newly minimized energy configuration allows for
# long wavelength modulations on top of the original 120° order.

for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1,1,[1,0,0]))
    noise = randn()/3
    set_exchange_at!(sys_inhom, 1.0 + noise, site1, site2; offset)
end

minimize_energy!(sys_inhom, maxiters=5_000)
plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], ndims=2)

# Traditional spin wave theory calculations become impractical for large system
# sizes. Significant acceleration is possible with the [kernel polynomial
# method](https://arxiv.org/abs/2312.08349). Enable it by selecting
# [`SpinWaveTheoryKPM`](@ref) in place of the traditional
# [`SpinWaveTheory`](@ref). Using KPM, the cost of an [`intensities`](@ref)
# calculation becomes linear in system size and scales inversely with the width
# of the line broadening `kernel`. Error tolerance is controlled through the
# dimensionless `tol` parameter. A relatively small value, `tol = 0.01`, helps
# to resolve the large intensities near the ordering wavevector. The alternative
# choice `tol = 0.1` would be twice faster, but would introduce significant
# numerical artifacts.
#
# Observe from the KPM calculation that disorder in the nearest-neighbor
# exchange serves to broaden the discrete excitation bands into a continuum.

swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.01)
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

# Now apply a magnetic field of magnitude 7.5 (energy units) along the global
# ``ẑ`` axis. This field fully polarizes the spins. Because gap opens, a larger
# tolerance of `tol = 0.1` can be used to accelerate the KPM calculation without
# sacrificing much accuracy. The resulting spin wave spectrum shows a sharp mode
# at the Γ-point (zone center) that broadens into a continuum along the K and M
# points (zone boundary).

set_field!(sys_inhom, [0, 0, 7.5])
randomize_spins!(sys_inhom)
minimize_energy!(sys_inhom)

energies = range(0.0, 9.0, 150)
swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.1)
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

# Add disorder to the ``z``-component of each magnetic moment ``g``-tensor. This
# further broadens intensities, now across the entire path. Some intensity
# modulation within the continuum is also apparent. This modulation is a
# finite-size effect, and would be mitigated by enlarging the system beyond
# 30×30 chemical cells.

for site in eachsite(sys_inhom)
    noise = randn()/6
    sys_inhom.gs[site] = [1 0 0; 0 1 0; 0 0 1+noise]
end

swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.1)
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

# For reference, the equivalent non-disordered system shows a single coherent
# mode.

set_field!(sys, [0, 0, 7.5])
randomize_spins!(sys)
minimize_energy!(sys)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, path; energies, kernel)
plot_intensities(res)
