# # 9. Spin waves in disordered magnets
#
# This example uses [`to_inhomogeneous`](@ref) and [`SpinWaveTheoryLarge`](@ref)
# to calculate the spectrum of Rb₂Mn``_c``Mg``_{1-c}``F₄, a disordered quasi-2D
# antiferromagnet. A Mn concentration ``c < 1`` introduces random magnetic
# vacancies. Inelastic neutron measurements were reported in [Cowley et al.,
# Phys. Rev. B **15**, 4292 (1977)](https://doi.org/10.1103/PhysRevB.15.4292).
# Sunny models this data using spin wave calculations on very large supercells.

using Sunny, GLMakie

# Model a single Mn layer using the lattice constants of the diluted crystal.
# Following Cowley et al., the model retains just the effective nearest-neighbor
# exchange ``J = 0.72`` meV. A weak easy-axis anisotropy orients the spins
# normal to the layer, and reproduces the zone-center gap of pure Rb₂MnF₄.

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(4.1457, 4.1457, 13.801, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])
sys = System(cryst, [1 => Moment(s=5/2, g=2)], :dipole; dims=(2, 2, 1), seed=0)
J = 0.72
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
set_onsite_coupling!(sys, S -> -0.008*S[3]^2, 1)
randomize_spins!(sys)
minimize_energy!(sys)

# A helper function that builds an 80×80 supercell, converts it
# [`to_inhomogeneous`](@ref), and calls [`set_vacancy_at!`](@ref) for a fraction
# ``1 - c`` of random sites. These vacancies model substitution of non-magnetic
# Mg. A more complicated model could also modify the local couplings via
# [`set_exchange_at!`](@ref), [`set_onsite_coupling_at!`](@ref), and related
# functions.

function diluted_system(sys, c)
    sys_inhom = to_inhomogeneous(resize_supercell(sys, (80, 80, 1)))
    for site in eachsite(sys_inhom)
        if rand(sys_inhom.rng) < 1 - c
            set_vacancy_at!(sys_inhom, site)
        end
    end
    minimize_energy!(sys_inhom)
    return sys_inhom
end

# Calculate ``S(𝐪, ω)`` at various Mn concentrations ``c``.
# [`SpinWaveTheoryLarge`](@ref) is optimized for large supercells such as this
# one. It uses the Lanczos method instead of direct matrix diagonalization to
# reduce computational cost to linear scaling in system size.
#
# Undiluted, there is a single sharp magnon. With dilution, the long-wavelength
# magnon near the ordering wavevector M stays relatively coherent, while modes
# near the zone boundary become strongly damped. Simultaneously, flat
# ``𝐪``-independent bands appear. These are localized excitations of spins with
# different numbers of magnetic neighbors. Near and below the percolation
# threshold ``c = 0.592…``, the low-energy weight is dominated by finite
# clusters.

qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 150; labels=["Γ", "X", "M", "Γ"])
energies = range(0, 8, 200)
kernel = gaussian(fwhm=0.3)

fig = Figure(size=(768, 600))
for (i, c) in enumerate([1.0, 0.9, 0.8, 0.54])
    sys_inhom = diluted_system(sys, c)
    swt = SpinWaveTheoryLarge(sys_inhom; measure=ssf_perp(sys_inhom))
    res = intensities(swt, path; energies, kernel, threaded=true, verbose=true)
    plot_intensities!(fig[fldmod1(i, 2)...], res; units, colorrange=(0, 8), title="c = $c")
end
fig

# Fig. 5 of Cowley et al. reports neutron intensity at Mn concentration ``c =
# 0.54`` and temperature 4.0 K. Their RLU notation uses a √2×√2 magnetic cell
# whereby ``(ξ, 0, l)`` maps to Sunny's ``(ξ/2, ξ/2, l)``. At the magnetic zone
# boundary ``ξ = 1/2``, the three peaks represent Mn²⁺ ions with 2, 3, and 4
# magnetic neighbors. Toward the magnetic zone center ``ξ = 1`` the fine
# structure fades, and the intensity grows and moves to lower energies.

sys_inhom = diluted_system(sys, 0.54)
formfactors = [1 => FormFactor("Mn2")]
swt = SpinWaveTheoryLarge(sys_inhom; measure=ssf_perp(sys_inhom; formfactors))
ξs = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
energies = range(1, 8, 141)
kernel = gaussian(fwhm=0.6)
res = intensities(swt, [[ξ/2, ξ/2, 2.9] for ξ in ξs]; energies, kernel, kT=4.0*units.K, threaded=true)

# Experimental intensities digitized from Fig. 5 of Cowley et al. for each ξ. As
# in the paper, a single overall scale factor matches the integrated intensity
# between 1 and 7 meV at the zone boundary ξ = 1/2.

exp_energies = 1.0:0.2:8.0
exp_data = [
    [2.0, 1.77, 1.39, 1.9, 1.7, 1.52, 1.43, 1.62, 1.61, 1.8, 1.98, 1.99, 1.93, 1.44, 1.35, 1.19, 1.56, 1.6, 2.02, 2.1, 1.95, 1.93, 1.08, 0.67, 0.64, 0.68, 1.06, 0.99, 0.78, 0.7, 0.41, 0.21, 0.13, -0.1, 0.06, 0.03],
    [2.87, 2.8, 2.61, 2.46, 2.1, 2.89, 2.54, 2.26, 2.84, 3.0, 3.07, 2.6, 2.54, 2.49, 2.03, 2.38, 1.88, 2.08, 2.73, 2.45, 2.32, 1.61, 1.26, 0.96, 0.8, 1.14, 1.04, 1.31, 1.04, 0.73, 0.43, 0.36, 0.08, -0.01, -0.03, 0.09],
    [4.61, 3.86, 3.41, 3.92, 3.53, 3.21, 3.93, 3.94, 3.78, 3.84, 3.54, 3.04, 3.03, 2.65, 2.57, 2.39, 2.55, 2.67, 2.47, 1.93, 1.54, 1.13, 0.97, 0.67, 0.84, 0.7, 0.81, 0.62, 0.8, 0.46, 0.28, 0.04, 0.01, -0.19, -0.11, -0.08],
    [10.68, 9.58, 8.34, 7.64, 7.61, 8.0, 7.32, 6.65, 6.44, 3.98, 3.37, 2.83, 3.0, 3.06, 1.98, 2.33, 2.77, 2.59, 1.38, 1.3, 1.13, 0.86, 0.34, 0.62, 0.38, 0.25, 0.66, 0.5, 0.23, 0.55, 0.05, -0.01, 0.02, -0.04, -0.05, 0.05],
    [19.28, 15.79, 13.03, 10.23, 8.72, 7.91, 6.51, 5.37, 3.62, 3.54, 3.32, 3.02, 2.48, 2.77, 2.25, 2.45, 1.82, 1.71, 1.71, 1.03, 0.85, 0.75, 0.79, 1.05, 0.86, 0.8, 0.79, 0.81, 0.87, 0.79, 0.59, 0.57, 0.54, 0.54, 0.44],
    [23.06, 14.46, 10.88, 9.01, 7.95, 6.19, 5.57, 3.86, 3.19, 2.88, 2.98, 2.67, 2.76, 2.89, 2.37, 2.23, 2.03, 1.72, 1.51, 1.46, 1.18, 1.26, 1.06, 1.16, 0.95, 1.13, 1.11, 1.12, 1.05, 0.93, 0.92, 0.84, 1.01, 0.77, 0.84, 0.99],
]
scale = sum(res.data[1:4:121, 1]) / sum(exp_data[1][1:31])

fig = Figure(size=(768, 650))
for (i, ξ) in enumerate(ξs)
    ax = Axis(fig[fldmod1(i, 2)...]; xlabel="Energy (meV)", ylabel="Intensity", title="ξ = $ξ")
    lines!(ax, energies, res.data[:, i], label="Spin wave theory")
    scatter!(ax, exp_energies[eachindex(exp_data[i])], scale * exp_data[i]; marker=:xcross, color=:black, label="Cowley et al.")
    ξ == 0.5 && axislegend(ax; position=:lb)
end
fig
