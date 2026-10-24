using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.10.0"

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(4.1457, 4.1457, 13.801, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]])
sys = System(cryst, [1 => Moment(s=5/2, g=2)], :dipole; dims=(2, 2, 1), seed=0)
J = 0.72
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
set_onsite_coupling!(sys, S -> -0.008*S[3]^2, 1)
randomize_spins!(sys)
minimize_energy!(sys)

function diluted_system(sys, x)
    sys_inhom = to_inhomogeneous(resize_supercell(sys, (80, 80, 1)))
    for site in eachsite(sys_inhom)
        if rand(sys_inhom.rng) < x
            set_vacancy_at!(sys_inhom, site)
        end
    end
    minimize_energy!(sys_inhom)
    return sys_inhom
end

qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path = q_space_path(cryst, qs, 150; labels=["Γ", "X", "M", "Γ"])
energies = range(0, 8, 200)
kernel = gaussian(fwhm=0.2)

fig = Figure(size=(768, 600))
for (i, x) in enumerate([0.0, 0.1, 0.2, 0.46])
    sys_inhom = diluted_system(sys, x)
    swt = SpinWaveTheoryLarge(sys_inhom; measure=ssf_perp(sys_inhom))
    res = intensities(swt, path; energies, kernel, threaded=true)
    plot_intensities!(fig[fldmod1(i, 2)...], res; units, colorrange=(0, 8), title="x = $x")
end
fig

sys_inhom = diluted_system(sys, 0.46)
formfactors = [1 => FormFactor("Mn2")]
swt = SpinWaveTheoryLarge(sys_inhom; measure=ssf_perp(sys_inhom; formfactors))
hs = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
energies = range(1, 8, 141)
kernel = gaussian(fwhm=0.6)
res = intensities(swt, [[h, h, 2.9] for h in hs]; energies, kernel, kT=4.0*units.K, threaded=true)

fig = Figure(size=(768, 650))
for (i, h) in enumerate(hs)
    ax = Axis(fig[fldmod1(i, 2)...]; xlabel="Energy (meV)", ylabel="Intensity", title="𝐪 = ($h, $h, 2.9)")
    lines!(ax, energies, res.data[:, i])
    h == 0.25 && vlines!(ax, [3.15, 4.85, 6.40]; color=:gray, linestyle=:dash, linewidth=3, label="Measured peaks")
    h == 0.25 && axislegend(ax; position=:lb)
end
fig
