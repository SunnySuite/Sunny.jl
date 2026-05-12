using Sunny, GLMakie

a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], 227, setting="1")

view_crystal(cryst)

latsize = (1, 1, 1)
S = 3/2
J = 7.5413*meV_per_K # (~ 0.65 meV)
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=0)
set_exchange!(sys, J, Bond(1, 3, [0,0,0]))

randomize_spins!(sys)
minimize_energy!(sys)

@assert energy_per_site(sys) ≈ -2J*S^2

s0 = sys.dipoles[1,1,1,1]
plot_spins(sys; color=[s'*s0 for s in sys.dipoles])

shape = [0 1 1;
         1 0 1;
         1 1 0] / 2
sys_prim = reshape_supercell(sys, shape)
@assert energy_per_site(sys_prim) ≈ -2J*S^2
plot_spins(sys_prim; color=[s'*s0 for s in sys_prim.dipoles])

swt = SpinWaveTheory(sys_prim)
kernel = lorentzian(fwhm=0.8)
formfactors = [FormFactor("Co2")]
formula = intensity_formula(swt, :perp; kernel, formfactors)

qpoints = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]
path, xticks = reciprocal_space_path(cryst, qpoints, 50)
energies = collect(0:0.01:6)
is = intensities_broadened(swt, path, energies, formula)

fig = Figure()
ax = Axis(fig[1,1]; aspect=1.4, ylabel="ω (meV)", xlabel="𝐪 (r.l.u.)",
          xticks, xticklabelrotation=π/10)
heatmap!(ax, 1:size(is, 1), energies, is, colormap=:gnuplot2)
fig

radii = 0.01:0.02:3 # (1/Å)
output = zeros(Float64, length(radii), length(energies))
for (i, radius) in enumerate(radii)
    n = 300
    qs = reciprocal_space_shell(cryst, radius, n)
    is = intensities_broadened(swt, qs, energies, formula)
    output[i, :] = sum(is, dims=1) / size(is, 1)
end

fig = Figure()
ax = Axis(fig[1,1]; xlabel="Q (Å⁻¹)", ylabel="ω (meV)")
heatmap!(ax, radii, energies, output, colormap=:gnuplot2)
fig
