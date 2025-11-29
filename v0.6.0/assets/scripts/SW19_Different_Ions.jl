using Sunny, GLMakie

a = 3.0
b = 8.0
c = 4.0
latvecs = lattice_vectors(a, b, c, 90, 90, 90)
positions = [[0, 0, 0], [0, 1/2, 0]]
types = ["Cu2", "Fe2"]
cryst = Crystal(latvecs, positions, 1; types)
view_crystal(cryst)

J_Cu_Cu = 1.0
J_Fe_Fe = 1.0
J_Cu_Fe = -0.1
sys = System(cryst, (2,1,1), [SpinInfo(1,S=1/2,g=2), SpinInfo(2,S=2,g=2)], :dipole; seed=0)
set_exchange!(sys, J_Cu_Cu, Bond(1, 1, [-1, 0, 0]))
set_exchange!(sys, J_Fe_Fe, Bond(2, 2, [-1, 0, 0]))
set_exchange!(sys, J_Cu_Fe, Bond(2, 1, [0, 1, 0]))
set_exchange!(sys, J_Cu_Fe, Bond(1, 2, [0, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)

swt = SpinWaveTheory(sys)
q_points = [[0,0,0], [1,0,0]]
path, xticks = reciprocal_space_path(cryst, q_points, 200)
formula = intensity_formula(swt, :perp; kernel=delta_function_kernel)
disp, intensity = intensities_bands(swt, path, formula)
fig = Figure()
ax = Axis(fig[1,1]; xlabel="Momentum (RLU)", ylabel="Energy (meV)",
          title="Spin wave dispersion", xticks)
ylims!(ax, 0.0, 4.5)
xlims!(ax, 1, size(disp, 1))
lines!(ax, disp[:, 1]; color=vec(sum(intensity[:, 1:2]; dims=2)), colorrange=(0,2))
lines!(ax, disp[:, 3]; color=vec(sum(intensity[:, 3:4]; dims=2)), colorrange=(0,2))
fig

function plot_intensities(formfactors, title)
    formula = intensity_formula(swt, :trace; kernel=gaussian(fwhm=0.2), formfactors)
    energies = collect(0:0.02:140)
    is = intensities_broadened(swt, path, energies, formula)
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="Momentum (RLU)", ylabel="Energy (meV)", title, xticks)
    ylims!(ax, 0.0, 4.5)
    heatmap!(ax, 1:size(is, 1), energies, is; colorrange=(0.1,50),
             colormap=Reverse(:thermal), lowclip=:white)
    for d in eachcol(disp)
        lines!(ax, d; color=:pink)
    end
    fig
end

plot_intensities([FormFactor("Cu2"), FormFactor("Fe2")], "All correlations")

plot_intensities([FormFactor("Cu2"), zero(FormFactor)], "Cu-Cu correlations")

plot_intensities([zero(FormFactor), FormFactor("Fe2")], "Fe-Fe correlations")
