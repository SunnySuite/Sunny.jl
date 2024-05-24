# # SW18 - Distorted kagome
#
# This is a Sunny port of [SpinW Tutorial
# 18](https://spinw.org/tutorials/19tutorial), originally authored by Goran
# Nilsen and Sandor Toth. This tutorial illustrates Sunny's **experimental** for
# studying incommensurate, single-``Q`` structures. The test system is
# KCu₃As₂O₇(OD)₃. The Cu ions are arranged in a distorted kagome lattice, and
# exhibit helical magnetic order, as described in [G. J. Nilsen, et al., Phys.
# Rev. B **89**, 140412 (2014)](https://doi.org/10.1103/PhysRevB.89.140412).

using Sunny, GLMakie, Test, LinearAlgebra

# Build a crystal

latvecs = lattice_vectors(10.2, 5.94, 7.81, 90, 117.7, 90)
positions = [[0, 0, 0], [1/4, 1/4, 0]]
types = ["Cu1","Cu2"]
cryst = Crystal(latvecs,positions,12;types,setting="b1")
view_crystal(cryst)

# Define the interactions

sys = System(cryst,(1,1,1), [SpinInfo(1, S=1/2, g=2), SpinInfo(3, S=1/2, g=2)], :dipole, seed=0)
J   = -2
Jp  = -1
Jab = 0.75
Ja  = -J/.66 - Jab
Jip = 0.01
set_exchange!(sys,J,Bond(1, 3, [0, 0, 0]))
set_exchange!(sys,Jp,Bond(3, 5, [0, 0, 0]))
set_exchange!(sys,Ja,Bond(3, 4, [0, 0, 0]))
set_exchange!(sys,Jab, Bond(1, 2, [0,0,0]))
set_exchange!(sys,Jip,Bond(3, 4, [0, 0, 1]))

# Optimize the spiral structure

axis = [0, 0, 1]
randomize_spins!(sys)
k = Sunny.minimize_energy_spiral!(sys, axis; k_guess=randn(3))

# If convergence succeeds, there are two possible results ±k_ref. These each
# correspond to a different chirality of the magnetic order.

k_ref = [0.785902495, 0.0, 0.107048756]
@assert isapprox(k, k_ref; atol=1e-6) || isapprox(k, [1, 0, 1] - k_ref ; atol=1e-6)
@assert Sunny.spiral_energy_per_site(sys, axis, k) ≈ -0.5869788384 # 3/4 the SpinW result?


# Define a path in reciprocal space.

q_points = [[0,0,0], [1,0,0]]
density = 200
path, xticks = reciprocal_space_path(cryst, q_points, density);
swt = SpinWaveTheory(sys)
formula = Sunny.intensity_formula_SingleQ(swt,k,axis, :perp; kernel=delta_function_kernel)
disp, intensity = Sunny.intensities_bands_SingleQ(swt, path, formula);

γ = 0.01
energies = collect(0:0.01:5.5)
broadened_formula = Sunny.intensity_formula_SingleQ(swt,k,axis, :perp; kernel=lorentzian(γ), formfactors=nothing)
is = Sunny.intensities_broadened_SingleQ(swt, path, energies, broadened_formula);

fig = Figure()
ax = Axis(fig[1,1]; xlabel="Momentum (r.l.u.)", ylabel="Energy (meV)", title="Spin wave dispersion: ", xticks)
ylims!(ax, 0, 5)
xlims!(ax, 1, size(disp, 1))
for i in axes(disp, 2)
    lines!(ax, 1:length(disp[:,i]), disp[:,i];color="black",colorrange = (0,0.03))
end
fig



radii = 0.01:0.02:3 
output = zeros(Float64, length(radii), length(energies))
for (i, radius) in enumerate(radii)
    axis = 300
    qs = reciprocal_space_shell(cryst, radius, axis)
    is1 = Sunny.intensities_broadened_SingleQ(swt, qs, energies, broadened_formula);
    is2=dropdims(sum(is1[:,:,:,:],dims=[3,4]),dims=(3,4))
    output[i, :] = sum(is2, dims=1) / size(is2, 1)
end

fig = Figure()
ax = Axis(fig[1,1]; xlabel="Q (Å⁻¹)", ylabel="ω (meV)",title="Convoluted powder spectra:")
ylims!(ax, 0.0, 5)
heatmap!(ax, radii, energies, output, colormap=:gnuplot2,colorrange=(0.0,1))
fig
