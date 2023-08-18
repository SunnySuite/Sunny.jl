# # SpinW parity

using Sunny, GLMakie

# ## SpinW Tutorial 1

latvecs = lattice_vectors(3,8,8,90,90,90)
cryst = Crystal(latvecs,[[0,0,0]]; types = ["Cu"])
sys = System(cryst, (1,1,1), [SpinInfo(1,S=1,g=2)], :SUN, seed = 0)
set_exchange!(sys,-1,Bond(1,1,(1,0,0)))

@assert energy(sys) == -1.0

swt = SpinWaveTheory(sys)
formula = intensity_formula(swt,:perp;kernel = delta_function_kernel,mode_fast = true)
qxs = range(0,1,length=100)
qs = [[x,0,0] for x in qxs]
dispersion, intensities = intensities_bands(swt,qs)

f = Figure()
ax = Axis(f[1,1],title = "Spin wave dispersion: ω(Q), T = 0.0K", xlabel = "(qx,0,0)", ylabel = "Energy [meV]")
colorrange = (minimum(intensities),maximum(intensities))
for band_ix in axes(dispersion)[2]
    lines!(ax, qxs, dispersion[:,band_ix]; color=intensities[:,band_ix], linewidth = 2.0, colorrange)
end
Colorbar(f[1,2],ax.scene.plots[2])
display(f)

#@assert all(intensities[:,2] .== 0.5)


r_params = (0,2.5,2.5/100)
ω_params = (0,4.5,4.5/250)
is, counts = Sunny.powder_average_bin_centers(swt, r_params, ω_params, 1000, formula)
rbcs = axes_bincenters(r_params...)[1]
wbcs = axes_bincenters(ω_params...)[1]
heatmap(rbcs,wbcs,is ./ counts)

