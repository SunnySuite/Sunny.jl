using Sunny, LinearAlgebra, GLMakie

a = b = 4.05012#hide
c = 6.75214#hide
latvecs = lattice_vectors(a, b, c, 90, 90, 120)#hide
positions = [[0,0,0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]#hide
types = ["Fe", "I", "I"]#hide
FeI2 = Crystal(latvecs, positions; types)#hide
cryst = subcrystal(FeI2, "Fe")#hide
sys = System(cryst, (4,4,4), [SpinInfo(1,S=1,g=2)], :SUN, seed=2)#hide
J1pm   = -0.236#hide
J1pmpm = -0.161#hide
J1zpm  = -0.261#hide
J2pm   = 0.026#hide
J3pm   = 0.166#hide
J′0pm  = 0.037#hide
J′1pm  = 0.013#hide
J′2apm = 0.068#hide
J1zz   = -0.236#hide
J2zz   = 0.113#hide
J3zz   = 0.211#hide
J′0zz  = -0.036#hide
J′1zz  = 0.051#hide
J′2azz = 0.073#hide
J1xx = J1pm + J1pmpm#hide
J1yy = J1pm - J1pmpm#hide
J1yz = J1zpm#hide
set_exchange!(sys, [J1xx 0.0 0.0; 0.0 J1yy J1yz; 0.0 J1yz J1zz], Bond(1,1,[1,0,0]))#hide
set_exchange!(sys, [J2pm 0.0 0.0; 0.0 J2pm 0.0; 0.0 0.0 J2zz], Bond(1,1,[1,2,0]))#hide
set_exchange!(sys, [J3pm 0.0 0.0; 0.0 J3pm 0.0; 0.0 0.0 J3zz], Bond(1,1,[2,0,0]))#hide
set_exchange!(sys, [J′0pm 0.0 0.0; 0.0 J′0pm 0.0; 0.0 0.0 J′0zz], Bond(1,1,[0,0,1]))#hide
set_exchange!(sys, [J′1pm 0.0 0.0; 0.0 J′1pm 0.0; 0.0 0.0 J′1zz], Bond(1,1,[1,0,1]))#hide
set_exchange!(sys, [J′2apm 0.0 0.0; 0.0 J′2apm 0.0; 0.0 0.0 J′2azz], Bond(1,1,[1,2,1]))#hide
D = 2.165#hide
set_onsite_coupling!(sys, S -> -D*S[3]^2, 1)#hide
sys

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[s[3] for s in sys.dipoles])

kT = 0.2  # Temperature in meV
langevin = Langevin(; damping=0.2, kT)

suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.027;

for _ in 1:10_000
    step!(sys, langevin)
end

suggest_timestep(sys, langevin; tol=1e-2)

plot_spins(sys; color=[s[3] for s in sys.dipoles])

sys_large = resize_supercell(sys, (16,16,4)) # 16x16x4 copies of the original unit cell
plot_spins(sys_large; color=[s[3] for s in sys_large.dipoles])

kT = 3.5 * meV_per_K     # 3.5K ≈ 0.30 meV
langevin.kT = kT
for _ in 1:10_000
    step!(sys_large, langevin)
end

suggest_timestep(sys_large, langevin; tol=1e-2)
langevin.dt = 0.040;

dt = 2*langevin.dt
ωmax = 7.5  # Maximum energy to resolve (meV)
nω = 120    # Number of energies to resolve
sc = dynamical_correlations(sys_large; dt, nω, ωmax)

add_sample!(sc, sys_large)

for _ in 1:2
    for _ in 1:1000               # Enough steps to decorrelate spins
        step!(sys_large, langevin)
    end
    add_sample!(sc, sys_large)
end

sc

formula = intensity_formula(sc, :trace; kT)

qs = [[0, 0, 0], [0.5, 0.5, 0.5]]
is = intensities_interpolated(sc, qs, formula; interpolation = :round)

ωs = available_energies(sc)
fig = lines(ωs, is[1,:]; axis=(xlabel="meV", ylabel="Intensity"), label="(0,0,0)")
lines!(ωs, is[2,:]; label="(π,π,π)")
axislegend()
fig

formfactors = [FormFactor("Fe2"; g_lande=3/2)]
new_formula = intensity_formula(sc, :perp; kT, formfactors = formfactors)

points = [[0,   0, 0],  # List of wave vectors that define a path
          [1,   0, 0],
          [0,   1, 0],
          [1/2, 0, 0],
          [0,   1, 0],
          [0,   0, 0]]
density = 40
path, xticks = reciprocal_space_path(cryst, points, density);

is_interpolated = intensities_interpolated(sc, path, new_formula;
    interpolation = :linear,       # Interpolate between available wave vectors
);
# Add artificial broadening
is_interpolated_broadened = broaden_energy(sc, is, (ω, ω₀) -> lorentzian(fwhm=0.1)(ω-ω₀));

cut_width = 0.3
density = 15
paramsList, markers, ranges = reciprocal_space_path_bins(sc,points,density,cut_width);

total_bins = ranges[end][end]
energy_bins = paramsList[1].numbins[4]
is_binned = zeros(Float64,total_bins,energy_bins)
integrated_kernel = integrated_lorentzian(fwhm=0.1) # Lorentzian broadening
for k in eachindex(paramsList)
    bin_data, counts = intensities_binned(sc,paramsList[k], new_formula;
        integrated_kernel = integrated_kernel
    )
    is_binned[ranges[k],:] = bin_data[:,1,1,:] ./ counts[:,1,1,:]
end

fig = Figure()
ax_top = Axis(fig[1,1],ylabel = "meV",xticklabelrotation=π/8,xticklabelsize=12;xticks)
ax_bottom = Axis(fig[2,1],ylabel = "meV",xticks = (markers, string.(points)),xticklabelrotation=π/8,xticklabelsize=12)

heatmap!(ax_top,1:size(is_interpolated,1), ωs, is_interpolated;
    colorrange=(0.0,0.07),
)

heatmap!(ax_bottom,1:size(is_binned,1), ωs, is_binned;
    colorrange=(0.0,0.05),
)

fig

ωidx = 60
target_ω = ωs[ωidx]

params = unit_resolution_binning_parameters(sc)
params.binstart[1:2] .= -1 # Expand plot range slightly

# Set energy integration range
omega_width = 0.3
params.binstart[4] = target_ω - (omega_width/2)
params.binend[4] = target_ω # `binend` should be inside (e.g. at the center) of the range
params.binwidth[4] = omega_width

integrate_axes!(params, axes = 3) # Integrate out z direction entirely

params

is, counts = intensities_binned(sc,params,new_formula)

fig = Figure()
ax = Axis(fig[1,1];
    title="Δω=0.3 meV (Binned)", aspect=true,
    xlabel = "[H, 0, 0]",
    ylabel = "[0, K, 0]"
)
bcs = axes_bincenters(params)
hm = heatmap!(ax,bcs[1],bcs[2],is[:,:,1,1] ./ counts[:,:,1,1])
function add_lines!(ax,params)#hide
  bes = Sunny.axes_binedges(params)#hide
  hrange = range(-2,2,length=17)#hide
  linesegments!(ax,[(Point2f(params.covectors[1,1:3] ⋅ [h,-10,0],params.covectors[2,1:3] ⋅ [h,-10,0]),Point2f(params.covectors[1,1:3] ⋅ [h,10,0],params.covectors[2,1:3] ⋅ [h,10,0])) for h = hrange],linestyle=:dash,color=:black)#hide
  krange = range(-2,2,length=17)#hide
  linesegments!(ax,[(Point2f(params.covectors[1,1:3] ⋅ [-10,k,0],params.covectors[2,1:3] ⋅ [-10,k,0]),Point2f(params.covectors[1,1:3] ⋅ [10,k,0],params.covectors[2,1:3] ⋅ [10,k,0])) for k = krange],linestyle=:dash,color=:black)#hide
  xlims!(ax,bes[1][1],bes[1][end])#hide
  ylims!(ax,bes[2][1],bes[2][end])#hide
end#hide
add_lines!(ax,params)
Colorbar(fig[1,2], hm);
fig

latvecs = sys.crystal.latvecs
metric = latvecs' * I(3) * latvecs

(latvecs * [1/2,1,0]) ⋅ (latvecs * [1,0,0]) == 0

f = Figure()#hide
ax = Axis(f[1,1])#hide
arrows!(ax,[Point2f(0,0),Point2f(latvecs[1:2,1] ./ 2)],[Vec2f(latvecs[1:2,1] ./ 2), Vec2f(latvecs[1:2,2])],arrowcolor = :blue,arrowsize = 30.,linewidth = 5.,linecolor = :blue)#hide
arrows!(ax,[Point2f(0,0)],[Vec2f(latvecs[1:2,:] * [1/2,1,0])],arrowcolor = :red,arrowsize = 30.,linewidth = 5.,linecolor = :red, linestyle = :dash)#hide
scatter!(ax,[Point2f(latvecs[1:2,:] * [a,b,0]) for a in -1:1, b in -1:1][:],color = :black)#hide
annotations!(ax,["0","0+b","0+a", "a/2", "b"],[Point2f(0,-0.3),Point2f(latvecs[1:2,2]) .- Vec2f(0,0.3),Point2f(latvecs[1:2,1]) .- Vec2f(0,0.3),Point2f(latvecs[1:2,1] ./ 4) .- Vec2f(0,0.3),Point2f(latvecs[1:2,1] ./ 2) .+ Vec2f(latvecs[1:2,2] ./ 2) .+ Vec2f(0.3,0.3)],color=[:black,:black,:black,:blue,:blue])#hide
f#hide

params.covectors[2,1:3] = [1/2,1,0] # [1/2,1,0] times [a;b;c] is (a/2 + b)
params#hide

# Zoom out horizontal axis
params.binstart[1], params.binend[1] = -2, 2

# Adjust vertical axis bounds to account for
# length of a/2 + b
params.binstart[2], params.binend[2] = -2 * sqrt(3/4), 2 * sqrt(3/4)

# Re-compute in the new coordinate system
is, counts = intensities_binned(sc,params,new_formula)

fig = Figure(; size=(600,250))#hide
ax_right = Axis(fig[1,3];#hide
    title="ω≈$(round(target_ω, digits=2)) meV with Δω=0.3 meV (Binned)", aspect=true,#hide
    xlabel = "[H, -1/2H, 0]"#hide
)#hide
bcs = axes_bincenters(params)#hide
hm_right = heatmap!(ax_right,bcs[1],bcs[2],is[:,:,1,1] ./ counts[:,:,1,1])#hide
add_lines!(ax_right,params)
Colorbar(fig[1,4], hm_right);#hide

# New basis matrix
A = [1    0 0
     -1/2 1 0
     0    0 1]

# Define our grid of wave vectors
npoints = 60
as = range(-2, 2, npoints)
bs = range(-3/√3, 3/√3, npoints)
qs_ortho = [[a, b, 0] for a in as, b in bs]

# Convert to original RLU system for input to Sunny
qs = [A * q for q in qs_ortho]

# Use interpolation to get intensities
is = intensities_interpolated(sc, qs, new_formula; interpolation=:linear)

ax_left = Axis(fig[1,2];#hide
    title="ω≈$(round(ωs[ωidx], digits=2)) meV (Interpolated)", aspect=true,#hide
    xlabel = "[H, -1/2H, 0]", ylabel = "[0, K, 0]"#hide
)#hide
hm_left = heatmap!(ax_left, as, bs, is[:,:,ωidx])#hide
add_lines!(ax_left,params)
Colorbar(fig[1,1], hm_left);#hide
fig

metric = (latvecs * inv(A'))' * I(3) * (latvecs * inv(A'))

is_static = instant_intensities_interpolated(sc, qs, new_formula; interpolation = :linear)

hm = heatmap(as, bs, is_static;
    axis=(
        title="Instantaneous Structure Factor",
        xlabel = "[H, -1/2H, 0]",
        ylabel = "[0, K, 0]",
        aspect=true
    )
)
Colorbar(hm.figure[1,2], hm.plot)
hm
