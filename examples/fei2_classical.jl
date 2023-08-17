# # Structure Factors with Classical Dynamics 

using Sunny, GLMakie#hide

# In our previous [Case Study: FeI$_{2}$](@ref), we used linear spin wave theory
# (LSWT) to calculate the dynamical structure factor. Here, we perform a similar
# calculation using classical spin dynamics. Because we are interested in the
# coupled dynamics of spin dipoles and quadrupoles, we employ a [classical
# dynamics of SU(3) coherent states](https://arxiv.org/abs/2209.01265) that
# generalizes the Landau-Lifshitz equation.
#
# Compared to LSWT, simulations using classical dynamics are much slower, and
# are limited in $k$-space resolution. However, they make it is possible to
# capture nonlinear effects associated with finite temperature fluctuations.
# Classical dynamics are also appealing for studying out-of-equilibrium systems
# (e.g., relaxation of spin glasses), or systems with quenched inhomogeneities
# that require large simulation volumes.
#
# In this tutorial, we show how to study the finite temperature dynamics of
# FeI$_2$ using the classical approach. It is important to stress that the
# estimation of ``S(𝐪,ω)`` with classical dynamics is fundamentally a Monte
# Carlo calculation: sample spin configurations are drawn from thermal
# equilibrium and used as initial conditions for generating dissipationless
# trajectories. The correlations of these trajectories are then averaged and
# used to calculate scattering intensities. It is therefore important to ensure
# that the initial spin configurations are sampled appropriately and that
# sufficient statistics are collected. We will demonstrate one approach here.
#
# As an overview, we will:
#
# 1. Identify the ground state
# 2. Measure correlation data describing the excitations around that ground state
# 3. Use the correlation data to compute scattering intensities
#
# As the implementation of the FeI$_2$ model is already covered in detail in the
# LSWT tutorial, we will not repeat it below. Instead, we will assume that you
# already have defined a `sys` in the same way with lattice dimensions $4×4×4$. 

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
S = spin_operators(sys, 1)#hide
set_onsite_coupling!(sys, -D*S[3]^2, 1)#hide

# ## Finding a ground state

# Sunny uses the [Langevin dynamics of SU(_N_) coherent
# states](https://arxiv.org/abs/2209.01265) to sample spin configurations
# from the thermal equlibrium. One first constructs a
# [`Langevin`](@ref) integrator. This requires a time step, temperature, and a
# phenomenological damping parameter $λ$ that sets the coupling to the thermal
# bath.

Δt = 0.05/D    # Should be inversely proportional to the largest energy scale
               ## in the system. For FeI2, this is the easy-axis anisotropy,
               ## `D = 2.165` (meV). The prefactor 0.05 is relatively small,
               ## and achieves high accuracy.
kT = 0.2       # Temperature of the thermal bath (meV).
λ = 0.1        # This value is typically good for Monte Carlo sampling,
               ## independent of system details.

langevin = Langevin(Δt; kT, λ);

# Sampling with a large number of Langevin time-steps should hopefully
# thermalize the system to a target temperature. For demonstration purposes, we
# will here run for a relatively short period of 1,000 timesteps.

randomize_spins!(sys)
for _ in 1:1_000
    step!(sys, langevin)
end

# To inspect the structure of the spin configuration, we use the function
# [`minimize_energy!`](@ref) to find a nearby _local_ minimum. Then
# [`print_wrapped_intensities`](@ref) provides information about the Fourier
# modes in reciprocal lattice units (RLU).

minimize_energy!(sys)
print_wrapped_intensities(sys)

# The wide distribution of intensities indicates an imperfect magnetic order.
# The defects are immediately apparent in the real-space spin configuration.

plot_spins(sys)

# In this case, we can find the correct ground state simply by running the
# Langevin dynamics for longer.

for _ in 1:10_000
    step!(sys, langevin)
end
minimize_energy!(sys)
print_wrapped_intensities(sys)

# This is the perfect single-$𝐐$ ground state. This worked because the
# temperature `kT = 0.2` was carefully selected. It is below the magnetic
# ordering temperature, but still large enough that the Langevin dynamics could
# quickly overcome local energy barriers. More complicated magnetic orderings
# will frequently require more sophisticated sampling schemes. One possibility
# is simulated annealing, where `kT` is slowly lowered over the course of the
# sampling procedure.

# ## Calculating Thermal-Averaged Correlations $\langle S^{\alpha\beta}(𝐪,ω)\rangle$
#
# Now that we have identified an appropriate ground state, we will adjust
# the temperature so that the system still remains near the ground state, but still
# has thermal fluctuations.

kT = 3.5 * meV_per_K     # 3.5K ≈ 0.30 meV
langevin.kT = kT;

# Additionally, since these classical simulations are conducted on a finite-sized lattice,
# obtaining acceptable resolution in momentum space requires the use of a larger
# system size. We resize the magnetic supercell to a much larger
# simulation volume, provided as multiples of the original unit cell.

sys_large = resize_supercell(sys, (16,16,4)) # 16x16x4 copies of the original unit cell
plot_spins(sys_large)

# As stressed above, we need to sample multiple spin configurations
# from the thermal equilibrium distribution.
# We therefore begin by using Langevin dynamics to
# bring the system into thermal equilibrium at the new temperature.

## At the new temperature
for _ in 1:10_000
    step!(sys_large, langevin)
end

# Now that we are in a state drawn from thermal equilibrium, we are ready to begin
# collecting correlation data ``S^{\alpha\beta}``.
#
# To collect one sample of spin-spin correlation data, we integrate the
# [Hamiltonian dynamics of SU(_N_) coherent states](https://arxiv.org/abs/2204.07563).
# This generates trajectories ``S^\alpha(\vec r,t)``.
# However, note that the spins are only defined at the finitely-many lattice
# sites, so ``\vec r`` is discrete.
# From the trajectories, Sunny computes fourier-transformed correlations ``S^{\alpha\beta}(q,\omega)``,
# where ``q`` is discrete for the same reason.
#
# The correlation data from this sample is now ready to be averaged together with data
# from other samples to form the thermal average.
# Samples of correlation data are accumulated and averaged into a
# `SampledCorrelations`, which is initialized by calling
# [`dynamical_correlations`](@ref):

sc = dynamical_correlations(sys_large; Δt=2Δt, nω=120, ωmax=7.5)

# `dynamical_correlations` takes a `System`
# and three keyword parameters that determine the ω information that will be
# available: an integration step size, the number of ωs to resolve, and the
# maximum ω to resolve. For the time step, twice the value used for the Langevin
# integrator is usually a good choice.

# `sc` currently contains no data. The first sample can be accumulated into it by
# calling [`add_sample!`](@ref).

add_sample!(sc, sys_large)

# Additional samples can be added after re-sampling new 
# spin configurations from the thermal distribution.
# The new samples are generated in the same way as the first sample, by running
# the Langevin dynamics.
# The dynamics should be run long enough that consecutive samples are uncorrelated, or
# else the thermal average will converge only very slowly.

for _ in 1:2
    for _ in 1:1000               # Enough steps to decorrelate spins
        step!(sys_large, langevin)
    end
    add_sample!(sc, sys_large)    # Accumulate the sample into `sc`
end

# Now, `sc` has more samples included:

sc

# ## Computing Scattering Intensities

# With the thermally-averaged correlation data ``\langle S^{\alpha\beta}(q,\omega)\rangle``
# in hand, we now need to specify how to extract a scattering intensity from this information.
# This is done by constructing an [`intensity_formula`](@ref).
# By way of example, we will use a formula which computes the trace of the structure
# factor and applies a classical-to-quantum temperature-dependent rescaling `kT`.

formula = intensity_formula(sc, :trace; kT = kT)

# Recall that ``\langle S^{\alpha\beta}(q,\omega)\rangle`` is only available at certain discrete
# ``q`` values, due to the finite lattice size.
# There are two basic approaches to handling this discreteness.
# The first approach is to interpolate between the available
# data using [`intensities_interpolated`](@ref). For example, we can plot single-$q$ slices
# at (0,0,0) and (π,π,π) using this method:

qs = [[0, 0, 0], [0.5, 0.5, 0.5]]
is = intensities_interpolated(sc, qs, formula; interpolation = :round)

ωs = available_energies(sc)
fig = lines(ωs, is[1,:]; axis=(xlabel="meV", ylabel="Intensity"), label="(0,0,0)")
lines!(ωs, is[2,:]; label="(π,π,π)")
axislegend()
fig

# The resolution in energy can be improved by increasing `nω` (and decreasing `Δt`),
# and the general accuracy can be improved by collecting additional samples from the thermal
# equilibrium.
#
# For real calculations, one often wants to apply further corrections and more
# accurate formulas. Here, we apply [`FormFactor`](@ref) corrections appropriate
# for `Fe2` magnetic ions, and a dipole polarization correction `:perp`.

formfactors = [FormFactor("Fe2"; g_lande=3/2)]
new_formula = intensity_formula(sc, :perp; kT = kT, formfactors = formfactors)

# Frequently, one wants to extract energy intensities along lines that connect
# special wave vectors--a so-called "spaghetti plot". The function 
# [`reciprocal_space_path`](@ref) creates an appropriate horizontal axis for this plot
# by linearly sampling between provided $q$-points, with a given sample density.

points = [[0,   0, 0],  # List of wave vectors that define a path
          [1,   0, 0],
          [0,   1, 0],
          [1/2, 0, 0],
          [0,   1, 0],
          [0,   0, 0]] 
density = 40
path, xticks = reciprocal_space_path(cryst, points, density);

# Again using [`intensities_interpolated`](@ref), we can evaluate the (interpolated) intensity
# at each point on the `path`.
# Since scattering intensities are only available at a certain discrete ``(Q,\omega)``
# points, the intensity on the path can be calculated by interpolating between these
# discrete points:

is_interpolated = intensities_interpolated(sc, path, new_formula;
    interpolation = :linear,       # Interpolate between available wave vectors
);
## Add artificial broadening
is_interpolated_broadened = broaden_energy(sc, is, (ω, ω₀)->lorentzian(ω-ω₀, 0.05));

# The second approach to handle the discreteness of the data is to bin the intensity
# at the discrete points into the bins of a histogram.
# First, the five sub-histograms are set up using [`reciprocal_space_path_bins`](@ref) in
# analogy to `reciprocal_space_path`.

cut_width = 0.3
density = 15
paramsList, markers, ranges = reciprocal_space_path_bins(sc,points,density,cut_width);

# Then, the intensity data is computed using [`intensities_binned`](@ref) for each sub-histogram:

total_bins = ranges[end][end]
energy_bins = paramsList[1].numbins[4]
is_binned = zeros(Float64,total_bins,energy_bins)
integrated_kernel = integrated_lorentzian(0.05) # Lorentzian broadening
for k in eachindex(paramsList)
    bin_data, counts = intensities_binned(sc,paramsList[k], new_formula;
        integrated_kernel = integrated_kernel
    )
    is_binned[ranges[k],:] = bin_data[:,1,1,:] ./ counts[:,1,1,:]
end

# The graph produced by interpolating (top) is similar to the one produced by binning (bottom):

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


# Note that we have clipped the colors in order to make the higher-energy
# excitations more visible.

# Often it is useful to plot cuts across multiple wave vectors but at a single
# energy. 
ωidx = 60
target_ω = ωs[ωidx]

## Binning
fig = Figure()
ax_left = Axis(fig[1,2],title="Δω=0.3 meV (Binned)", aspect=true)

params = unit_resolution_binning_parameters(sc)
params.covectors[1:3,1:3] .= cryst.recipvecs  # We will express bin limits in absolute units
omega_width = 0.3
params.binstart[4] = target_ω
params.binend[4] = target_ω + (sc.Δω/2)
params.binwidth[4] = omega_width
params.binwidth[1:2] .*= 1.78  # Increase bin size to avoid empty bins

params.binstart[1], params.binstart[2] = -3, -3
params.binend[1], params.binend[2] = 3, 3

integrate_axes!(params, axes = [3,4])

integrate_axes!(params,axes = [3,4])
is, counts = intensities_binned(sc,params,new_formula)
bcs = axes_bincenters(params)
hm_left = heatmap!(ax_left,bcs[1],bcs[2],is[:,:,1,1] ./ counts[:,:,1,1])
Colorbar(fig[1,1], hm_left);

## Interpolating
ax_right = Axis(fig[1,3],title="ω≈$(ωs[ωidx]) meV (Interpolated)", aspect=true)
npoints = 60
qvals = range(-3, 3, length=npoints)
qs_absolute = [[a, b, 0] for a in qvals, b in qvals]
qs = [cryst.recipvecs \ q for q in qs_absolute]

is = intensities_interpolated(sc, qs, new_formula; interpolation=:linear)

hm_right = heatmap!(ax_right,is[:,:,ωidx])
Colorbar(fig[1,4], hm_right)
hidedecorations!(ax_right); hidespines!(ax_right)
fig

# Finally, we note that instantaneous structure factor data, ``𝒮(𝐪)``, can be
# obtained from a dynamic structure factor with [`instant_intensities_interpolated`](@ref).

is_static = instant_intensities_interpolated(sc, qs, new_formula; interpolation = :linear)

hm = heatmap(is_static; axis=(title="Instantaneous Structure Factor", aspect=true))
Colorbar(hm.figure[1,2], hm.plot)
hidedecorations!(hm.axis); hidespines!(hm.axis)
hm