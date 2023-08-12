# # Structure Factors with Classical Dynamics 

using Sunny, GLMakie#hide

# In our previous [Case Study: FeI$_{2}$](@ref), we used linear spin wave theory
# (LSWT) to calculate the dynamical structure factor. Here, we perform a similar
# calculation using classical spin dynamics. Because we are interested in the
# coupled dynamics of spin dipoles and quadrupoles, we must employ a [classical
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
# equilibrium and used as initial conditions for generating a dissipationless
# trajectories. The correlations of these trajectories are then averaged and
# used to calculate the structure factor. It is therefore important to ensure
# that the initial spin configurations are sampled appropriately and that
# sufficient statistics are collected. We will demonstrate one approach here.
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

# The [Langevin dynamics of SU(_N_) coherent
# states](https://arxiv.org/abs/2209.01265) can be used to sample spin
# configurations in thermal equlibrium. One first constructs a
# [`Langevin`](@ref) integrator. This requires a time step, temperature, and a
# phenomenological damping parameter $λ$ that sets the coupling to the thermal
# bath.

Δt = 0.05/D    # A good rule of thumb is `Δt = 0.05/E₀`, where E₀ is the
               ## largest energy scale in the system. In FeI2, this is the
               ## easy-axis anisotropy scale, `D = 2.165` (1/meV). 
kT = 0.2       # Temperature for spin fluctuations (meV).
λ = 0.1        # A good value is 0.1, independent of system details.
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

plot_spins(sys; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# In this case, we can find the correct ground state simply by running the
# Langevin sampler for longer.

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

# # Calculating $S(𝐪,ω)$

# Because classical simulations are conducted on a finite-sized lattice,
# obtaining acceptable resolution in momentum space requires the use of a larger
# system size. We can now resize the magnetic supercell to a much larger
# simulation volume, provided as multiples of the original unit cell.

sys_large = resize_supercell(sys, (16,16,4))
plot_spins(sys_large; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# As stressed above, the calculation process depends on sampling equilibrium
# spin configurations. We therefore begin by using Langevin dynamics to
# thermalize the system to a target temperature.

kT = 3.5 * meV_per_K     # 3.5K, converted to meV
langevin.kT = kT

for _ in 1:10_000
    step!(sys_large, langevin)
end

# To estimate the dynamic structure factor, we can collect spin-spin correlation
# data by first generating an initial condition at thermal equilibrium and then
# integrating the [Hamiltonian dynamics of SU(_N_) coherent
# states](https://arxiv.org/abs/2204.07563). Samples are accumulated into a
# `SampledCorrelations`, which is initialized by calling
# [`dynamical_correlations`](@ref).  `dynamical_correlations` takes a `System`
# and three keyword parameters that determine the ω information that will be
# available: an integration step size, the number of ωs to resolve, and the
# maximum ω to resolve. For the time step, twice the value used for the Langevin
# integrator is usually a good choice.

sc = dynamical_correlations(sys_large; Δt=2Δt, nω=120, ωmax=7.5)

# `sc` currently contains no data. A sample can be accumulated into it by
# calling [`add_sample!`](@ref).

add_sample!(sc, sys_large)

# Additional samples can be added after generating new spin configurations,
# which will again be achieved by running the Langevin dynamics. It is important
# that enough steps are chosen so that the resulting samples are decorrelated.

for _ in 1:2
    for _ in 1:1000               # Enough steps to decorrelate spins
        step!(sys_large, langevin)
    end
    add_sample!(sc, sys_large)    # Accumulate the sample into `sc`
end

# ## Accessing structure factor data

# The basic functions for accessing intensity data are [`intensities_interpolated`](@ref)
# and [`intensities_binned`](@ref). Both functions accept an [`intensity_formula`](@ref)
# to specify how to combine the correlations recorded in the `SampledCorrelations` into
# intensity data. By default, a formula computing the unpolarized intensity is used,
# but alternative formulas can be specified.
#
# By way of example, we will use a formula which computes the trace of the structure
# factor and applies a classical-to-quantum temperature-dependent rescaling `kT`.

formula = intensity_formula(sc, :trace; kT = kT)

# Using the formula, we plot single-$q$ slices at (0,0,0) and (π,π,π):

qs_rlu = [[0, 0, 0], [0.5, 0.5, 0.5]]
qs = [cryst.recipvecs*q for q in qs_rlu]  # Convert from RLU to absolute units
is = intensities_interpolated(sc, qs; interpolation = :round, formula = formula)

fig = lines(ωs(sc), is[1,:]; axis=(xlabel="meV", ylabel="Intensity"), label="(0,0,0)")
lines!(ωs(sc), is[2,:]; label="(π,π,π)")
axislegend()
fig

# To get smoother results, it will be necessary to collect many more statistics.
#
# For real calculations, one often wants to apply further corrections and more
# accurate formulas Here, we apply [`FormFactor`](@ref) corrections appropriate
# for `Fe2` magnetic ions, and a dipole polarization correction `:perp`.

formfactors = [FormFactor("Fe2"; g_lande=3/2)]
new_formula = intensity_formula(sc, :perp; kT = kT, formfactors = formfactors)

# Frequently one wants to extract energy intensities along lines that connect
# special wave vectors. The function [`connected_path_from_rlu`](@ref) linearly
# samples between provided $q$-points, with a given sample density.

points = [[0,   0, 0],  # List of wave vectors that define a path
          [1,   0, 0],
          [0,   1, 0],
          [1/2, 0, 0],
          [0,   1, 0],
          [0,   0, 0]] 
density = 40
path, xticks = connected_path_from_rlu(cryst, points, density);

# Since scattering intensities are only available at a certain discrete ``(Q,\omega)``
# points, the intensity on the path can be calculated by interpolating between these
# discrete points:

is = intensities_interpolated(sc, path;
    interpolation = :linear,       # Interpolate between available wave vectors
    formula = new_formula
)
is = broaden_energy(sc, is, (ω, ω₀)->lorentzian(ω-ω₀, 0.05))  # Add artificial broadening
labels = string.(points)
heatmap(1:size(is,1), ωs(sc), is;
    axis = (
        ylabel = "meV",
        xticklabelrotation=π/8,
        xticklabelsize=12,
        xticks,
    ),
    colorrange=(0.0,0.7),
)

# Note that we have clipped the colors in order to make the higher-energy
# excitations more visible.
#
# Whereas [`intensities_interpolated`](@ref) either rounds or linearly
# interpolates between the discrete ``(Q,\omega)`` points Sunny calculates
# correlations at, [`intensities_binned`](@ref) performs histogram binning
# analgous to what is done in experiments. The graphs produced by each method
# are similar.
cut_width = 0.3
density = 15
paramsList, markers, ranges = connected_path_bins(sc,points,density,cut_width)

total_bins = ranges[end][end]
energy_bins = paramsList[1].numbins[4]
is = zeros(Float64,total_bins,energy_bins)
integrated_kernel = integrated_lorentzian(0.05) # Lorentzian broadening
for k in eachindex(paramsList)
    h,c = intensities_binned(sc,paramsList[k];formula = new_formula,integrated_kernel = integrated_kernel)
    is[ranges[k],:] = h[:,1,1,:] ./ c[:,1,1,:]
end

heatmap(1:size(is,1), ωs(sc), is;
    axis = (
        ylabel = "meV",
        xticks = (markers, labels),
        xticklabelrotation=π/8,
        xticklabelsize=12,
    ),
    colorrange=(0.0,0.05),
)


# Often it is useful to plot cuts across multiple wave vectors but at a single
# energy. 

npoints = 60
qvals = range(-3, 3, length=npoints)
qs = [[a, b, 0] for a in qvals, b in qvals]

is = intensities_interpolated(sc, qs; formula = new_formula,interpolation = :linear);

ωidx = 60
hm = heatmap(is[:,:,ωidx]; axis=(title="ω=$(ωs(sc)[ωidx]) meV", aspect=true))
Colorbar(hm.figure[1,2], hm.plot)
hidedecorations!(hm.axis); hidespines!(hm.axis)
hm

# Finally, we note that instantaneous structure factor data, ``𝒮(𝐪)``, can be
# obtained from a dynamic structure factor with [`instant_intensities_interpolated`](@ref).

is_static = instant_intensities_interpolated(sc, qs; formula = new_formula, interpolation = :linear)

hm = heatmap(is_static; axis=(title="Instantaneous Structure Factor", aspect=true))
Colorbar(hm.figure[1,2], hm.plot)
hidedecorations!(hm.axis); hidespines!(hm.axis)
hm