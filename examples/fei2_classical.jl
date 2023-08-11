# # Using Classical Dynamics to Calcualte $S(\mathbf{q},\omega)$.

using Sunny, GLMakie

# The FeI$_2$ tutorial shows how to use traditional linear spin wave theory to
# calculate dynamical information in Sunny. Here we show how similar information
# may be calculated using classical dynamics, specifically by employing a
# generalization of Landau-Lifshitz dynamics to coherent states of SU($3$). This
# approach is generallly computationally more expensive than the traditional
# approach, and it can only be performed for finite-sized systems. The classical
# approach nonetheless provides a number of complementary advantages: it is
# possible perform simulations at finite temperature while retaining
# nonlinearities; out-of-equilibrium behavior may be examined directly; and it
# is straightforward to incorporate inhomogenties, chemical or otherwise.  
#
# In this tutorial, we show how to calculate the finite temperature dynamics of
# FeI$_2$ using the classical approach. It is important to stress that
# estimating ``S(𝐪,ω)`` with classical dynamics is fundamentally a Monte Carlo
# process: sample spin configurations are drawn from thermal equilibrium and
# used as initial conditions for generating a dissipationless trajectories. The
# correlations of these trajectories are then averaged and used to calculate the
# structure factor. It is therefore important to ensure that the initial spin
# configurations are sampled appropriately and that sufficient statistics are
# collected. We will demonstrate one approach here.
#
# As the implementation of FeI$_2$ is already covered in detail in the previous
# tutorial, we will not repeat it below. Instead, we will assume that you
# already have a $4\times 4\times 4$ `sys` available and proceed to the
# calculation.

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
# 
# The [Langevin dynamics of SU(_N_) coherent
# states](https://arxiv.org/abs/2209.01265) can be used to sample spin
# configurations in thermal equlibrium. Begin by constructing a
# [`Langevin`](@ref) integrator.

Δt = 0.05/D    # Single-ion anisotropy is the strongest interaction, so 1/D is
               # a natural dynamical time-scale (in units with ħ=1).
λ = 0.1        # Dimensionless magnitude of coupling to thermal bath
langevin = Langevin(Δt; kT=0, λ);

# Attempt to find a low-energy spin configuration by lowering the temperature kT
# from 2 to 0 using 20,000 Langevin time-steps.
randomize_spins!(sys)
for kT in range(2, 0, 20_000)
    langevin.kT = kT
    step!(sys, langevin)
end

# Because the quench was relatively fast, it is expected to find defects in the
# magnetic order. These can be visualized.

plot_spins(sys; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# If we had used a slower annealing procedure, involving 100,000 or more
# Langevin time-steps, it would very likely find the correct ground state.
# Instead, for purposes of illustration, let's analyze the imperfect spin
# configuration currently stored in `sys`.
#
# An experimental probe of magnetic order order is the 'instantaneous' or
# 'static' structure factor intensity, available via
# [`instant_correlations`](@ref) and related functions. To infer periodicities
# of the magnetic supercell, however, it is sufficient to look at the structure
# factor weights of spin sublattices individually, _without_ phase averaging.
# This information is provided by [`print_wrapped_intensities`](@ref) (see the
# API documentation for a physical interpretation).

print_wrapped_intensities(sys)

# The above is consistent with known results. The zero-field energy-minimizing
# magnetic structure of FeI$_2$ is single-$q$. If annealing were perfect, then
# spontaneous symmetry breaking would select one of $±q = [0, -1/4, 1/4]$,
# $[1/4, 0, 1/4]$, or $[-1/4,1/4,1/4]$.

# Let's break the symmetry by hand. Given a list of $q$ modes, Sunny can suggest
# a magnetic supercell with appropriate periodicity. The result is printed in
# units of the crystal lattice vectors.

suggest_magnetic_supercell([[0, -1/4, 1/4]], sys.latsize)

# The function [`reshape_supercell`](@ref) allows an arbitrary reshaping of the
# system. After selecting the supercell geometry, it becomes much easier to find
# the energy-minimizing spin configuration.

sys_supercell = reshape_supercell(sys, [1 0 0; 0 1 -2; 0 1 2])

langevin.kT = 0
for i in 1:10_000
    step!(sys_supercell, langevin)
end

plot_spins(sys_supercell; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# Instead of using Langevin integration to perform the final minimization, we
# could also have used [`minimize_energy`](@ref).

## Calculating ``S(𝐪,ω)``
# Because classical simulations are conducted on a finite-sized lattice,
# obtaining acceptable resolution in momentum space requires the use of a larger
# system size. We can now resize the magnetic supercell to a much larger
# simulation volume, provided as multiples of the original unit cell.

sys_large = resize_supercell(sys_supercell, (16,16,4))
plot_spins(sys_large; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# As stressed above, the calculation process depends on sampling equilibrium
# spin configurations. We therefore begin by using Langevin dynamics to
# thermalize the system to a target temperature.

kT = 3.5 * meV_per_K     # 3.5K in units of meV
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
# special wave vectors. The function [`connected_path`](@ref) linearly samples
# between provided $q$-points, with a given sample density.

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