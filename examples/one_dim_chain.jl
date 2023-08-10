using Sunny, LinearAlgebra, GLMakie, Optim

# 1D FM chain (J1=0.1D) with single-ion anisotropy D. Data created with Landau-Lifshitz at finite but small temperature. Fit using SWT kernel to extract J and D.

# In this Example, we consider a 1D chain of spin-1 sites.
# The sites along the chain interact via a ferromagnetic nearest-neighbor
# interaction ``J\sum_{\langle i,j\rangle} \mathbf{S}_i \cdot \mathbf{S}_j``, with ``J < 0``.
# By default, the ground state would be ferromagnetic and highly degenerate, since the spins
# can align in any direction.
# An on-site interaction, ``D\sum_i (S^z_i)^2`` breaks this isotropy by making it easier for
# the spins to align in the ``\pm z`` direction than in any other orientation.
# Thus, the entire Hamiltonian is:
#
# ``\mathcal{H} = \overbrace{J\sum_{\langle i,j\rangle} \mathbf{S}_i \cdot \mathbf{S}_j}^{\text{Ferromagnetic}}\;\;\;\;\;\; \overbrace{-D\sum_i (S^z_i)^2}^{\text{Easy-axis single-ion anisotropy}}``
#
# The goal of this Example is to illustrate how to determine the parameters ``J`` and ``D``
# from "experiment" data by fitting using Sunny's implementation of Linear Spin Wave Theory.
# In our case, the "experiment" data will actually be simulation data produced
# using Landau-Lifschitz dynamics.

# # Creating simulated "experiment" data using Landau-Lifschitz dynamics

# Our simulated data will use ground truth values ``J_0 = -1\,\text{meV}`` and ``D_0 = 10\,\text{meV}`` with a lattice spacing ``a = 10`` angstrom.

# We begin with a 1D chain of spin-1 sites along the ``x`` direction.

## Establish geometry of the unit cell.
## "P1" is required due to the rotational symmetry about the
## x-axis being broken.
chain_spacing = 10. # Angstrom
latvecs = chain_spacing * I(3)
one_dimensional_chain = Crystal(latvecs,[[0,0,0]],"P1")

## Establish geometry of the whole chain.
chain_length = 64 # Number of atoms
latsize = (chain_length,1,1) # 1D chain is Nx1x1 lattice
spin_one_chain = System(one_dimensional_chain, latsize, [SpinInfo(1,S=1)], :SUN)

# Configure the nearest-neighbor interaction:

## Scalar J indicates J*(Sᵢ⋅Sⱼ)
J_groundtruth = -1.

## Interaction is with the left and right next neighbor along the chain (x-direction)
nearest_neighbor_right = Bond(1,1,(1,0,0))
nearest_neighbor_left = Bond(1,1,(-1,0,0))

set_exchange!(spin_one_chain,J_groundtruth,nearest_neighbor_right)
set_exchange!(spin_one_chain,J_groundtruth,nearest_neighbor_left)

# Configure the symmetry-breaking easy-axis term:
D_groundtruth = 10.
Sz = spin_operators(spin_one_chain, 1)[3]
set_onsite_coupling!(spin_one_chain, -D_groundtruth*Sz^2, 1)

# With the ground-truth hamiltonian in place, we use Sunny's classical dynamics
# to generate ficticious experiment data at temperature `kT = 0.1`.

Δt = 0.05/D_groundtruth
λ = 0.1
kT = 0.1
langevin = Langevin(Δt; kT, λ);

function viz_chain(sys;kwargs...)#hide
  ups = map(x -> abs2(x[1]), sys.coherents)[:];#hide
  zs = map(x -> abs2(x[2]), sys.coherents)[:];#hide
  downs = map(x -> abs2(x[3]), sys.coherents)[:];#hide
#hide
  f = Figure()#hide
  ax = LScene(f[1,1];show_axis = false)#hide
  _ = Makie.cam3d!(ax.scene, projectiontype=Makie.Orthographic)#hide
#hide
  linewidth = 5.#hide
  arrowsize = 10.#hide
  lengthscale = 15.#hide
  pts = [Point3f(Sunny.global_position(sys,site)) for site in eachsite(sys)][:]#hide
#hide
  ## Ups#hide
  vecs = [Vec3f([0,0,1]) for site in eachsite(sys)][:]#hide
  cols = map(x -> (:blue,x), ups)#hide
  Makie.arrows!(ax, pts .+ 0.5 .* vecs, vecs;#hide
        linecolor = cols, arrowcolor = cols,#hide
        lengthscale, arrowsize, linewidth, kwargs...)#hide
#hide
  ## Downs#hide
  vecs = [Vec3f([0,0,-1]) for site in eachsite(sys)][:]#hide
  cols = map(x -> (:red,x), downs)#hide
  Makie.arrows!(ax, pts .+ 0.5 .* vecs, vecs;#hide
        linecolor = cols, arrowcolor = cols,#hide
        lengthscale, arrowsize, linewidth, kwargs...)#hide
#hide
  cols = map(x -> (:green,x), zs)#hide
  meshscatter!(ax,pts, markersize = 7., color = cols)#hide
  f#hide
end#hide
randomize_spins!(spin_one_chain)
viz_chain(spin_one_chain)

# First, we thermalize the chain, and then take several samples
# in order get reasonably good "experiment" data.

nStep = 50_000#hide
for _ in 1:nStep#hide
    step!(spin_one_chain, langevin)#hide
end#hide
## ... thermalize ...
viz_chain(spin_one_chain)

#

sc = dynamical_correlations(spin_one_chain; Δt, nω = 80, ωmax = 20.);#hide

for _ in 1:10_000#hide
    step!(spin_one_chain, langevin)#hide
end#hide
add_sample!(sc, spin_one_chain)#hide
## ... some time later ...
viz_chain(spin_one_chain)

#

for _ in 1:10_000#hide
    step!(spin_one_chain, langevin)#hide
end#hide
add_sample!(sc, spin_one_chain)#hide
## ... some time later ...
viz_chain(spin_one_chain)

#

for _ in 1:20#hide
    for _ in 1:10_000#hide
        step!(spin_one_chain, langevin)#hide
    end#hide
    add_sample!(sc, spin_one_chain)#hide
end#hide
## ... some time later ...
viz_chain(spin_one_chain)

# Now that we have collected several samples,

sc

# we are ready to generate the intensity data.
# Since this is supposed to represent an experiment, the intensity data will go in a histogram:

SIMULATED_EXPERIMENT_HISTOGRAM_PARAMS = unit_resolution_binning_parameters(sc)

# Here's what the experiment data looks like:

formula = intensity_formula(sc;kT)
is, counts = intensities_binned(sc,SIMULATED_EXPERIMENT_HISTOGRAM_PARAMS;formula)

SIMULATED_EXPERIMENT_DATA = (is ./ counts)[:,1,1,:]

bcs = axes_bincenters(SIMULATED_EXPERIMENT_HISTOGRAM_PARAMS)
f = Figure()#hide
ax = Axis(f[1,1])#hide
heatmap!(ax,bcs[1],bcs[4],log10.(SIMULATED_EXPERIMENT_DATA))
f#hide

# # Fitting to the experiment data

# To fit this data, we first model the known aspects of the system in Sunny.
# The first steps are the same whether we are simulating a known system or modelling an
# unknown system:

## Same as before
chain_spacing = 10. # Angstrom
latvecs = chain_spacing * I(3)
one_dimensional_chain = Crystal(latvecs,[[0,0,0]],"P1")
chain_length = 64 # Number of atoms
latsize = (chain_length,1,1) # 1D chain is Nx1x1 lattice
spin_one_chain = System(one_dimensional_chain, latsize, [SpinInfo(1,S=1)], :SUN)



# Originally, the next step would have been to configure the hamiltonian by
# specifying the ``J`` and ``D`` values. However, since these are unknowns,
# we will avoid using them as long as possible, and instead proceed to set up the
# bonds, spin operators, and `Langevin` integrator--none of which require the values:

Δt = 0.05
λ = 0.1
kT = 0. # LSWT uses zero temperature
langevin = Langevin(Δt; kT, λ);

nearest_neighbor_right = Bond(1,1,(1,0,0))
nearest_neighbor_left = Bond(1,1,(-1,0,0))

Sz = spin_operators(spin_one_chain, 1)[3]


# After this setup work is done *once*, we create a function `forward_problem(J_trial,D_trial)`
# which will compute the Linear Spin Wave Theoretic spectrum at the trial values of the
# ``J`` and ``D`` fitting parameters.
# In other words, the part of the original calculation which depends on the fitting
# parameters gets wrapped into a function:

function forward_problem(J_trial, D_trial)

  ## Ensure there is no phase transition (or else LSWT will throw errors)
  J_trial = min(J_trial,0)
  D_trial = max(D_trial,0)

  ## Uses J_trial
  set_exchange!(spin_one_chain,J_trial,nearest_neighbor_right)
  set_exchange!(spin_one_chain,J_trial,nearest_neighbor_left)

  ## Uses D_trial
  set_onsite_coupling!(spin_one_chain, -D_trial*Sz^2, 1)

  ## Perform spin wave calculation, continued below...

#+
# Note that `forward_problem` refers to variables defined outside
# of the scope of the function. This allows us to reuse those variables in each call to
# `forward_problem`, without reconstructing them each time. In general, the more that is known
# about the system you are modelling, the later in the code `function forward_problem(...)` can be
# inserted, and the more setup work can be re-used.
#
# !!! tip "`forward_problem` is a closure"
#     In computer progrogramming parlance, `forward_problem` is said to 'capture' variables such
#     as `spin_one_chain` from the enviroment. Since the result of calling `forward_problem` depends
#     not only on `J_trial` and `D_trial`, but also on `spin_one_chain`, it's no longer a function
#     of only its arguments.
#
#     Since `forward_problem` is not a closed system, but
#     `forward_problem + (captured variables)` _is_ a closed system, the latter is called 
#     the 'closure' of the former.

# ## Spin wave calculation

# We can leverage our knowledge that the ground state should be ferromagnetic
# to simplify the spin wave calculation. Since the ferrommagnetic unit cell is just one site,
# the simplified system is extremely simple:

  ## ... perform spin wave calculation, continued from above.
  one_site_system = reshape_geometry(spin_one_chain,[1 0 0; 0 1 0; 0 0 1])

#+
# After restricting to a single site, it's best to re-thermalize the system
# at zero temperature to ensure a good classical ground state for LSWT:

  langevin.kT = 0.
  nStep = 1_000
  for _ in 1:nStep
      step!(one_site_system, langevin)
  end

#+
# The spin wave intensity data must be placed in a histogram with the same parameters
# as the experiment data, in order to ensure a good comparision.
#
# The kernel and `intensities_bin_centers` used here are temporary, until a better
# binning method is written.

  swt = SpinWaveTheory(one_site_system)
  formula = intensity_formula(swt; kernel = lorentzian(0.5))
  params = SIMULATED_EXPERIMENT_HISTOGRAM_PARAMS
  is_swt = Sunny.intensities_bin_centers(swt, params, formula)

  return is_swt[:,1,1,:]
end # end of forward_problem

# We can see the different possible results from LSWT by plotting the dispersion:

function plot_forward(J,D)
  is_swt = forward_problem(J,D)
  bcs = axes_bincenters(SIMULATED_EXPERIMENT_HISTOGRAM_PARAMS)
  heatmap(bcs[1],bcs[4],log10.(is_swt))
end

plot_forward(-1,10)

#

plot_forward(-6,2)

#

plot_forward(-0.01,15)

# Now, we can easily define a least-squares loss function comparing the "experiment" data to the LSWT result:

function get_loss(parameters)
  J,D = parameters
  is_swt = forward_problem(J,D)
  sqrt(sum(abs2.(SIMULATED_EXPERIMENT_DATA .- is_swt)))
end

# Sweeping the parameters over a range containing the true value reveals
# that the loss is minimized near the true parameters (dot).
# The minimum loss is not exactly at the ground truth parameters in this case.
# Gradient descent (finite-differenced) can be used to find the actual minimizer:

nJ = 30
nD = 35
loss_landscape = zeros(Float64,nJ,nD)
Js = range(-2,0,length=nJ)
Ds = range(8,12,length=nD)
for (ij,J) in enumerate(Js)
  for (id,D) in enumerate(Ds)
    loss_landscape[ij,id] = get_loss([J,D])
  end
end

fig = Figure()
ax = Axis(fig[1,1],xlabel = "J [meV]", ylabel = "D [meV]")
contourf!(ax,Js,Ds,loss_landscape)

x0 = [-2,9.5]
opt_result = optimize(get_loss,x0,method=GradientDescent(alphaguess=1e-3),store_trace=true,extended_trace = true,time_limit=10.)
lines!(ax,Point2f.(Optim.x_trace(opt_result)))
scatter!(ax,-1,10)
fig

# The fit can be verified by plotting the LSWT band structure over top of the experiment data:

bcs = axes_bincenters(SIMULATED_EXPERIMENT_HISTOGRAM_PARAMS)
f = Figure()#hide
ax = Axis(f[1,1]; xlabel="Q [R.L.U.]", ylabel="Energy (meV)")#hide
heatmap!(ax,bcs[1],bcs[4],log10.(SIMULATED_EXPERIMENT_DATA), colormap = :deepsea)
f#hide


J_trial, D_trial = opt_result.minimizer
set_exchange!(spin_one_chain,J_trial,nearest_neighbor_right)#hide
set_exchange!(spin_one_chain,J_trial,nearest_neighbor_left)#hide

set_onsite_coupling!(spin_one_chain, -D_trial*Sz^2, 1)#hide
one_site_system = reshape_geometry(spin_one_chain,[1 0 0; 0 1 0; 0 0 1])#hide

langevin.kT = 0.#hide
nStep = 1_000#hide
for _ in 1:nStep#hide
    step!(one_site_system, langevin)#hide
end#hide

swt = SpinWaveTheory(one_site_system)#hide
params = SIMULATED_EXPERIMENT_HISTOGRAM_PARAMS

path = [[q,0,0] for q in bcs[1]]
disp, intensity = intensities_bands(swt, path)

for i in axes(disp)[2]
    lines!(ax, bcs[1], disp[:,i]; color=intensity[:,i], colormap = :turbo,linewidth = 5,colorrange = (0.,1.))
end
Colorbar(f[1,2],colormap = :turbo, limits = (0.,1.))
Colorbar(f[1,3],colormap = :deepsea, limits = (0.,1.))
f



