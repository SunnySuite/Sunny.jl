# # 4. Generalized spin dynamics of FeI₂ at finite *T*

using Sunny, LinearAlgebra, GLMakie

# In the [previous FeI₂ tutorial](@ref "1. Multi-flavor spin wave simulations of
# FeI₂ (Showcase)"), we used multi-flavor spin wave theory to calculate the
# dynamical structure factor. Here, we perform a similar calculation using a
# [generalized classical spin dynamics](https://arxiv.org/abs/2209.01265) that
# captures the coupled dynamics of spin dipoles and quadrupoles for
# configurations sampled at finite temperature.
#
# Compared to spin wave theory, simulations using classical dynamics will be
# slower and limited in ``𝐪``-space resolution. However, they make it is
# possible to study [temperature driven phase
# transitions](https://arxiv.org/abs/2310.19905). They may also be used to study
# out-of-equilibrium systems (e.g., relaxation of spin glasses), or systems with
# quenched inhomogeneities that require large simulation volumes.
#
# In this tutorial, we show how to study the finite temperature dynamics of FeI₂
# using the classical approach. It is important to stress that the estimation of
# ``S(𝐪,ω)`` with classical dynamics is fundamentally a Monte Carlo
# calculation: sample spin configurations are drawn from thermal equilibrium and
# used as initial conditions for generating dissipationless trajectories. The
# correlations of these trajectories are then averaged and used to calculate
# scattering intensities. It is therefore important to ensure that the initial
# spin configurations are sampled appropriately and that sufficient statistics
# are collected. We will demonstrate one approach here.
#
# As an overview, we will:
#
# 1. Identify the ground state.
# 2. Measure correlation data describing the excitations around that ground
#    state.
# 3. Use the correlation data to compute scattering intensities.
#
# To begin, please follow our [previous tutorial](@ref "1. Multi-flavor spin
# wave simulations of FeI₂ (Showcase)") to initialize a FeI₂ `sys` with lattice
# dimensions ``4×4×4``. 

a = b = 4.05012#hide 
c = 6.75214#hide
latvecs = lattice_vectors(a, b, c, 90, 90, 120)#hide
positions = [[0,0,0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]#hide
types = ["Fe", "I", "I"]#hide
FeI2 = Crystal(latvecs, positions; types)#hide
cryst = subcrystal(FeI2, "Fe")#hide
units = Units(:meV)#hide
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

# ## Finding a ground state

# As [previously observed](@ref "1. Multi-flavor spin wave simulations of FeI₂
# (Showcase)"), direct energy minimization is susceptible to trapping in a local
# energy minimum.

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[s[3] for s in sys.dipoles])

# Alternatively, one can search for the ordered state by sampling spin
# configurations from thermal equilibrium. Sunny supports this via a
# [`Langevin`](@ref) dynamics of SU(_N_) coherent states. This dynamics involves
# a dimensionless `damping` magnitude and target temperature `kT` for thermal
# fluctuations.

kT = 0.2  # Temperature in meV
langevin = Langevin(; damping=0.2, kT)

# Use [`suggest_timestep`](@ref) to select an integration timestep for the given
# error tolerance, e.g. `tol=1e-2`. The spin configuration in `sys` should
# ideally be relaxed into thermal equilibrium, but the current, energy-minimized
# configuration will also work reasonably well.

suggest_timestep(sys, langevin; tol=1e-2)
langevin.dt = 0.027;

# Sample spin configurations using Langevin dynamics. We have carefully selected
# a temperature of 0.2 eV that is below the ordering temperature, but large
# enough to that the dynamics can overcome local energy barriers and annihilate
# defects.

for _ in 1:10_000
    step!(sys, langevin)
end

# Calling [`suggest_timestep`](@ref) shows that thermalization has not
# substantially altered the suggested `dt`.

suggest_timestep(sys, langevin; tol=1e-2)

# Although thermal fluctuations are present, the correct antiferromagnetic order
# (2 up, 2 down) has been found.

plot_spins(sys; color=[s[3] for s in sys.dipoles])

# For other phases, it can be much harder to find thermal equilibrium, and more
# complicated sampling procedures may be necessary.

# ## Calculating Thermal-Averaged Correlations $⟨S^{αβ}(𝐪,ω)⟩$
#
# Our aim is to study the classical spin dynamics for states sampled in thermal
# equilibrium. To minimize finite size effects, and achieve sufficient momentum
# space resolution, we should significantly enlarge the system volume. The
# function [`resize_supercell`](@ref) takes new dimensions as multiples of the
# unit cell lattice vectors.

sys_large = resize_supercell(sys, (16,16,4)) # 16x16x4 copies of the original unit cell
plot_spins(sys_large; color=[s[3] for s in sys_large.dipoles])

# Now we will re-thermalize the system to a configuration just above the
# ordering temperature.

kT = 3.5 * units.K # 3.5 K ≈ 0.30 meV
langevin.kT = kT
for _ in 1:10_000
    step!(sys_large, langevin)
end

# With this increase in temperature, the suggested timestep has increased slightly.

suggest_timestep(sys_large, langevin; tol=1e-2)
langevin.dt = 0.040;

# The next step is to collect correlation data ``S^{αβ}`` into a
# [`SampledCorrelations`](@ref) object. This will involve sampling spin
# configurations from thermal equilibrium, and then integrating [an
# energy-conserving generalized classical spin
# dynamics](https://arxiv.org/abs/2204.07563) to collect Fourier-space
# information about normal modes. Quantization of these modes yields the
# magnons, and the associated dynamical spin-spin correlations can be compared
# with neutron scattering intensities ``S^{αβ}(q,ω)``. Because this is a
# real-space calculation, data is only available for discrete ``q`` modes (the
# resolution scales like inverse system size).
#
# The `SampledCorrelations` object requires specification of an integration
# timestep `dt`, an energy range for the intensity measurements, and a choice of
# pair correlation measurement. A rule of thumb is that the `dt` can be twice
# larger than the Langevin integration timestep.

dt = 2*langevin.dt
energies = range(0, 7.5, 120)
sc = SampledCorrelations(sys_large; dt, energies, measure=ssf_perp(sys_large))

# The function [`add_sample!`](@ref) will collect data by running a dynamical
# trajectory starting from the current system configuration. 

add_sample!(sc, sys_large)

# To collect additional data, it is required to re-sample the spin configuration
# from the thermal distribution. Statistical error is reduced by fully
# decorrelating the spin configurations between calls to `add_sample!`.

for _ in 1:2
    for _ in 1:1000               # Enough steps to decorrelate spins
        step!(sys_large, langevin)
    end
    add_sample!(sc, sys_large)
end

# Now, `sc` has more samples included:

sc

# ## Computing Scattering Intensities

# Extract the thermally-averaged correlation data ``⟨S^{αβ}(q,ω)⟩`` for a given
# set of q-points using [`intensities`](@ref). A classical-to-quantum rescaling
# of normal mode occupations will be performed according to the temperature
# `kT`.

qs = Sunny.QPoints([[0, 0, 0], [0.5, 0.5, 0.5]])
res = intensities(sc, qs; kT, energies) 

fig = lines(res.energies, res.data[:,1]; axis=(xlabel="meV", ylabel="Intensity"), label="(0,0,0)")
lines!(res.energies, res.data[:,2]; label="(π,π,π)")
axislegend()
fig

# Significant fluctuations indicate a large stochastic error. This data can be
# made smooth by iterating over more calls to `add_sample!` iterations.
#
# Next, we will measure intensities along the [`q_space_path`](@ref) that
# connects high symmetry points. Here we will also apply a [`FormFactor`](@ref)
# appropriate to Fe²⁺.

points = [[0,   0, 0],  # List of wave vectors that define a path
          [1,   0, 0],
          [0,   1, 0],
          [1/2, 0, 0],
          [0,   1, 0],
          [0,   0, 0]] 
qpath = q_space_path(cryst, points, 500)
formfactors = [FormFactor("Fe2"; g_lande=3/2)]
res = intensities(sc, qpath; kT, energies, formfactors)
plot_intensities(res; colormap=:viridis, colorrange=(0.0, 1.0))

# On can also view the intensity along a ``𝐪``-space slice at a fixed energy
# value.

axis1 = [1, -1/2, 0]
axis2 = [0, 1, 0]
grid = q_space_grid(cryst, axis1, range(-1.5, 1.5, 300), axis2, (-1.5, 1.5))
res = intensities(sc, grid; energies=[3.88], kT)
plot_intensities(res; colorrange=(0.0, 0.3))

# Significant ``𝐪``-space discretization is apparent. This is a consequence of
# the finite system size used to simulate the classical dynamics. The basic
# intensity data stored inside `SampledCorrelations` is only available on a
# discrete grid of ``𝐪``-points, with resolution that scales inversely to
# linear system size. Use [`available_wave_vectors`](@ref) to obtain this grid
# as a 3D array.

size(available_wave_vectors(sc))
