# # Case Study: FeI$_{2}$
# 
# FeI$_2$ is an effective spin-1 material with strong single-ion anisotropy.
# Quadrupolar fluctuations give rise to a single-ion bound state that cannot be
# described by a dipole-only model. This tutorial illustrates how to use the
# linear spin wave theory of SU(3) coherent states (i.e. 2-flavor bosons) to
# model the magnetic behavior in FeI$_2$. The original study was performed in
# [Bai et al., Nature Physics 17, 467–472
# (2021)](https://doi.org/10.1038/s41567-020-01110-1).
#
# ```@raw html
# <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/docs/src/assets/FeI2_crystal.jpg" style="float: left;" width="400">
# ```
#
# The Fe atoms are arranged in stacked triangular layers. The effective spin
# interactions include various anisotropic exchange interactions, and a strong
# single-ion anisotropy:
# 
# ```math
# \mathcal{H}=\sum_{(i,j)} J^{\alpha\beta}_{ij} S^{\alpha}_i S^{\beta}_j - D\sum_i \left(S^z\right)^2
# ```
# 
# We will formulate this Hamiltonian in Sunny and then calculate its dynamic
# structure factor.
#
# Begin by importing `Sunny` and `GLMakie`, a plotting package.

using Sunny, GLMakie
#nb Sunny.offline_viewers()  # Inject Javascript code for additional plotting capabilities 

# If you see an error `Package <X> not found in current path`, add the package
# by typing `] add <X>` at the Julia prompt.

# ## Crystals and symmetry analysis
# A [`Crystal`](@ref) describes the crystallographic unit cell and will usually
# be loaded from a `.cif` file. Here, we instead build a crystal by listing all
# atoms and their types.

a = b = 4.05012  # Lattice constants for triangular lattice
c = 6.75214      # Spacing in the z-direction

latvecs = lattice_vectors(a, b, c, 90, 90, 120) # A 3x3 matrix of lattice vectors that
                                                 ## define the conventional unit cell
positions = [[0,0,0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]  # Positions of atoms in fractions
                                                         ## of lattice vectors
types = ["Fe", "I", "I"]
FeI2 = Crystal(latvecs, positions; types)

# Observe that Sunny inferred the space group, 'P -3 m 1' (164) and labeled the
# atoms according to their point group symmetries.

#nb # An interactive viewer of the crystal and its bonds is available for Jupyter notebooks.
#nb view_crystal(FeI2, 8.0)

# Only the Fe atoms are magnetic, so we discard the I ions using
# [`subcrystal`](@ref).

cryst = subcrystal(FeI2, "Fe")

# Importantly, `cryst` retains the spacegroup symmetry of the full FeI$_2$
# crystal. This information will be used, for example, to propagate exchange
# interactions between symmetry-equivalent bonds.

# ## Spin systems
# To simulate a system of many spins, construct a [`System`](@ref).

sys = System(cryst, (4,4,4), [SpinInfo(1, S=1, g=2)], :SUN, seed=2)

# The system includes $4×4×4$ unit cells, i.e. 64 Fe atoms, each with spin $S=1$
# and a $g$-factor of 2. Quantum mechanically, spin $S=1$ involves a
# superposition of $2S+1=3$ distinct angular momentum states. In `:SUN` mode,
# this superposition will be modeled explicitly using the formalism of SU(3)
# coherent states, which captures both dipolar and quadrupolar fluctuations. For
# the more traditional dipole dynamics, use `:dipole` mode instead.

# ## Interactions and anisotropies

# ### Symmetry analysis
# The next step is to add interactions to the system. The command
# [`print_symmetry_table`](@ref) shows all symmetry-allowed interactions up to a
# cutoff distance.

print_symmetry_table(cryst, 8.0)

# The allowed $g$-tensor is expressed as a 3×3 matrix in the free coefficients
# `A`, `B`, ... The allowed single-ion anisotropy is expressed as a linear
# combination of Stevens operators. The latter correspond to polynomials of the
# spin operators, as we will describe below.
# 
# The allowed exchange interactions are given as a 3×3 matrix for representative
# bonds. The notation `Bond(i, j, n)` indicates a bond between atom indices `i`
# and `j`, with cell offset `n`. In the general case, it will be necessary to
# associate atom indices with their positions in the unit cell; these can be
# viewed with `display(cryst)`. Note that the order of the pair $(i, j)$ is
# significant if the exchange tensor contains antisymmetric
# Dzyaloshinskii–Moriya (DM) interactions.
# 
# In the case of FeI$_2$, `Bond(1, 1, [1,0,0])` is one of the 6 nearest-neighbor
# Fe-Fe bonds on a triangular lattice layer, and `Bond(1, 1, [0,0,1])` is an
# Fe-Fe bond between layers. 

# ### Assigning interactions and anisotropies

# The function [`set_exchange!`](@ref) assigns an exchange interaction to a
# bond, and will propagate the interaction to all symmetry-equivalent bonds in
# the unit cell. The FeI$_2$ interactions below follow [Bai et al](https://doi.org/10.1038/s41567-020-01110-1).

J1pm   = -0.236 
J1pmpm = -0.161
J1zpm  = -0.261
J2pm   = 0.026
J3pm   = 0.166
J′0pm  = 0.037
J′1pm  = 0.013
J′2apm = 0.068

J1zz   = -0.236
J2zz   = 0.113
J3zz   = 0.211
J′0zz  = -0.036
J′1zz  = 0.051
J′2azz = 0.073

J1xx = J1pm + J1pmpm 
J1yy = J1pm - J1pmpm
J1yz = J1zpm

set_exchange!(sys, [J1xx   0.0    0.0;
                    0.0    J1yy   J1yz;
                    0.0    J1yz   J1zz], Bond(1,1,[1,0,0]))
set_exchange!(sys, [J2pm   0.0    0.0;
                    0.0    J2pm   0.0;
                    0.0    0.0    J2zz], Bond(1,1,[1,2,0]))
set_exchange!(sys, [J3pm   0.0    0.0;
                    0.0    J3pm   0.0;
                    0.0    0.0    J3zz], Bond(1,1,[2,0,0]))
set_exchange!(sys, [J′0pm  0.0    0.0;
                    0.0    J′0pm  0.0;
                    0.0    0.0    J′0zz], Bond(1,1,[0,0,1]))
set_exchange!(sys, [J′1pm  0.0    0.0;
                    0.0    J′1pm  0.0;
                    0.0    0.0    J′1zz], Bond(1,1,[1,0,1]))
set_exchange!(sys, [J′2apm 0.0    0.0;
                    0.0    J′2apm 0.0;
                    0.0    0.0    J′2azz], Bond(1,1,[1,2,1]))

# The function [`set_onsite_coupling!`](@ref) assigns a single-ion anisotropy
# operator. It can be constructed, e.g., from the matrices given by
# [`spin_operators`](@ref) or [`stevens_operators`](@ref). Here we construct an
# easy-axis anisotropy along the direction $\hat{z}$.

D = 2.165
S = spin_operators(sys, 1)
set_onsite_coupling!(sys, -D*S[3]^2, 1)

# Any anisotropy operator can be converted to a linear combination of Stevens
# operators with [`print_stevens_expansion`](@ref).

# # Calculating structure factor intensities

# In the remainder of this tutorial, we will examine Sunny's tools for
# calculating the dynamical structure factor using a generalization of linear
# spin wave theory that captures quasi-particle excitations that include coupled
# dipolar and quadrupolar fluctuations.
# 
# A review of the quantum theory is provided in https://arxiv.org/abs/1307.7731.
# The same theory can alternatively be viewed as the quantization of a classical
# dynamics on SU(_N_) coherent states, https://arxiv.org/abs/2106.14125.
# Although outside the scope of this tutorial, Sunny provides support for
# simulating the classical SU(_N_) spin dynamics via a [`Langevin`](@ref)
# equation that generalizes the stochastic Landau-Lifshitz equation,
# https://arxiv.org/abs/2209.01265. This Langevin equation can be used to study
# equilibrium thermal fluctuations of classical SU(_N_) coherent states, or
# their non-equilibrium dynamics. Examples of the latter include relaxation of a
# spin glass, or the driven-dissipative dynamics of topological magnetic
# textures. To learn more, see the other Sunny documentation examples.

# ## Finding the ground state

# Begin with a random configuration and use [`minimize_energy!`](@ref) to find a
# configuration of the SU(3) coherent states (i.e. spin dipoles and quadrupoles)
# that locally minimizes energy.

randomize_spins!(sys)
minimize_energy!(sys);

# The expected ground state for FeI$_2$ is an antiferrogmanetic striped phase
# with a period of four spins (two up, two down). Visualizing the result of
# optimization, however, may indicate the system got stuck in a local minimum
# with defects.

plot_spins(sys)

# A better understanding of the magnetic ordering can often be obtained by
# moving to Fourier space. The 'instant' structure factor $𝒮(𝐪)$ is an
# experimental observable. To investigate $𝒮(𝐪)$ as true 3D data, Sunny
# provides [`instant_correlations`](@ref) and related functions. Here, however,
# we will use the lighter weight function [`print_wrapped_intensities`](@ref) to
# get a quick understanding of the periodicities present in the spin
# configuration.

print_wrapped_intensities(sys)

# The precise output may vary with Sunny version due to, e.g., different
# floating point roundoff effects. Very likely, however, the result will be
# approximately consistent with the known zero-field energy-minimizing magnetic
# structure of FeI$_2$, which is single-$Q$. Mathematically, spontaneous
# symmetry breaking should select one of $±Q = [0, -1/4, 1/4]$, $[1/4, 0, 1/4]$,
# or $[-1/4,1/4,1/4]$, associated with the three-fold rotational symmetry of the
# crystal spacegroup. In practice, however, one will frequently encounter
# competing "domains" associated with the three possible orientations of the
# ground state.

# If the desired ground state is already known, as with FeI$_2$, it could be
# entered by hand using [`set_dipole!`](@ref). Alternatively, in the case of
# FeI$_2$, we could repeatedly employ the above randomization and minimization
# procedure until a defect-free configuration is found. Some systems will have
# more complicated ground states, which can be much more challenging to find.
# For this, Sunny provides experimental support for powerful simulated annealing
# via [parallel tempering](https://en.wikipedia.org/wiki/Parallel_tempering),
# but that is outside the scope of this tutorial.

# Here, let's break the three-fold symmetry of FeI$_2$ by hand. Given one or
# more desired $Q$ modes, Sunny can suggest a magnetic supercell with
# appropriate periodicity. Let's arbitrarily select one of the three possible
# ordering wavevectors, $Q = [0, -1/4, 1/4]$. Sunny suggest a corresponding
# magnetic supercell in units of the crystal lattice vectors.

suggest_magnetic_supercell([[0, -1/4, 1/4]], sys.latsize)

# The function [`reshape_supercell`](@ref) allows an arbitrary reshaping of the
# system's supercell. We select the supercell appropriate to the broken-symmetry
# ground-state, which makes optimization much easier.

sys_min = reshape_supercell(sys, [1 0 0; 0 1 -2; 0 1 2])
randomize_spins!(sys_min)
minimize_energy!(sys_min)
plot_spins(sys_min; ghost_radius=3)

# ## Linear spin wave theory
#
# Now that we have found the ground state for a magnetic supercell, we can
# immediately proceed to perform zero-temperature calculations using linear spin
# wave theory. We begin by instantiating a `SpinWaveTheory` type using the
# supercell.

swt = SpinWaveTheory(sys_min)

# Select a sequence of wavevectors that will define a piecewise linear
# interpolation in reciprocal lattice units (RLU).

q_points = [[0,0,0], [1,0,0], [0,1,0], [1/2,0,0], [0,1,0], [0,0,0]];

# The function [`reciprocal_space_path`](@ref) will linearly sample a `path`
# between the provided $q$-points with a given `density`. The `xticks` return
# value provides labels for use in plotting.

density = 50
path, xticks = reciprocal_space_path(cryst, q_points, density);

# The [`dispersion`](@ref) function defines the quasiparticle excitation
# energies $ω_i(𝐪)$ for each point $𝐪$ along the reciprocal space path.

disp = dispersion(swt, path);

# In addition to the band energies $ω_i(𝐪)$, Sunny can calculate the inelastic
# neutron scattering intensity $I_i(𝐪)$ for each band $i$ according to an
# [`intensity_formula`](@ref). We choose to apply a polarization correction
# $(1 - 𝐪⊗𝐪)$ by setting the mode argument to `:perp`. Selecting
# `delta_function_kernel` specifies that we want the energy and intensity of
# each band individually.

formula = intensity_formula(swt, :perp; kernel=delta_function_kernel)

# The function [`intensities_bands`](@ref) uses linear spin wave theory to
# calculate both the dispersion and intensity data for the provided path.

disp, intensity = intensities_bands(swt, path; formula);

# These can be plotted in GLMakie.

fig = Figure()
ax = Axis(fig[1,1]; xlabel="𝐪", ylabel="Energy (meV)", xticks, xticklabelrotation=π/6)
ylims!(ax, 0.0, 7.5)
xlims!(ax, 1, size(disp, 1))
for i in axes(disp)[2]
    lines!(ax, 1:length(disp[:,i]), disp[:,i]; color=intensity[:,i])
end
fig

# To make comparisons with inelastic neutron scattering (INS) data, it is
# helpful to employ an empirical broadening kernel, e.g., a
# [`lorentzian`](@ref).

γ = 0.15 # width in meV
broadened_formula = intensity_formula(swt, :perp; kernel=lorentzian(γ))

# The [`intensities_broadened`](@ref) function requires an energy range in
# addition to the $𝐪$-space path.

energies = collect(0:0.01:10)  # 0 < ω < 10 (meV).
is1 = intensities_broadened(swt, path, energies, broadened_formula);

# A real FeI$_2$ sample will exhibit competing magnetic domains associated with
# spontaneous symmetry breaking of the 6-fold rotational symmetry of the
# triangular lattice. Note that the wavevectors $𝐪$ and $-𝐪$ are equivalent in
# the structure factor, which leaves three distinct domain orientations, which
# are related by 120° rotations about the $ẑ$-axis. Rather than rotating the
# spin configuration directly, on can rotate the $𝐪$-space path. Below, we use
# [`rotation_in_rlu`](@ref) to average the intensities over all three possible
# orientations.

R = rotation_in_rlu(cryst, [0, 0, 1], 2π/3)
is2 = intensities_broadened(swt, [R*q for q in path], energies, broadened_formula)
is3 = intensities_broadened(swt, [R*R*q for q in path], energies, broadened_formula)
is_averaged = (is1 + is2 + is3) / 3

fig = Figure()
ax = Axis(fig[1,1]; xlabel="(H,0,0)", ylabel="Energy (meV)", xticks, xticklabelrotation=π/6)
heatmap!(ax, 1:size(is_averaged, 1), energies, is_averaged)
fig

# This result can be directly compared to experimental neutron scattering data
# from [Bai et al.](https://doi.org/10.1038/s41567-020-01110-1)
# ```@raw html
# <img src="https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/docs/src/assets/FeI2_intensity.jpg">
# ```
#
# (The publication figure accidentally used a non-standard coordinate system to
# label the wave vectors.)
# 
# To get this agreement, the use of SU(3) coherent states is essential (in other
# words, we needed a theory of multi-flavored bosons). The lower band has large
# quadrupolar character, and arises from the strong easy-axis anisotropy of
# FeI$_2$. By setting `mode = :SUN`, the calculation captures this coupled
# dipole-quadrupole dynamics.
#
# An interesting exercise is to repeat the same study, but using `mode =
# :dipole` instead of `:SUN`. That alternative choice would constrain the
# coherent state dynamics to the space of dipoles only.

# The full dynamical spin structure factor (DSSF) can be retrieved as a $3×3$
# matrix with the [`dssf`](@ref) function, for a given path of $𝐪$-vectors.

disp, is = dssf(swt, path);

# The first output `disp` is identical to that obtained from `dispersion`. The
# second output `is` contains a list of $3×3$ matrix of intensities. For
# example, `is[q,n][2,3]` yields the $(ŷ,ẑ)$ component of the structure factor
# intensity for `nth` mode at the `q`th wavevector in the `path`.

# ## What's next?
#
# Sunny provides much more than just linear spin wave theory. Please explore the
# other examples in the documentation for more details. Sunny's capability to
# simulate the classical [`Langevin`](@ref) dynamics of SU(_N_) coherent states
# allows to incorporate finite-temperature spin fluctuations into the dynamical
# structure factor calculation, or to study systems under highly non-equilibrium
# conditions. Sunny also supports the construction of inhomogeneous systems
# (e.g., systems with quenched disorder) via [`to_inhomogeneous`](@ref).
