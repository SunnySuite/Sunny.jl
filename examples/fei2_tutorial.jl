# # Case Study: FeI$_{2}$
# 
# FeI$_{2}$ is an effective spin-1 material with strong single-ion anisotropy,
# making it an excellent candidate for treatment with SU(3) spin dynamics. In
# particular, one of the elementary excitations of the system can only be
# captured clasically with an SU(3) treatment. A magnon, clearly visible in the
# experimental data, would simply be absent if we were to employ traditional
# Landau-Lifshitz dynamics or SU(2) spin wave theory. Full details about the
# model can be found in reference [1].
# 
# The model contains a number of competing, anisotropic exchange interactions
# together with a strong single-ion anisotropy. Writing the exchange terms in
# the most general way, the Hamiltonian has the form:
# 
# ```math
# \mathcal{H}=\sum_{(i,j)} J^{\alpha\beta}_{ij} S^{\alpha}_i S^{\beta}_j - D\sum_i \left(S^z\right)^2
# ```
# 
# We will calculate a dynamic structure factor using this model. We begin by
# importing the required packages, starting with `Sunny`. We will also add
# `GLMakie`, a plotting package. If you see `Package X not found in current
# path`, you can install the package by entering `using Pkg; pkg"add X"` in the
# Julia REPL.

using Sunny, LinearAlgebra, GLMakie
#nb Sunny.offline_viewers()  # Inject Javascript code for additional plotting capabilities 


# ## Crystals and symmetry analysis
# One begins by constructing a [`Crystal`](@ref) that describes the
# crystallographic unit cell. Frequently, the crystal will be loaded from a
# `.cif` file. Below, we construct the crystal by providing a full list of all
# atoms in the conventional unit cell.

a = b = 4.05012  # Lattice constants for triangular lattice
c = 6.75214      # Spacing in the z-direction

lat_vecs = lattice_vectors(a, b, c, 90, 90, 120) # A 3x3 matrix of lattice vectors that
                                                 ## define the conventional unit cell
basis_vecs = [[0,0,0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]  # Positions of atoms in fractions
                                                          ## of lattice vectors
types = ["Fe", "I", "I"]
FeI2 = Crystal(lat_vecs, basis_vecs; types)

# Observe that Sunny indentified the correct space group, 'P -3 m 1' (164), in
# agreement with Ref. \[1\].

#nb # An interactive viewer of the crystal and its bonds is available for Jupyter notebooks.
#nb view_crystal(FeI2, 8.0)

# Only the Fe atoms are magnetic, so we discard the I ions using
# [`subcrystal`](@ref).

cryst = subcrystal(FeI2, "Fe")

# Importantly, `cryst` retains the spacegroup symmetry information for the full
# FeI2 crystal. This information will be used, for example, to propagate
# exchange interactions by symmetry.

# ## Spin systems
# The FeI2 unit cell contains only a single Fe atom. To simulate a system of
# many spins, construct a [`System`](@ref).

sys = System(cryst, (4,4,4), [SpinInfo(1,S=1)], :SUN, seed=0)

# This system includes $4×4×4$ unit cells, i.e. 64 spin moments, each with spin
# $S=1$. The default $g$-factor is 2, but this could be overriden with an
# additional argument to [`SpinInfo`](@ref). Recall that spin-1 has $N=2S+1=3$
# distinct angular momentum states. Because we selected `:SUN` mode, Sunny will
# simulate the dynamics of SU(3) coherent states, which includes both dipolar
# and quadrupolar fluctuations. For the more traditional dipole dynamics, use
# mode `:dipole` instead.

# ## Interactions and anisotropies

# ### Symmetry analysis
# The next step is to add interactions to the system. The command
# [`print_symmetry_table`](@ref) shows all symmetry-allowed interactions up to a
# distance cutoff.

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
# As an example, `Bond(1, 1, [1,0,0])` involves an offset of one unit cell in
# the direction of the first lattice vector. In the case of FeI2, this is one of
# the 6 nearest-neighbor Fe-Fe bonds on the triangular lattice.

# ### Assigning interactions and anisotropies

# The function [`set_exchange!`](@ref) assigns an exchange interaction to a
# bond, and will propagate the interaction to all symmetry-equivalent bonds in
# the unit cell. Below we define the FeI2 interactions following Ref. \[1\].

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

# The function [`set_anisotropy!`](@ref) assigns a single-ion anisotropy. It
# takes an abstract operator and an atom index. The operator may be a polynomial
# of spin operators or a linear combination of Stevens operators. Sunny provides
# special symbols for their construction: [`𝒮`](@ref) is a vector of the three
# spin operators and [`𝒪`](@ref) are the symbolic Stevens operators. Here we
# construct an easy-axis anisotropy.

D = 2.165
set_anisotropy!(sys, -D*𝒮[3]^2, 1)

# The function [`print_anisotropy_as_stevens`](@ref) will convert an anisotropy
# operator to a linear combination of Stevens operators.

# # Calculating a dynamical spin structure factor
# In the remainder of this tutorial, we will examine Sunny's tools for
# calculating structure factors using generalized SU(_N_) classical dynamics.
# This is a Monte Carlo calculation and will require the sampling of spin
# configurations from the Boltzmann distribution at a particular temperature.
# These samples are then used to generate dynamical trajectories that are
# analyzed to produce correlation information, i.e., a dynamical structure
# factor $\mathcal{S}^{\alpha\beta}(\mathbf{q},\omega)$.
# 
# ## Finding a ground state
# 
# The [`Langevin`](@ref) dynamics can be used to sample spin configurations in
# thermal equlibrium.

E0 = 2.165        # Largest energy scale in the Hamiltonian
Δt = 0.05/E0      # Safe choice for integration step size
λ = 0.1           # Magnitude of coupling to thermal bath
langevin = Langevin(Δt, 0, λ);

# Below we use a simulated annealing to find a low-energy configuration. We lower
# the temperature from `E0` to `0` over `nsteps` of the Langevin dynamics. 

randomize_spins!(sys)
nsteps = 10_000
for kT in range(E0, 0, nsteps)
    langevin.kT = kT
    step!(sys, langevin)
end

# The annealed configuration can be visualized with `plot_spins`.

plot_spins(sys; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# There are likely to be defects in the magnetic order because `nsteps = 10_000`
# is a relatively fast quench. Repeating the annealing procedure with, say,
# `nsteps = 100_000` would likely yield the correct magnetic order. Let's
# instead continue with the imperfect one.
#
# The zero-field energy-minimizing magnetic structure of FeI$_2$ is known to be
# single-$q$. If annealing were perfect, then spontaneous symmetry breaking
# would select one of $±q = [0, -1/4, 1/4]$, $[1/4, 0, 1/4]$, or
# $[-1/4,1/4,1/4]$. Let us check which modes appear in the static structure
# factor intensity (for each sublattice independently).

print_dominant_wavevectors(sys)

# Let's break the symmetry by hand. Given a list of $q$ modes, Sunny can suggest
# a magnetic supercell with appropriate periodicity. The result is printed in
# units of the crystal lattice vectors.

suggest_magnetic_supercell([[0, -1/4, 1/4]], sys.latsize)

# After reshaping the system to the suggested supercell volume, it becomes much
# easier to relax the spin configuration.

A = [1 0 0; 0 1 -2; 0 1 2]
sys_supercell = reshape_volume(sys, A)

randomize_spins!(sys_supercell)
langevin.kT = 0
for i in 1:nsteps
    step!(sys_supercell, langevin)
end

plot_spins(sys_supercell; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# We can now extend this magnetic supercell to a much larger simulation volume.

sys_large = reshape_volume(sys_supercell, diagm([16,16,4]))
plot_spins(sys_large; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)

# ## Calculating the structure factor
# Apply the Langevin dynamics to thermalize the system to 0.5K. Note that the
# number of time-steps required to decorrelate may vary significantly for
# different systems and thermodynamic conditions.

nsteps_decorr = round(Int, 2.0/Δt)     # System-dependent estimate
kT = 0.5 * meV_per_K            # 0.5K in units of meV
langevin.kT = kT

for _ in 1:5nsteps_decorr
    step!(sys_large, langevin)
end

# With our equilibrated system, we can begin to measure the
# [`DynamicStructureFactor`](@ref). Three keyword parameters are required to
# determine the ω information that will be calculated: an integration step size,
# the number of ωs to resolve, and the maximum ω to resolve. For the time step,
# twice the value used for the Langevin integrator is usually a good choice.

sf = DynamicStructureFactor(sys_large; Δt=2Δt, nω=120, ωmax=7.5);

# `sf` currently contains dynamical structure data generated from a single
# sample. Additional samples can be added by generating a new spin configuration
# and calling `add_sample!`:

for _ in 1:2
    for _ in 1:nsteps_decorr        # Generate a new sample spin configuration
        step!(sys_large, langevin)
    end
    add_sample!(sf, sys_large)      # Accumulate the sample into `sf`
end

# ## Accessing structure factor data 
# The basic function for accessing intensity data is [`intensities`](@ref),
# which, in addition to the structure factor data itself, takes a list of wave
# vectors and a mode parameter. The options for the mode parameter are `:trace`,
# `:perp` and `:full` which return, respectively, the trace, the unpolarized
# intensity, and the full set of matrix elements (correlations of spin
# components) at the specified wave vectors. For example, we can plot two
# single-$q$ slices as follows. 

qs = [[0, 0, 0], [0.5, 0.5, 0.5]]
is = intensities(sf, qs, :trace; kT)

fig = Figure()
ax = Axis(fig[1,1]; xlabel="meV", ylabel="Intensity")
l1 = lines!(ax, ωs(sf), is[1,:])
l2 = lines!(ax, ωs(sf), is[2,:])
Legend(fig[1,2], [l1, l2], ["(0,0,0)", "(π,π,π)"])
fig

# Note that we provided the optional keyword `kT` to `intensities` to enable
# Sunny to apply a classical-to-quantum rescaling of intensities. 
#
# Frequently we want to extract energy intensities along lines that connect
# special wave vectors. Sunny provides a function `connected_path` to makes this
# easy. The density of sample points can be tuned with a density argument.

points = [[0.0, 0.0, 0.0],  # List of wave vectors that define a path
          [1.0, 0.0, 0.0],
          [0.0, 1.0, 0.0],
          [0.5, 0.0, 0.0],
          [0.0, 1.0, 0.0],
          [0.0, 0.0, 0.0]] 
formfactors = [FormFactor(1, "Fe2"; g_lande=3/2)]  # Ion information for each site to
                                                   ## retrieve form factor correction parameters 
density = 40
path, markers = connected_path(points, density)

is = intensities(sf, path, :perp; 
    interpolation = :linear,       # Interpolate between available wave vectors
    kT,                            # Temperature for intensity correction
    formfactors,                   # Form factor information 
)

fig = Figure()
labels = ["($(p[1]),$(p[2]),$(p[3]))" for p in points]
ax = Axis(fig[1,1];
    ylabel = "meV",
    xticks = (markers, labels),
    xticklabelrotation=π/8,
    xticklabelsize=12,
)
heatmap!(ax, 1:size(is,1), ωs(sf), is; colorrange=(0.0, 0.5))
fig

# Often it is useful to plot cuts across multiple wave vectors but at a single
# energy. 

npoints = 60
qvals = range(-2.0, 2.0, length=npoints)
qs = [[a, b, 0.0] for a in qvals, b in qvals]
forfactors = [FormFactor(1, "Fe2")]

is = intensities(sf, qs, :perp;
    interpolation = :linear,
    kT,
    formfactors,
);

ωidx = 30
ω = ωs(sf)[ωidx]
fig = Figure()
ax = Axis(fig[1,1]; title="ω=$ω meV", aspect=true)
hidedecorations!(ax); hidespines!(ax)
hm = heatmap!(ax, is[:,:,ωidx])
Colorbar(fig[1,2], hm)
fig

# Note that Brillouin zones appear "skewed". This is a consequence
# of the fact that our reciprocal lattice vectors are not orthogonal. It is
# often useful to express our wave vectors in terms of an orthogonal basis,
# where each basis element is specified as a linear combination of reciprocal
# lattice vectors. For our crystal, with reciprocal vectors $a^*$, $b^*$ and
# $c^*$, we can define an orthogonal basis by taking $\hat{a}^* = 0.5(a^* +
# b^*)$, $\hat{b}^*=a^* - b^*$, and $\hat{c}^*=c^*$. Below, we map `qs` to
# wavevectors `ks` in the new coordinate system and get their intensities.

A = [0.5  1.0 0.0;
     0.5 -1.0 0.0;
     0.0  0.0 1.0]
ks = [A*q for q in qs]

is_ortho = intensities(sf, ks, :perp;
    interpolation = :linear,
    kT,
    formfactors,
);

fig = Figure()
ax = Axis(fig[1,1]; title="ω=$ω meV", aspect=true)
hidedecorations!(ax); hidespines!(ax)
hm = heatmap!(ax, is_ortho[:,:,ωidx])
Colorbar(fig[1,2], hm)
fig

# Finally, we note that static structure factor data can be obtained from a
# dynamic structure factor with `static_intensities`:

is_static = static_intensities(sf, ks, :perp;
    interpolation = :linear,
    kT,
    formfactors,
)

fig = Figure()
ax = Axis(fig[1,1]; title="Static Structure Factor", aspect=true)
hidedecorations!(ax); hidespines!(ax)
hm = heatmap!(ax, is_static)
Colorbar(fig[1,2], hm)
fig