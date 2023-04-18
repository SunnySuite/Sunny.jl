# > ![](https://raw.githubusercontent.com/SunnySuite/Sunny.jl/main/assets/sunny_logo.jpg)
# _This is a
#  [tutorial](https://github.com/SunnySuite/SunnyTutorials/)
#  for the [Sunny](https://github.com/SunnySuite/Sunny.jl/) package, which
#  enables dynamical simulations of ordered and thermally disordered spins with
#  dipole and higher order moments._
#
# # Spin Dynamics of the Heisenberg pyrochlore antiferromagnet and applications to MgCr2O4
# **Author**: Martin Mourigal <br>
# **Date**: September 9, 2022 (Updated by David Dahlbom on April 17, 2023) <br>
#
# 
# In this tutorial, we will walk through a first example in Sunny and calculate
# the spin dynamical properties of the Heisenberg pyrochlore antiferromagnet and
# apply this knowledge to MgCr2O4 and ZnCr2O4, which are known to approximate
# this model. Relevant publications include:
# 
# [1] P. H. Conlon and J. T. Chalker, Phys. Rev. Lett. 102, 237206 (2009)
# (https://doi.org/10.1103/PhysRevLett.102.237206)
# 
# [2] P. H. Conlon and J. T. Chalker, Phys. Rev. B 81, 224413 (2010)
# (https://doi.org/10.1103/PhysRevB.81.224413)
# 
# [3] X. Bai, J. A. M. Paddison, _et al._ Phys. Rev. Lett. 122, 097201 (2019)
# (https://doi.org/10.1103/PhysRevLett.122.097201)
#
# ## Setting up Julia
# To run the examples in the tutorial, you will need a working installation of
# the Julia programming language and the Sunny package. Some useful references
# for getting started are:
#
# https://github.com/SunnySuite/Sunny.jl/wiki/Getting-started-with-Julia 
# 
# https://sunnysuite.github.io/Sunny.jl/dev/
#
# We will begin by loading the relevant packages.

using Sunny # The main package
using GLMakie # Plotting package
#nb Sunny.offline_viewers(); # Enable Javascript plotting
#nb Makie.inline!(true); # Embed Makie plots in notebook
#md Makie.inline!(true); # Embed Makie plots in markdown document
cif_path = joinpath("..", "Sunny.jl", "examples", "longer_examples");
 
# ## Setting up the crystal structure
# 
# Before specifying the interactions of our system, we first must set up the
# crystal. We will begin by specifying the pyrochlore lattice (illustrated
# below) in the manner that is typical of theorists.

# ![](https://www.oist.jp/sites/default/files/photos/20200110-diagram-of-pyrochlore-lattice.jpg)
# 
# _Picture Credits: Theory of Quantum Matter Unit, OIST_

# ### "Theorist" Method

## === Define the crystal structure of a B-site pyrochlore lattice with only B-site atoms === 
lat_vecs = lattice_vectors(8.3342, 8.3342, 8.3342, 90.0, 90.0, 90.0);  
bas_vecs = [  [0.87500000, 0.62500000, 0.37500000],
              [0.62500000, 0.12500000, 0.62500000],
              [0.87500000, 0.87500000, 0.12500000],  
              [0.62500000, 0.87500000, 0.37500000],  
              [0.87500000, 0.12500000, 0.87500000],  
              [0.62500000, 0.62500000, 0.12500000],  
              [0.87500000, 0.37500000, 0.62500000],  
              [0.62500000, 0.37500000, 0.87500000],  
              [0.37500000, 0.62500000, 0.87500000],  
              [0.12500000, 0.12500000, 0.12500000],  
              [0.37500000, 0.87500000, 0.62500000],  
              [0.12500000, 0.87500000, 0.87500000],  
              [0.37500000, 0.12500000, 0.37500000],  
              [0.12500000, 0.62500000, 0.62500000],  
              [0.37500000, 0.37500000, 0.12500000],  
              [0.12500000, 0.37500000, 0.37500000]];
bas_typs = ["B","B","B","B", "B","B","B","B", "B","B","B","B", "B","B","B","B"]
xtal_pyro   = Crystal(lat_vecs, bas_vecs; types=bas_typs) # We will call this crystal the Theoretical Pyrochlore

## === Return crystalographic information gathered by Sunny about the crystal structure === 
view_crystal(xtal_pyro, 3.2) # Sunny can draw the crystal structure.
display(xtal_pyro) # Sunny can identify the space-group and point-group for this system. Here Fd-3m

# ### "Experimentalist" Method #1 (Incorrect)
# A real crystal is more complicated than this, however, and we will now
# construct the system using the actual CIF file of MgCr2O4 from ICSD. This can
# be done by directly importing the data from a CIF file by hand, but this can
# be dangerous, as shown below.

## === Define the crystal structure of MgCr2O4 by copying the info from a .cif file === 
lat_vecs     = lattice_vectors(8.3342, 8.3342, 8.3342, 90.0, 90.0, 90.0)
bas_vecs     = [ [0.12500, 0.12500, 0.12500],
                 [0.50000, 0.50000, 0.50000],
                 [0.26070, 0.26070, 0.26070]]
bas_typs     = ["Mg","Cr","O"]
xtal_mgcro_1 = Crystal(lat_vecs, bas_vecs; types=bas_typs)

# Sunny returned a valid crystal, but it did get right space group for MgCr2O4.
# This can be fixed by modifying the input to include the space group and the
# setting.
# 
# ### "Experimentalist" Method #2 (Correct)

## === Define the crystal structure of MgCr2O4 by copying the info from a CIF file INCLUDING space group and setting === 
lat_vecs     = lattice_vectors(8.3342, 8.3342, 8.3342, 90.0, 90.0, 90.0)
bas_vecs     = [ [0.12500, 0.12500, 0.12500],
                 [0.50000, 0.50000, 0.50000],
                 [0.26070, 0.26070, 0.26070]]
bas_typs     = ["Mg","Cr","O"]
lat_spg      = 227 # Space Group Number 
lat_set      = "2" # Space Group setting
xtal_mgcro_2 = Crystal(lat_vecs, bas_vecs, lat_spg; types=bas_typs, setting=lat_set)

# This result is correct, but at this point we might as well import the
# CIF file directly, which we now proceed to do.

# ### "Experimentalist Method #3 (Correct -- if your CIF file is)

# To import a CIF file, simply give the path to `Crystal`. One may optionally
# specify a precision parameter to apply to the symmetry analysis.

## === Define the crystal structure of MgCr2O4 by importing a cif file === 
xtal_mgcro_3 = Crystal(joinpath(cif_path, "MgCr2O4_160953_2009.cif"); symprec=0.001)

# Finally, we wish to restrict attention to the magnetic atoms in the unit cell
# while maintaining symmetry information for the full crystal, which is required
# to determine the correct exchange and g-factor anisotropies. This can be
# achieved with the `subcrystal` function.

## === Define the crystal structure of MgCr2O4 by importing a cif file === 
xtal_mgcro = subcrystal(xtal_mgcro_2,"Cr")

# ## Making a `System` and assigning interactions 
# ### Making a `System`
# Before assigning and interactions, we first have to set up a `System` using
# our `Crystal`.

dims = (6, 6, 6)  # Supercell dimensions 
spininfos = [SpinInfo(1; S=3/2)]  # Specify spin information, note that all sites are symmetry equivalent 
sys_pyro = System(xtal_pyro, dims, spininfos, :dipole)    # Make a system in dipole (Landau-Lifshitz) mode on pyrochlore lattice
sys_mgcro = System(xtal_mgcro, dims, spininfos, :dipole); # Same on MgCr2O4 crystal

# To understand what interactions we may assign to this system, we have to
# understand the the symmetry properties of the crystal, which we turn to next.
#
# ### Symmetry analysis for exchange and single-ion anisotropies
# 
# `print_symmetry_table` reports all the allowed exchange interactions,
# single-site anisotropies, and g-factors on our crystal. It takes a `Cyrstal`
# and a distance. We'll look at both the "theoriest's" pyrochlore lattice,

display("======== Pyrochlore Lattice Exchange Interactions to Third Neighbor ========")
print_symmetry_table(xtal_pyro, 5.9) 

# and for the the MgCrO4 crystal,

display("======== MgCr2O4 Exchange Interactions to Third Neighbor ========")
print_symmetry_table(xtal_mgcro, 6.0) 


# Note that the exchange anisotropies allowed on the the pyrochlore lattice are
# slightly different from those on the MgCr2O4 cyrstal, due to the role of Mg
# and O in the bonds. Also note that Sunny has correctly identified the two
# inequivalent bonds 3a and 3b having the same length. A question may arises to
# know which bond is J3a and which is J3b, let's plot the structure.
#
#nb view_crystal(xtal_mgcro, 5.9) # Sunny can draw the crystal structure. 
#
# The above plot shows that the second interaction (cyan color) with the
# distance of 5.89Å is in fact the one hopping through a chromium site, meaning
# it is J3a! We will need to be careful with that later.
#
# ### Building the exchange interactions for our system
#
# We begin by setting the scale of our exchange interactions on each bond. 
## === Define Values of Exchange Interactions ===
J1      = 3.27  # value of J1 in meV from Bai's PRL paper
J_pyro  = [1.00,0.000,0.000,0.000]*J1    # pure nearest neighbor pyrochlore
J_mgcro = [1.00,0.0815,0.1050,0.085]*J1; # further neighbor pyrochlore relevant for MgCr2O4
## val_J_mgcro = [1.00,0.000,0.025,0.025]*val_J1; # this is a funny setting from Conlon-Chalker

# These values are then assigned to their corresponding bonds with `set_exchange!`.

## === Assign exchange interactions to pyrochlore system ===
set_exchange!(sys_pyro, J_pyro[1], Bond(1, 3, [0,0,0])) # J1
set_exchange!(sys_pyro, J_pyro[2], Bond(1, 2, [0,0,0])) # J2
set_exchange!(sys_pyro, J_pyro[3], Bond(2, 6, [0,0,0])) # J3a
set_exchange!(sys_pyro, J_pyro[4], Bond(1, 5, [0,0,0])) # J3b

## === Assign exchange interactions to MgCr2O4 system ===
set_exchange!(sys_mgcro, J_mgcro[1], Bond(1, 2, [0,0,0]))  # J1
set_exchange!(sys_mgcro, J_mgcro[2], Bond(1, 7, [0,0,0]))  # J2
set_exchange!(sys_mgcro, J_mgcro[3], Bond(1, 3, [0,0,0]))  # J3a -- Careful here!  
set_exchange!(sys_mgcro, J_mgcro[4], Bond(1, 3, [1,0,0])); # J3b -- And here!

# We will not be assigning any single-ion anisotropies, so we have finished
# specifying our models. For good measure, we will set both systems to a random
# (infinite temperature) initial condition.

randomize_spins!(sys_pyro)
randomize_spins!(sys_mgcro)

# ## Cooling our `System` amd calculating the instantaneous and dynamic structure factors at the final temperature ##
# 
# We begin by thermalizing our system at a particular temperature. We will
# accomplish this by running Langevin dynamics. To do this, we must set up a
# Langevin integrator.

Δt = 0.05  # Integration time step in meV^-1
λ  = 0.1   # Phenomenological damping parameter
kT = 1.8   # Desired temperature in meV
langevin = Langevin(Δt; λ, kT); # Construct integrator

# We can now thermalize our systems by running the integrator.
for _ in 1:2000
    step!(sys_pyro, langevin)
    step!(sys_mgcro, langevin)
end

# As a sanity check, we'll plot the real-space spin configurations of both
# systems after themalization. First the pyrochlore,

plot_spins(sys_pyro; arrowlength=0.5, linewidth=0.2, arrowsize=0.5)

# and then the MgCr2O4,

plot_spins(sys_mgcro; arrowlength=0.5, linewidth=0.2, arrowsize=0.5)

# ## Instantaneous Structure Factor
# Next we can examine the instantaneous structure factor.

isf_pyro  = InstantStructureFactor(sys_pyro)
isf_mgcro = InstantStructureFactor(sys_mgcro);

# This generates a single sample. Let's add 10 more.

for _ in 1:10
    ## Run dynamics to decorrelate
    for _ in 1:500
        step!(sys_pyro, langevin)
        step!(sys_mgcro, langevin)
    end
    ## Add samples
    add_sample!(isf_pyro, sys_pyro)
    add_sample!(isf_mgcro, sys_mgcro)
end

# To retrieve the intensities, we call `intensities` on an array of wave vectors.
qvals = -4.0:0.025:4.0
qs = [(qa, qb, 0) for qa in qvals, qb in qvals]      # Wave vectors to query

Sq_pyro  = instant_intensities(isf_pyro, qs, :perp)  # :perp applies polarization factor 
Sq_mgcro = instant_intensities(isf_mgcro, qs, :perp);

# Finally we can plot the results.

fig = Figure(; resolution=(1200,500))
ax_pyro  = Axis(fig[1,1]; aspect=true, title="Pyrochlore")
ax_mgcro = Axis(fig[1,3]; aspect=true, title="MgCr2O4")
hm = heatmap!(ax_pyro, qvals, qvals, Sq_pyro)
Colorbar(fig[1,2], hm)
hm = heatmap!(ax_mgcro, qvals, qvals, Sq_mgcro)
Colorbar(fig[1,4], hm)
fig

# ## Dynamical Structure Factor
# We can also calculate the dynamical structure factor. 

dsf_pyro  = DynamicStructureFactor(sys_pyro;  Δt, ωmax = 10.0, nω = 100)
dsf_mgcro = DynamicStructureFactor(sys_mgcro; Δt, ωmax = 10.0, nω = 100);

# Let's add a couple more samples. 

for _ in 1:2
    ## Run dynamics to decorrelate
    for _ in 1:500
        step!(sys_pyro, langevin)
        step!(sys_mgcro, langevin)
    end
    ## Add samples
    add_sample!(dsf_pyro, sys_pyro)
    add_sample!(dsf_mgcro, sys_mgcro)
end

# We can now examine the structure factor intensities along a path in momentum space. 

## Plot pyrochlore results at a few slices
fig = Figure(; resolution=(1200,900))
qbs = 0.0:0.5:1.5
for (i, qb) in enumerate(qbs)
    path, labels = connected_path([[-4.0, qb, 0.0],[4.0, qb, 0.0]], 40)  # Generate a path of wave
                                                                         ## vectors through the BZ
    Sqω_pyro  = intensities(dsf_pyro, path, :perp; )  # Temperature keyword enables intensity rescaling

    ax = Axis(fig[fldmod1(i,2)...]; xlabel = "q = (x, $qb, 0)", ylabel="E (meV)")
    heatmap!(ax, [p[1] for p in path], ωs(dsf_pyro), Sqω_pyro; colorrange=(0.0, 4.0))
end
fig

# And let's take a look at the same slices for MgCr2O4.
fig = Figure(; resolution=(1200,900))
qbs = 0.0:0.5:1.5
for (i, qb) in enumerate(qbs)
    path, labels = connected_path([[-4.0, qb, 0.0],[4.0, qb, 0.0]], 40)  # Generate a path of wave vectors through the BZ
    Sqω_mgcro  = intensities(dsf_mgcro, path, :perp; )  # Temperature keyword enables intensity rescaling

    ax = Axis(fig[fldmod1(i,2)...]; xlabel = "q = (x, $qb, 0)", ylabel="E (meV)")
    heatmap!(ax, [p[1] for p in path], ωs(dsf_mgcro), Sqω_mgcro; colorrange=(0.0, 4.0))
end
fig

# Finally, we note that the instant structure factor can be calculated from the
# dynamical structure factor. We simply call `instant_intensities` rather than
# `intensities`. An advantage of doing this (as opposed to using
# `InstantStructureFactor`) is that Sunny is able to apply a temperature- and
# energy-dependent intensities rescaling when working with the dynamic structure
# factor. 
qvals = -4.0:0.05:4.0
qs = [(qa, qb, 0) for qa in qvals, qb in qvals]      # Wave vectors to query

Sq_pyro  = instant_intensities(dsf_pyro, qs, :perp; kT)  
Sq_mgcro = instant_intensities(dsf_mgcro, qs, :perp; kT);

# Finally we can plot the results. It is useful to compare these to our earlier results
# using `InstantStructureFactor`.

fig = Figure(; resolution=(1200,500))
ax_pyro  = Axis(fig[1,1]; aspect=true, title="Pyrochlore")
ax_mgcro = Axis(fig[1,3]; aspect=true, title="MgCr2O4")
hm = heatmap!(ax_pyro, qvals, qvals, Sq_pyro)
Colorbar(fig[1,2], hm)
hm = heatmap!(ax_mgcro, qvals, qvals, Sq_mgcro)
Colorbar(fig[1,4], hm)
fig