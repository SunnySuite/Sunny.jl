# # Classical Ising model
# 
# This tutorial illustrates how to simulate the classical 2D Ising model.

using Sunny, Plots

# Sunny expects a 3D [`Crystal`](@ref) unit cell. To model a square lattice, we
# create an orthogonal unit cell where the $z$-spacing is distinct from the $x$
# and $y$ spacing.
a = 1
latvecs = lattice_vectors(a,a,10a,90,90,90)
crystal = Crystal(latvecs, [[0,0,0]])

# Now create a system of spins with linear size `L = 128` in the $x$ and $y$
# directions, and only one layer in the $z$ direction. Following the Ising model
# convention, set the spin magnitude to $S = 1$. The option `:dipole` means that
# the system will store Heisenberg spins, as opposed to SU($N$) coherent states.
# By restricting these spins to the $z$-axis, they become Ising-like.
# 
# By default, Sunny uses physical units for the external field (Tesla) and for
# the magnetic moment of a spin ($S g μ_B$, where $g=2$ by default and $μ_B$ is
# the Bohr magneton). The options below effectively set $g = μ_B = 1$.
# Consequently, the Zeeman coupling of spin $s$ to external field $H$ has the
# dimensionless form $H s$.
L = 128
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=1)], :dipole, units=Units.theory, seed=0)
polarize_spins!(sys, (0,0,1))

# Select an arbitrary nearest-neighbor bond (here, two sites displaced by one
# lattice constant in the $x$-direction) and set a ferromagnetic exchange. This
# interaction will be symmetry-propagated to all nearest-neighbors bonds in the
# system.
set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

# If an external field is desired, it can be set as below.
H = 0
set_external_field!(sys, (0,0,H))

# The critical temperature for the Ising model is known analytically.
Tc = 2/log(1+√2)

# Perform `nsweeps` Monte Carlo sweeps. A sweep consists of, on average, one
# trial spin flip per site in the system, which is accepted or rejected using
# the Metropolis acceptance criterion.
nsweeps = 4000
sampler = LocalSampler(kT=Tc, propose=propose_flip)
for i in 1:nsweeps
    step!(sys, sampler)
end

# Plot the Ising spins by extracting the $z$-component of the dipoles
heatmap(reshape([s.z for s in sys.dipoles], (L,L)))
