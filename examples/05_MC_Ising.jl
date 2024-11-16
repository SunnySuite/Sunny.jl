# # 5. Monte Carlo sampling of the Ising model
# 
# This tutorial illustrates simulation of the classical 2D Ising model.

using Sunny, GLMakie

# To model the 2D square lattice, create an elongated tetragonal cell with one
# atom.

latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
crystal = Crystal(latvecs, [[0, 0, 0]])

# Create a [`System`](@ref) of spin dipoles. Following the Ising convention, we
# will restrict the dipoles to ``±1`` along the global ``\hat{z}``-axis. Select
# ``g=-1`` so that the Zeeman coupling between external field ``𝐁`` and spin
# dipole ``𝐒`` is ``-𝐁⋅𝐒``. The system size is 128×128.

L = 128
sys = System(crystal, [1 => Moment(s=1, g=-1)], :dipole; dims=(L, L, 1))
polarize_spins!(sys, (0, 0, 1))

# Use [`set_exchange!`](@ref) to include a ferromagnetic Heisenberg interaction
# along nearest-neighbor bonds. The [`Bond`](@ref) below connects two spins
# displaced by the lattice vector ``𝐚₁``. This interaction will be propagated
# to all nearest-neighbors bonds in the system, consistent with the symmetries
# of the square lattice.

set_exchange!(sys, -1.0, Bond(1, 1, (1, 0, 0)))

# If an external field is desired, it can be set using [`set_field!`](@ref).

B = 0
set_field!(sys, (0, 0, B))

# The critical temperature for the Ising model is known analytically.

Tc = 2/log(1+√2)

# Use a [`LocalSampler`](@ref) to perform `nsweeps` Monte Carlo sweeps. A sweep
# consists of, on average, one trial update per spin in the system. Each
# proposed update is accepted or rejected according to the Metropolis acceptance
# probability. As its name suggests, the [`propose_flip`](@ref) function will
# only propose pure spin flips, ``𝐒 \rightarrow -𝐒``.

nsweeps = 4000
sampler = LocalSampler(kT=Tc, propose=propose_flip)
for i in 1:nsweeps
    step!(sys, sampler)
end

# Plot the Ising spins by extracting the ``z``-component of the dipoles

heatmap(reshape([S[3] for S in sys.dipoles], (L, L)))
