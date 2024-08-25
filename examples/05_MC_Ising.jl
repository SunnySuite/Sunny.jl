# # 5. Monte Carlo sampling of the Ising model
# 
# This tutorial illustrates simulation of the classical 2D Ising model.

using Sunny, GLMakie

# [`Crystal`](@ref) unit cell are always 3D. To model a square lattice, we
# create a tetragonal cell with one atom and an elongated lattice constant
# ``c``.

a = 1
latvecs = lattice_vectors(a, a, 10a, 90, 90, 90)
crystal = Crystal(latvecs, [[0, 0, 0]])

# Create a [`System`](@ref) of spin dipoles. Following the Ising convention, we
# will restrict the dipoles to ``±1`` along the global ``\hat{z}``-axis. Select
# ``g=-1`` so that the Zeeman coupling between external field ``𝐁`` and spin
# dipole ``𝐬`` is ``-𝐁⋅𝐬``. The initial supercell size is ``L×L``.

L = 128
sys = System(crystal, [SpinInfo(1, S=1, g=-1)], :dipole; latsize=(L, L, 1), seed=0)
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
# only propose pure spin flips, ``𝐬 \rightarrow -𝐬``.
nsweeps = 4000
sampler = LocalSampler(kT=Tc, propose=propose_flip)
for i in 1:nsweeps
    step!(sys, sampler)
end

# Plot the Ising spins by extracting the ``z``-component of the dipoles
heatmap(reshape([s[3] for s in sys.dipoles], (L, L)))
