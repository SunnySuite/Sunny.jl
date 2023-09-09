# # Monte Carlo and Skyrmions
#
# This example considers a simple model on the triangular lattice that gives
# rise to skyrmions at finite temperature [Okubo et al, Phys. Rev. Lett. 108,
# 017206 (2012)](https://doi.org/10.1103/PhysRevLett.108.017206). The
# frustration in this system poses a challenge for Monte Carlo simulation, and
# provides an opportunity to demonstrate some of the advanced sampling
# algorithms provided by Sunny.

using Sunny, GLMakie

# The Okubo et al. model is defined on a triangular lattice.

a, c = 1.0, 10.0
latvecs = lattice_vectors(a, a, c, 90, 90, 120)
positions = [[0, 0, 0]]
cryst = Crystal(latvecs, positions)

# The model involves spin dipoles, and requires dimensionless ("theory") units
# to get the correct Zeeman coupling magnitude (Bohr magneton of one).

dims = (36,36,1)
sys = System(cryst, dims, [SpinInfo(1; S=1, g=1)], :dipole, units=Units.theory, seed=0)

# The J₁-J₃ model involves a ferromagnetic nearest neighbor interaction, no
# next-nearest neighbor interaction, and an antiferrogmanetic next-next-nearest
# neighbor interaction. The initial study employed `J₃/J₁ = -1` yielding a
# preferred wavevector with periodicity √3/2 lattice constants.

J₁ = -1
J₂ = 0.0
J₃ = -3J₁

set_exchange!(sys, J₁, Bond(1, 1, [1,0,0]))
set_exchange!(sys, J₂, Bond(1, 1, [1,2,0]))
set_exchange!(sys, J₃, Bond(1, 1, [2,0,0]))

# The phase diagram includes single-Q, double-Q, and triple-Q regions. The field
# and temperature conditions below are within the 3Q, Skyrmion crystal phase.

h = 2.0*J₃
kT = 0.3*J₃
set_external_field!(sys, [0, 0, h])


# We  will use Monte Carlo to generate spin configurations sampled from the
# Boltzmann equilibrium distribution. Two options provided by Sunny are
# [`LocalSampler`](@ref) and [`Langevin`](@ref). The local sampling method
# proposes random updates to individual spins, and accepts or rejects the
# proposal according to the [Metropolis
# algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm).
# This method is easier to use, and works well even for Ising-like spins.
#
# By contrast, the Langevin method samples spins using a continuous stochastic
# differential equation. The method begins with the physical Landau-Lifshitz
# spin dynamics, and adds phenomenological noise and damping terms to model
# coupling to a thermal bath. A [fluctuation-dissipation
# theorem](https://en.wikipedia.org/wiki/Fluctuation-dissipation_theorem)
# guarantees good samples for the target temperature `kT`. A tuneable,
# dimemnsionless parameter $λ$ controls the strength of coupling to the thermal
# bath (a good choice is $λ = 0.1$, while the the theoretical limit $λ → ∞$
# produces pure Brownian motion). The integration timestep `Δt` must be
# carefully selected based on the energy scales in the problem (this will be
# discussed in detail below). 

Δt = 0.02
λ = 0.1
langevin = Langevin(Δt; kT, λ)

# Starting from a random configuration, we may attempt to equilibrate into the
# 3Q, Skyrmion crystal phase.

randomize_spins!(sys)
for _ in 1:50_000
    step!(sys, langevin)
end

# Use custom plotting functions to observe the z component of spin

include(joinpath(pkgdir(Sunny), "examples", "extra", "plotting2d.jl"))
z_mean(s1, s2, s3) = (s1[3] + s2[3] + s3[3]) / 3
plot_triangular_plaquettes(z_mean, [sys.dipoles])

# Through tuning the ratio of $J₃/J₁$, it is possible to favor an arbitrary
# wavevector magnitude $Q$. 

x = 4
Q = 2π/x  # Wavelength of `x` lattice constants
J₃ = -J₁ / (8cos(Q/2)^2 - 4cos(Q/2)) # 0.2913...
set_exchange!(sys, J₃, Bond(1, 1, [2,0,0]))

h = 0.1*J₃
set_external_field!(sys, [0, 0, h])

langevin.Δt = 0.1 # Energy scales have decreased

randomize_spins!(sys)
for _ in 1:50_000
    step!(sys, langevin)
end

plot_triangular_plaquettes(z_mean, [sys.dipoles])

# TODO: Compare alternative sampling strategy.

sampler = LocalSampler(; kT, propose=propose_uniform)
randomize_spins!(sys)
for _ in 1:20_000
    step!(sys, sampler)
end
plot_spins(sys; arrowlength=1.0, linewidth=0.35, arrowsize=0.5)
