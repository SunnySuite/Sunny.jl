""" An example of using the FastDipole library to compute SU(2)
     spin dynamics on a triangular lattice under the Hamiltonian:

    H = ∑ᵢⱼ Jᵢⱼ Sᵢ ⋅ Sⱼ + B ∑ᵢ Sᵢ

    where Jᵢⱼ = 2 for nearest-neighbor i, j, and 0 otherwise,
    and B = 1.0.
"""

using FastDipole

# Specify lattice
const lattice = Lattice(
    [ 1.0  0.5 ;   # Matrix of lattice vectors
      0.0 √3/2 ],
    [[0.0, 0.0]],  # List of basis vectors (here, just [0.0, 0.0])
    [5, 5]         # Make a size 5x5 lattice
)

# Specify Hamiltonian -- interface can / should be modified to integrate
#   with David's code
const J = 2.0
const field = ExternalField([0.0, 0.0, 1.0])
const pair_int = PairInteraction(J, 1, nothing, lattice)
const interactions = [pair_int, field]

# Instantiate the system
system = SpinSystem(lattice, interactions)
rand!(system)

# ==== Dynamics! ====

const NITERS = 10000
const Δt     = 0.001
energies = Vector{Float64}()

integrator = HeunP(system)
for it in 1:NITERS
    evolve!(integrator, Δt)
    # Compute the energy
    push!(energies, energy(system))
end

# Check that the energy hasn't fluctuated much
ΔE = maximum(energies) - minimum(energies)
println("After $(NITERS) iterations with Δt = $(Δt), observed ΔE = $(ΔE).")

println("Animating further integration...")
live_integration(system, 100, Δt)

