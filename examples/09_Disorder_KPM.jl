# # 9. Disordered systems
#
# This example uses [`SpinWaveTheoryKPM`](@ref) to efficiently calculate the
# neutron scattering spectrum of a disordered triangular antiferromagnet. The
# model is inspired by YbMgGaO4, as studied in [Paddison et al, Nature Phys.,
# **13**, 117–122 (2017)](https://doi.org/10.1038/nphys3971) and [Zhu et al,
# Phys. Rev. Lett. **119**, 157201
# (2017)](https://doi.org/10.1103/PhysRevLett.119.157201). Disordered occupancy
# of non-magnetic Mg/Ga sites can be modeled as a stochastic distribution of
# exchange constants and ``g``-factors. Including this disorder introduces
# broadening of the spin wave spectrum.

using Sunny, GLMakie, LinearAlgebra
import Random

# Setup
i = 23
for i in 1:50
    @show i
    Random.seed!(i)

    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1))
    set_exchange!(sys, +1.0, Bond(1, 1, [1,0,0]))

    randomize_spins!(sys)
    minimize_energy!(sys; g_abstol=1e-10).iterations |> println
    # @show Sunny.minimize_energy2!(sys)

    sys_inhom = to_inhomogeneous(repeat_periodically(sys, (10, 10, 1)))

    for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1, 1, [1, 0, 0]))
        noise = randn()/3
        set_exchange_at!(sys_inhom, 1.0 + noise, site1, site2; offset)
    end

    # randomize_spins!(sys_inhom)
    minimize_energy!(sys_inhom; g_abstol=1e-10, maxiters=2_000).iterations |> println

    set_field!(sys_inhom, [0, 0, 7.5])
    randomize_spins!(sys_inhom)
    minimize_energy!(sys_inhom).iterations |> println

    for site in eachsite(sys_inhom)
        noise = randn()/6
        sys_inhom.gs[site] = [1 0 0; 0 1 0; 0 0 1+noise]
    end

    randomize_spins!(sys_inhom)
    minimize_energy!(sys_inhom).iterations |> println
    # @show Sunny.minimize_energy2!(sys_inhom; maxiters=5_000)
    # println(energy_per_site(sys_inhom))
end

minimize_energy!(sys_inhom; subiters=30_000, maxiters=30_000)
minimize_energy!(sys_inhom; maxiters=30_000)

s1 = reshape(copy(sys_inhom.dipoles), 30, 30)


energies = map(range(0, 1.1, 100)) do c
    @. sys_inhom.dipoles = 0.5 * normalize((1-c) * s0 + c * s1)
    energy_per_site(sys_inhom)
end
plot(energies)

heatmap(reshape([s[2] for s in s0], (30, 30)))
heatmap(reshape([s[2] for s in s1], (30, 30)))
heatmap(reshape([s[1] for s in (s1 - s0)], (30, 30)))
heatmap(reshape([s[2] for s in (s1 - s0)], (30, 30)))

polarize_spins!(sys_inhom, [0, 0, -1])
println(energy_per_site(sys_inhom))
