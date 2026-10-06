using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.5"

latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
crystal = Crystal(latvecs, [[0, 0, 0]])

L = 128
sys = System(crystal, [1 => Moment(s=1, g=-1)], :dipole; dims=(L, L, 1))
polarize_spins!(sys, (0, 0, 1))

set_exchange!(sys, -1.0, Bond(1, 1, (1, 0, 0)))

B = 0
set_field!(sys, (0, 0, B))

Tc = 2/log(1+√2)

nsweeps = 4000
sampler = LocalSampler(kT=Tc, propose=propose_flip)
for i in 1:nsweeps
    step!(sys, sampler)
end

heatmap(reshape([S[3] for S in sys.dipoles], (L, L)))
