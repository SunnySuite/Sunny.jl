using Sunny, GLMakie

a = 1
latvecs = lattice_vectors(a,a,10a,90,90,90)
crystal = Crystal(latvecs, [[0,0,0]])

L = 128
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=-1)], :dipole, seed=0)
polarize_spins!(sys, (0,0,1))

set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))

B = 0
set_field!(sys, (0, 0, B))

Tc = 2/log(1+√2)

nsweeps = 4000
sampler = LocalSampler(kT=Tc, propose=propose_flip)
for i in 1:nsweeps
    step!(sys, sampler)
end

heatmap(reshape([s.z for s in sys.dipoles], (L,L)))
