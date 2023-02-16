using Sunny, Plots

latvecs = lattice_vectors(1,1,2,90,90,90)
crystal = Crystal(latvecs, [[0,0,0]])

L = 128

# Set g=1 and use 'theory' units (i.e., set μB=1) so that the external field H
# becomes dimensionless, following convention.
sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=1)], :dipole, units=Units.theory, seed=0)
polarize_spins!(sys, (0,0,1))

# Ferromagnetic nearest-neighbor exchange
set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))
# Optional applied field
H = 0
set_external_field!(sys, (0,0,H))


Tc = 2/log(1+√2) # Critical temperature 2.269...
sampler = LocalSampler(kT=Tc, propose=propose_flip)

for i in 1:2_000
    step!(sys, sampler)
end

# Plot the z-component of dipole
heatmap(reshape([s.z for s in sys.dipoles], (L,L)))
