# Demo of `MixedIntegrator`: a thermostatted "bath" region at the two ends of
# a spin chain, coupled to a "bulk" region that evolves under exact
# (symplectic) Hamiltonian dynamics. This is the kind of setup used when
# simulating S(q,ω): the bath injects/removes energy to hold the system near
# a target temperature, while the measured bulk region is free of the
# artificial damping and numerical drift that an explicit integrator would
# otherwise introduce into the dynamics.

using Sunny, Printf

### Build a simple 1D Heisenberg spin chain

cryst = Crystal(lattice_vectors(1, 1, 1, 90, 90, 90), [[0, 0, 0]])
L = 40
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(L, 1, 1), seed=0)
set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))

### Thermalize the whole system first with ordinary `Langevin`, so the
### starting configuration is a sensible thermal state rather than pure noise.

kT_target = 0.3
randomize_spins!(sys)
thermalizer = Langevin(0.02; damping=0.2, kT=kT_target)
for _ in 1:2000
    step!(sys, thermalizer)
end

### Designate the first and last 5 sites along the chain (indices 1:5 and
### 36:40) as a Langevin bath held at `kT_target`. Every other site gets zero
### damping and zero kT, so `MixedIntegrator` will integrate it with the
### exactly energy-conserving `ImplicitMidpoint` method instead.

is_bath = [x <= 5 || x > L-5 for x in 1:L, _ in 1:1, _ in 1:1, _ in 1:1]
damping = map(b -> b ? 0.2 : 0.0, is_bath)
kT      = map(b -> b ? kT_target : 0.0, is_bath)

thermal = Langevin(0.02; damping, kT)
mixed = Sunny.MixedIntegrator(thermal)
println(mixed)
println("Bath sites: ", findall(vec(is_bath)))
println("Bulk sites (exact Hamiltonian dynamics): ", findall(vec(.!is_bath)))

### Run the mixed dynamics for a while and report energy statistics. The bath
### continually exchanges energy with the bulk, so total energy need not stay
### fixed -- but running the bulk's `ImplicitMidpoint` half-steps is what
### prevents it from accumulating its own artificial numerical drift on top
### of the physical bath coupling.

nsteps = 4000
Es = zeros(nsteps)
for i in 1:nsteps
    step!(sys, mixed)
    Es[i] = Sunny.energy_per_site(sys)
end

@printf("\nEnergy per site: start = %.4f, end = %.4f\n", Es[1], Es[end])
@printf("Energy per site over last 1000 steps: mean = %.4f, std = %.4f\n",
        sum(Es[end-999:end])/1000, sqrt(sum(abs2, Es[end-999:end] .- sum(Es[end-999:end])/1000)/1000))

### Sanity check: freeze the bath entirely (damping = kT = 0 everywhere) and
### run only the conservative half-steps on the bulk to confirm the energy
### of the full chain barely drifts -- i.e., the `ImplicitMidpoint` portion
### of `MixedIntegrator` is doing its job.

sys_check = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(L, 1, 1), seed=0)
set_exchange!(sys_check, -1.0, Bond(1, 1, [1, 0, 0]))
sys_check.dipoles .= sys.dipoles
cons = Sunny.ImplicitMidpoint(0.02)
E0 = Sunny.energy_per_site(sys_check)
for _ in 1:500
    step!(sys_check, cons; mask=trues(size(sys.dipoles)))
end
E1 = Sunny.energy_per_site(sys_check)
@printf("\nPure ImplicitMidpoint drift over 500 steps (reference): %.3e\n", abs(E1 - E0))
