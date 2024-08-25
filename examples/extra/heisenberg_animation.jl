using Sunny, GLMakie

# Square lattice
latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]])

# Simple Heisenberg model
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(10, 10, 1), seed=1)
J = -1.0
set_exchange!(sys, J, Bond(1, 1, (1, 0, 0)))
randomize_spins!(sys)

fig = plot_spins(sys; colorfn=i->sys.dipoles[i][3], colorrange=(-1, 1), ndims=2)

dt = 0.1/abs(J)
integrator = Langevin(dt; damping=0.05, kT=0)

# View an animation in real time
for _ in 1:500
    for _ in 1:5
        step!(sys, integrator)
    end
    notify(fig)
    sleep(1/60)
end

# Save an animation to file
randomize_spins!(sys)
record(fig, "animation.mp4", 1:500; framerate=30) do _
    for _ in 1:5
        step!(sys, integrator)
    end
    notify(fig)
end
