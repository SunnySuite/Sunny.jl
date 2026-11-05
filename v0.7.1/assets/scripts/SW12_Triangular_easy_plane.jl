using Sunny, GLMakie

latvecs = lattice_vectors(3, 3, 4, 90, 90, 120)
cryst = Crystal(latvecs, [[0, 0, 0]])

s = 3/2
J1 = +1.0
sys = System(cryst, [1 => Moment(; s, g=2)], :dipole; dims=(3, 3, 1))
set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))

undo_classical_to_quantum_rescaling = 1 / (1 - 1/2s)
D = 0.2 * undo_classical_to_quantum_rescaling
set_onsite_coupling!(sys, S -> D*S[3]^2, 1)
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; ndims=2)

qs = [[0, 0, 0], [1, 1, 0]]
path = q_space_path(cryst, qs, 400)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities_bands(swt, path)
plot_intensities(res)

measure = ssf_custom(sys) do q, ssf
    return real(ssf[3, 3])
end
swt = SpinWaveTheory(sys; measure)
res = intensities_bands(swt, path)
plot_intensities(res)

measure = ssf_custom((q, ssf) -> ssf, sys)
swt = SpinWaveTheory(sys; measure)
res = intensities_bands(swt, path)
res.data[7, 10]
