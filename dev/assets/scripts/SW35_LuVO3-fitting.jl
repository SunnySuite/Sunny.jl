using Sunny, GLMakie, LinearAlgebra
@assert pkgversion(Sunny) >= v"0.9.1"

a, b, c = 5.2821, 5.6144, 7.5283
latvecs = lattice_vectors(a, b, c, 90, 90, 90)
positions = [[1/2, 0, 0]]
cryst = Crystal(latvecs, positions, "P b n m")

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole_uncorrected)

set_exchange!(sys, 1.0, Bond(1, 2, (0, 0, 0)), :Jab => 0)
set_exchange!(sys, 1.0,  Bond(1, 3, (0, 0, 0)), :Jc => 0)
set_onsite_coupling!(sys, S -> -S[1]^2,  1, :Kxx => 0)
set_onsite_coupling!(sys, S -> -S[2]^2,  1, :Kyy => 0)

labels = [:Jab, :Jc, :Kxx, :Kyy]
skoulatos_fit = [5.95, 4.24, -0.48, -0.06] # Typo fixed
set_params!(sys, labels, skoulatos_fit)

minimize_energy!(sys)
plot_spins(sys; ghost_radius=4)

Q1 = [0, 1.0, 2.0]
Q2 = [0, 1.0, 3.0]
Q3 = [0, 1.5, 3.5]
path = q_space_path(cryst, [Q1, Q2, Q3], 500)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities_bands(swt, path)
plot_intensities(res; ylims=(0, 35), title="Previous work")

qs = [
    [0.0, 1.0, 2.0], [0.0, 1.0, 2.1], [0.0, 1.0, 2.2], [0.0, 1.0, 2.3],
    [0.0, 1.0, 2.4], [0.0, 1.0, 2.5], [0.0, 1.0, 2.6], [0.0, 1.0, 2.7],
    [0.0, 1.0, 2.8], [0.0, 1.0, 2.9], [0.0, 1.0, 3.0], [0.0, 1.1, 3.1],
    [0.0, 1.2, 3.2], [0.0, 1.3, 3.3], [0.0, 1.4, 3.4], [0.0, 1.5, 3.5]
]
Es = [
    [28.311], [28.111], [27.395], [26.279], [24.876], [22.758], [20.296],
    [17.405], [13.884], [5.697, 10.391], [3.693, 8.674], [13.027], [20.683],
    [27.368], [31.798], [33.431]
]

loss = make_loss_fn(sys, labels) do sys
    # minimize_energy!(sys) # Uncomment if the spin state should vary with params
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    res = intensities_bands(swt, qs)
    squared_error_bands(Es, res)
end

guess = [3, 3, -0.2, -0.1] # Guess for [Jab, Jc, Kxx, Kyy]
loss(guess)

import Optim

method = Optim.NelderMead()
options = Optim.Options(; g_tol=1e-8)
fit = Optim.optimize(loss, guess, method, options)
@assert isapprox(fit.minimizer, [6.1035, 3.9893, -0.6294, -0.0899]; rtol=1e-3)
fit.minimizer # [Jab, Jc, Kxx, Kyy]

U = uncertainty_matrix(loss, fit.minimizer)
@assert isapprox(sqrt.(diag(U) / 2), [0.2120, 0.1924, 0.0812, 0.0425]; rtol=1e-3)
sqrt.(diag(U) / 2) # [ΔJab, ΔJc, ΔKxx, ΔKyy]

set_params!(sys, labels, fit.minimizer)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities_bands(swt, path)
fig = plot_intensities(res; ylims=(0, 35), title="Updated fit")
qinds = find_qs_along_path(qs, path)
data_pairs = [(q, Eq) for (q, Eqs) in zip(qinds, Es) for Eq in Eqs]
plot!(fig[1, 1], data_pairs; color=:transparent, strokecolor=:magenta, strokewidth=3)
fig
