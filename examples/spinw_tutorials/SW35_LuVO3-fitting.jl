# # SW35 - LuVO₃ fitting
#
# This is a Sunny adaption of [SpinW Tutorial
# 35](https://spinw.org/tutorials/35tutorial) (no author attribution). This
# tutorial fits to the dispersion curves of LuVO₃ in its phase III. Data is
# obtained from Fig. 1c of [Skoulatos et al., Phys. Rev. B **91**, 161104(R)
# (2015)](https://doi.org/10.1103/PhysRevB.91.161104).
#
# Construct the [`Crystal`](@ref) using the custom ITA setting "P b n m" for
# spacegroup 62. The V³⁺ ions live on Wyckoff 4b.

using Sunny, GLMakie

units = Units(:meV, :angstrom)
a, b, c = 5.2821, 5.6144, 7.5283
latvecs = lattice_vectors(a, b, c, 90, 90, 90)
positions = [[1/2, 0, 0]]
cryst = Crystal(latvecs, positions, "P b n m")

# Construct the [`System`](@ref) with mode `:dipole_uncorrected` to avoid
# renormalization of the single-ion anisotropy. This is less
# quantum-mechanically accurate, but facilitates numerical comparison with
# previously fitted values.

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole_uncorrected)

# Exchanges ``J_{ab}`` and ``J_{c}`` are in the x̂-ŷ plane and the ẑ
# direction, respectively. The single-ion anisotropy term has the assumed form
# ``- \sum_α K_{αα} S_α^2``. Note that a shift of all ``K_{xx}, K_{yy}, K_{zz}``
# by some constant ``c`` amounts to a physically irrelevant shift ``- c |S|^2``
# of the spin Hamiltonian. Without loss of generality, one can select ``K_{zz} =
# 0``. With this convention, negative ``K_{xx}`` and ``K_{yy}`` amount to an
# easy-axis anisotropy in ``ẑ``. 

set_exchange!(sys, 1.0, Bond(1, 2, (0, 0, 0)), :Jab => 0)
set_exchange!(sys, 1.0,  Bond(1, 3, (0, 0, 0)), :Jc => 0)
set_onsite_coupling!(sys, S -> -S[1]^2,  1, :Kxx => 0)
set_onsite_coupling!(sys, S -> -S[2]^2,  1, :Kyy => 0)

# Parameters for Phase III as reported in Table I of Skoulatos et al. Here, it
# seems the paper has a typo. To get agreement in the calculated spin wave
# spectrum, it is necessary to swap the reported exchange parameters, ``J_{ab} =
# 4.24`` and ``J_{c} = 5.95`` (meV).

labels = [:Jab, :Jc, :Kxx, :Kyy]
skoulatos_fit = [5.95, 4.24, -0.48, -0.06]
set_params!(sys, labels, skoulatos_fit)

# Energy minimization yields the expected (π, π, π) Néel order with easy axis
# ẑ.

minimize_energy!(sys)
plot_spins(sys; ghost_radius=4)

# Magnon bands reproduce Fig. 1d of Skoulatos et al.

Q1 = [0, 1.0, 2.0]
Q2 = [0, 1.0, 3.0]
Q3 = [0, 1.5, 3.5]
path = q_space_path(cryst, [Q1, Q2, Q3], 500)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities_bands(swt, path)
plot_intensities(res)

# This tutorial will refit the model parameters using points on the dispersion
# curves. The energies below were extracted from Fig. 1d of Skoulatos et al.
# using the [WebPlotDigitizer](https://automeris.io/) tool.

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
];

# Use [`make_loss_fn`](@ref) to define an optimization target. The system
# already has the correct Néel magnetic order; to keep it fixed, the loss
# function should _not_ call [`minimize_energy!`](@ref). Sunny provides
# `labeled_peaks_mismatch` for robust comparison of excitation energies. The
# parameter `σ` is uncertainty of the experimental energies.

loss = make_loss_fn(sys, labels) do sys
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    res = intensities_bands(swt, qs)
    Sunny.labeled_peaks_mismatch(res, Es; σ=0.2)
end;

# Select some relatively non-informative parameter guess. Measure its loss
# (fitting mismatch) as a baseline.

guess = [3, 3, -0.2, -0.1] # Guess for [Jab, Jc, Kxx, Kyy]
loss(guess)

# The [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package provides a
# variety of powerful optimization methods. The particle swarm method is
# gradient free, has some ability to overcome local minima, and supports
# specification of lower and upper bounds on the optimization variables. Use 10
# independent runs of particle swarm to search for a good model.

import Optim

lower = [0., 0, -4, -4]
upper = [10., 10, 0, 0]
method = Optim.ParticleSwarm(; lower, upper, n_particles=10)
options = Optim.Options(; iterations=100)

fits = map(1:10) do i
    fit = Optim.optimize(loss, guess, method, options)
    println("Iteration $i, loss $(fit.minimum), params $(fit.minimizer)")
    fit
end;

# Fine tune the best model to about 4 digits of accuracy in meV.

options = Optim.Options(; g_tol=1e-4/units.meV)
best_fit = argmin(fit -> fit.minimum, fits)
best_fit = Optim.optimize(loss, best_fit.minimizer, Optim.NelderMead(), options)
best_fit.minimizer # [Jab, Jc, Kxx, Kyy]

# Compare the fitted spectrum to the experimentally measured peaks.

set_params!(sys, labels, best_fit.minimizer)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities_bands(swt, path)
fig = plot_intensities(res, ylims=(0, 40))
qinds = Sunny.find_qs_along_path(qs, path)
data_pairs = [(q, Eq) for (q, Eqs) in zip(qinds, Es) for Eq in Eqs]
plot!(fig[1, 1], data_pairs; color=:transparent, strokecolor=:magenta, strokewidth=3)
fig
