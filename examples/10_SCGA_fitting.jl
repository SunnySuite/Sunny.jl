# # 10. Fitting to diffuse scattering data
#
# The self-consistent Gaussian approximation, [`SCGA`](@ref), calculates the
# static structure factor ``\mathcal{S}(𝐪)`` in the paramagnetic phase. These
# simulations can be directly compared to diffuse scattering data and offer a
# robust pathway to fitting model parameters.
# 
# Working in the paramagnetic phase is convenient because the smoothness of
# ``\mathcal{S}(𝐪)`` facilitates optimization. Also, unlike spin wave theory,
# it is not to solve for the magnetic ground state, which may change as a
# function of the model parameters.
#
# This tutorial uses SCGA to fit inelastic neutron scattering data for the
# frustrated pyrochlore antiferromagnet MgCr₂O₄. The fitted exchange
# interactions, up to third nearest neighbor, reproduce work in [Bai et al.,
# Phys. Rev. Lett. 122, 097201
# (2019)](https://doi.org/10.1103/PhysRevLett.122.097201).
#
# Sunny's SCGA calculator is inspired by the Spinteract code [J. Paddison,
# J.Phys.: Condens. Matter **35**, 495802
# (2023)](https://doi.org/10.1088/1361-648X/acf261).

using Sunny, GLMakie, LinearAlgebra

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(8.3342, 8.3342, 8.3342, 90, 90, 90)
positions = [[1/2, 1/2, 1/2]]
cryst = Crystal(latvecs, positions, 227)

sys = System(cryst, [1 => Moment(; s=3/2, g=2)], :dipole)
sys = reshape_supercell(sys, primitive_cell(cryst))
set_spin_rescaling_for_static_sum_rule!(sys)

# Exchanges from Bai's PRL
set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]), :J1 => 0)
set_exchange!(sys, 1.0, Bond(1, 7, [0, 0, 0]), :J2 => 0)
set_exchange!(sys, 1.0, Bond(1, 3, [1, 0, 0]), :J3a => 0)
set_exchange!(sys, 1.0, Bond(1, 3, [0, 0, 0]), :J3b => 0)
labels = [:J1, :J2, :J3a, :J3b]

formfactors = [1 => FormFactor("Cr3")]
measure = ssf_perp(sys; formfactors)
dq = 1/8

using Serialization: serialize, deserialize
(; refdata, centers) = open(joinpath(@__DIR__, "10_MgCr2O4_SQ.bin"), "r") do io
    deserialize(io)
end

grid = q_space_grid(cryst, [1, 0, 0], centers, [0, 1, 0], centers, [0, 0, 1], centers)
grid.qs[isnan.(refdata)] .*= NaN

loss = make_loss_fn(sys, labels) do sys
    scga = SCGA(sys; measure, kT=20*units.K, dq)
    res = intensities_static(scga, grid)
    diffuse_error = squared_error(res.data, refdata; rescale=true)

    χ = map([100, 200, 300] * units.K) do kT
        scga = SCGA(sys; measure, kT, dq)
        magnetic_susceptibility_per_site(scga)[1, 1] / units.cgs_molar_susceptibility
    end
    χ_ref = [0.00375, 0.00320, 0.00274]
    susceptibility_error = squared_error(χ, χ_ref)

    return 0.5 * diffuse_error + 0.5 * susceptibility_error
end

loss([0, 0, 0, 0])

import Optim
import Zygote
import DifferentiationInterface as DI

# Optim options / logging
opts = Optim.Options(
    iterations = 500,
    g_tol      = 1e-6,
    show_trace = true,
    show_every = 5,
)

# Run optimization with reverse-mode AD
guess = [0.0, 0.0, 0.0, 0.0]
# opt = Optim.optimize(loss, guess, Optim.LBFGS(), opts)
opt = Optim.optimize(loss, guess, Optim.LBFGS(), opts; autodiff=DI.AutoZygote())

# Comparison with previous fit
J1 = 3.27
previous_fit = [J1, 0.0815*J1, 0.1050*J1, 0.0085*J1]
@show loss(opt.minimizer)
@show loss(previous_fit)
@show opt.minimizer
@show previous_fit

set_params!(sys, labels, opt.minimizer)
grid = q_space_grid(cryst, [1, 0, 0], centers, [0, 1, 0], centers, [0, 0, 1], centers)
scga = SCGA(sys; measure, kT, dq)
res = intensities_static(scga, grid)

j1 = 1
j2 = 4
fig = Figure()
xlabel = "[H, 0, 0]"
ylabel = "[0, K, 0]"
title1 = "L=$(centers[j1])"
title2 = "L=$(centers[j2])"
heatmap(fig[1, 1], centers, centers, reverse(refdata[:, :, j1]; dims=1); colormap=:gnuplot2, axis=(; aspect=1, xlabel, ylabel, title=title1))
heatmap(fig[2, 1], centers, centers, reverse(refdata[:, :, j2]; dims=1); colormap=:gnuplot2, axis=(; aspect=1, xlabel, ylabel, title=title2))
heatmap(fig[1, 2], centers, centers, res.data[:, :, j1]; colormap=:gnuplot2, axis=(; aspect=1))
heatmap(fig[2, 2], centers, centers, res.data[:, :, j2]; colormap=:gnuplot2, axis=(; aspect=1))


using Chairmarks
import FiniteDiff
@b loss(guess)
@b DI.gradient(loss,  DI.AutoZygote(), guess)
@b DI.gradient(loss,  DI.AutoFiniteDiff(), guess)



# units = Units(:meV, :angstrom)
# J1 = 3.27
# previous_fit = [J1, 0.0815*J1, 0.1050*J1, 0.0085*J1]
# set_params!(sys, labels, previous_fit)
# measure = ssf_custom((q, ssf) -> ssf, sys)
# kT = 20*meV_per_K
# dq = 1/6
# scga = SCGA(sys; measure, kT, dq)
# res = intensities_static(scga, [[0, 0, 0]]).data[1] / (16 * kT)


# magnetic_susceptibility_per_site(scga) / units.cgs_molar_susceptibility



kTs = range(20, 400, length=100)

set_params!(sys, labels, opt.minimizer)
χ = map(kTs * units.K) do kT
    scga = SCGA(sys; measure, kT, dq)
    magnetic_susceptibility_per_site(scga)[1, 1] / units.cgs_molar_susceptibility
end

set_params!(sys, labels, previous_fit)
χ_xb = map(kTs * units.K) do kT
    scga = SCGA(sys; measure, kT, dq)
    magnetic_susceptibility_per_site(scga)[1, 1] / units.cgs_molar_susceptibility
end

lines(kTs, χ, label="New fit")
lines!(kTs, χ_xb, label="Bai fit")
lines!(kTs_ref, χ_ref, label="Experiment")
axislegend()
xlims!(0, 400)
ylims!(0, 6e-3)
