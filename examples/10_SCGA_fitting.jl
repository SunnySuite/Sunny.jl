# # 10. Fitting to diffuse scattering data
#
# The self-consistent Gaussian approximation, [`SCGA`](@ref), calculates the
# static structure factor ``\mathcal{S}(𝐪)`` in the paramagnetic phase. These
# simulations can be directly compared to diffuse scattering data and offer a
# robust pathway to fitting model parameters.
# 
# Working in the paramagnetic phase has two strong advantages. First, one can
# avoid complications with finding the true magnetic ground state. Second, the
# associated ``\mathcal{S}(𝐪)`` data is inherently smooth, which helps to find
# the _globally_ optimal model parameters.
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

using Sunny, LinearAlgebra
import Zygote
import DifferentiationInterface as DI

latvecs = lattice_vectors(8.3342, 8.3342, 8.3342, 90, 90, 90)
positions = [[1/2, 1/2, 1/2]]
cryst = Crystal(latvecs, positions, 227)

sys = System(cryst, [1 => Moment(; s=3/2, g=2)], :dipole)
sys = reshape_supercell(sys, primitive_cell(cryst))
set_spin_rescaling_for_static_sum_rule!(sys)

# Exchanges from Bai's PRL
J1 = 3.27
set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]), :J1 => J1)
set_exchange!(sys, 1.0, Bond(1, 7, [0, 0, 0]), :J2 => 0.0815*J1)
set_exchange!(sys, 1.0, Bond(1, 3, [1, 0, 0]), :J3a => 0.1050*J1)
set_exchange!(sys, 1.0, Bond(1, 3, [0, 0, 0]), :J3b => 0.0085*J1)
labels = [:J1, :J2, :J3a, :J3b]
refvals = get_params(sys, labels)

formfactors = [1 => FormFactor("Cr3")]
measure = ssf_trace(sys; apply_g=false, formfactors)
kT = 20*meV_per_K
dq = 1/6

set_params!(sys, labels, refvals)
scga = SCGA(sys; measure, kT, dq)
grids = map([0, 1]) do l
    q_space_grid(cryst, [1, 0, 0], range(0, 4, 20), [0, 1, 0], range(0, 4, 20); offset=[0, 0, l])
end
refs = intensities_static.(Ref(scga), grids)

# using GLMakie
# fig = Figure(; size=(300, 800))
# plot_intensities!(fig[1, 1], refs[1])
# plot_intensities!(fig[2, 1], refs[2])

# println(round.(intensities_static(scga, grids[1]).data, digits=5))

refdata1 = 1.3 * refs[1].data
refdata2 = 1.3 * refs[2].data
refdata1[1] = NaN

loss = Sunny.fitting_loss(sys, labels) do sys
    scga = SCGA(sys; measure, kT, dq)
    res1 = intensities_static(scga, grids[1])
    res2 = intensities_static(scga, grids[2])
    return Sunny.squared_error((res1.data, res2.data), (refdata1, refdata2); rescale=true)
end

import Optim

# Optim options / logging
opts = Optim.Options(
    iterations = 500,
    g_tol      = 1e-8,
    show_trace = true,
    show_every = 5,
)

# Run optimization with reverse-mode AD
guess = [0.0, 0.0, 0.0, 0.0]
res = Optim.optimize(loss, guess, Optim.LBFGS(), opts; autodiff=DI.AutoZygote())

refvals - res.minimizer




guess = [0.0, 0.0, 0.0, 0.0]
loss(guess)

using Chairmarks
@b loss(guess)
@b DI.gradient(loss,  DI.AutoZygote(), guess)

@profview DI.gradient(loss,  DI.AutoZygote(), guess)