# Reproduce the GNLSW panels of Fig. 3a,b of Bai et al., Nat. Commun. 14, 4199
# (2023) (https://arxiv.org/abs/2107.05694). Calculates the dynamical structure
# factor of FeI₂ in fields of 3 T and 4 T along c, including 1/s corrections. At
# 4 T, magnon decay broadens the tops of the E4 and E6 bands near H = -1/2 and
# -3/2.
#
# The reference calculation kept only cubic vertices. Sunny also includes the
# quartic (Hartree-Fock) shifts, which raise the bands by less than 0.1 meV.

using Sunny, LinearAlgebra
using GLMakie
load_blas_for_threading()

# The model of the FeI₂ tutorial. The paper quotes g = 3.8(5), but its band
# energies are reproduced with g = 4.0. At g = 3.8, the two lowest bands sit
# about 0.1 meV too high at 4 T.

units = Units(:meV, :angstrom)
latvecs = lattice_vectors(4.05012, 4.05012, 6.75214, 90, 90, 120)
cryst = Crystal(latvecs, [[0, 0, 0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]; types=["Fe", "I", "I"])
cryst = subcrystal(cryst, "Fe")

sys = System(cryst, [1 => Moment(s=1, g=4.0)], :SUN)
J1pm, J1pmpm, J1zpm, J1zz = -0.236, -0.161, -0.261, -0.236
set_exchange!(sys, [J1pm+J1pmpm 0 0; 0 J1pm-J1pmpm J1zpm; 0 J1zpm J1zz], Bond(1, 1, [1, 0, 0]))
for (Jpm, Jzz, d) in [(0.026, 0.113, [1, 2, 0]), (0.166, 0.211, [2, 0, 0]), (0.037, -0.036, [0, 0, 1]),
                      (0.013, 0.051, [1, 0, 1]), (0.068, 0.073, [1, 2, 1])]
    set_exchange!(sys, diagm([Jpm, Jpm, Jzz]), Bond(1, 1, d))
end
set_onsite_coupling!(sys, S -> -2.165 * S[3]^2, 1)
sys = reshape_supercell(sys, [1 0 0; 0 1 -2; 0 1 2])

# The path ``(H, 1/2 - H/2, 0)`` for ``-2 ≤ H ≤ 0``.

path = q_space_path(cryst, [[H, 1/2 - H/2, 0] for H in -2:0], 200; labels=string.(-2:0))

# For each field, minimize to the single-domain magnetic order and calculate
# intensities to order 1/s. The broadening η = 0.1 is a Lorentzian half-width,
# approximating the 0.2 meV instrumental resolution.

η = 0.1
energies = 0:0.025:8
measure = ssf_perp(sys; formfactors=[1 => FormFactor("Fe2")])
res = map([3, 4]) do B
    set_field!(sys, [0, 0, B * units.T])
    randomize_spins!(sys)
    minimize_energy!(sys)
    swt = SpinWaveTheory(sys; measure)
    Sunny.corrected_intensities(swt, path; energies, η, threaded=true, verbose=true)
end

# Compare with the right-most panels of Fig. 3a,b of Bai et al., each
# normalized to its maximum intensity.

parula = cgrad([RGBf(0.24, 0.15, 0.66), RGBf(0.28, 0.32, 0.96), RGBf(0.18, 0.53, 0.97), RGBf(0.07, 0.69, 0.84),
                RGBf(0.22, 0.78, 0.59), RGBf(0.67, 0.78, 0.22), RGBf(1.0, 0.77, 0.22), RGBf(0.98, 0.98, 0.08)])
fig = Figure(size=(700, 600))
axis = (; aspect=1/2, xlabel="H on (H, 1/2 - H/2, 0) (r.l.u.)", ylabel="Energy (meV)")
for (i, (B, r)) in enumerate(zip([3, 4], res))
    r.data ./= maximum(r.data)
    plot_intensities!(fig[1, i], r; colormap=parula, colorrange=(0, 1), title="μ₀H = $B T", axis)
end
fig
