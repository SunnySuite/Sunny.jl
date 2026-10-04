# Reproduce Fig. 4 of Veillette, James and Essler, PRB 72, 134429 (2005)
# (https://arxiv.org/abs/cond-mat/0506667). Calculates the magnon dispersion of
# Cs₂CuCl₄ at order 1/s.
#
# The paper defines the renormalized dispersion as the pole of the full Dyson
# equation, their Eq. (40), which is what `dyson=:nambu` solves. It is read off
# here as the peak of the spectral function, not from the first-order poles of
# `corrected_intensities_bands`, which differ at s = 1/2 by up to 0.02 meV.

using Sunny, LinearAlgebra
using GLMakie
load_blas_for_threading()

# One layer of the anisotropic triangular lattice: J along b (here x), and J′
# on the zig-zag bonds, which carry a DM vector along a (here z). The cycloid is
# made commensurate at k = 5/9 in a 9-cell supercell. This requires J′ = 0.1334
# meV rather than 0.128 meV, which makes k the classical minimum.

(J, D, k) = (0.374, 0.020, 5/9)
Jp = (-J*sin(2π*k) - D*cos(π*k)) / sin(π*k)
cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 90), [[0, 0, 0], [1/2, 1/2, 0]], 1)
sys = System(cryst, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :dipole)
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J, Bond(2, 2, [1, 0, 0]))
for b in (Bond(1, 2, [0, 0, 0]), Bond(1, 2, [0, -1, 0]), Bond(2, 1, [1, 1, 0]), Bond(2, 1, [1, 0, 0]))
    set_exchange!(sys, Jp*I(3) + dmvec([0, 0, -D]), b)
end
sys = repeat_periodically(sys, (9, 1, 1))
for site in eachsite(sys)
    θ = 2π * k * global_position(sys, site)[1]
    set_dipole!(sys, [cos(θ), sin(θ), 0], site)
end

# The supercell folds the dispersion into 18 bands. The principal mode, which
# Fig. 4 shows, is polarized along a, and so is the peak of S^aa.

measure = ssf_custom((q, ssf) -> real(ssf[3, 3]), sys; apply_g=false)
swt = SpinWaveTheory(sys; measure)

# The path of Fig. 4 in the paramagnetic Brillouin zone

path = q_space_path(cryst, [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 0, 0]], 150;
                    labels=["(000)", "(010)", "(011)", "(000)"])
η = 0.01
energies = 0:(η/4):0.55

# The brightest band of LSWT is the principal mode. Each corrected peak is
# sought within 0.1 meV of it, which keeps a folded band from taking over near
# a Goldstone mode, where the principal mode carries little weight.
ω0 = let b = intensities_bands(swt, path)
    [b.disp[argmax(b.data[:, iq]), iq] for iq in axes(b.data, 2)]
end
function peaks(res)
    return map(axes(res.data, 2)) do iq
        win = findall(e -> abs(e - ω0[iq]) < 0.1, energies)
        energies[win[argmax(res.data[win, iq])]]
    end
end

# A `:nambu` pole near ω = 0 at (010) is flagged as a breakdown of the
# resummation. It belongs to a folded band far from the principal mode, so the
# flag is disabled here.

corrected = peaks(Sunny.corrected_intensities(swt, path; energies, η, mark_breakdown=false, threaded=true))

fig = Figure()
ax = Axis(fig[1, 1]; ylabel="Energy (meV)", xticks=path.xticks)
lines!(ax, ω0; color=:blue, linestyle=:dash, label="LSWT")
lines!(ax, corrected; color=:red, label="1/s")
ylims!(ax, 0, 0.5)
axislegend(ax; position=(0.4, 1))
fig

# Against the solid line of Fig. 4, digitized, the median deviation is 0.005
# meV, about 1% of the bandwidth. Deviations up to 0.03 meV occur only on the
# steep slopes next to (000), where a small offset in wavevector is a large one
# in energy. Along (010)-(011) the result sits about 0.01 meV high, of which
# 0.004 meV is the commensurate approximation (k = 5/9 rather than 0.553). The
# paper also renormalizes the magnon energies inside the loops, iterating to
# self-consistency, to move the threshold of the two-magnon continuum. That is
# omitted here: it changes the dispersion by less than 0.01 meV.
