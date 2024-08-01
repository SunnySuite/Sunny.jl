using Sunny, GLMakie

units = Units(:meV)
a = b = 8.539 # (Å)
c = 5.2414
latvecs = lattice_vectors(a, b, c, 90, 90, 120)
types = ["Fe", "Nb", "Ba", "Si", "O", "O", "O"]
positions = [[0.24964,0,0.5], [0,0,0], [0.56598,0,0], [2/3,1/3,0.5220],
             [2/3,1/3,0.2162], [0.5259,0.7024,0.3536], [0.7840,0.9002,0.7760]]
langasite = Crystal(latvecs, positions, 150; types)
cryst = subcrystal(langasite, "Fe")
view_crystal(cryst)

latsize = (1,1,7)
S = 5/2
seed = 5
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole)

J₁ = 0.85
J₂ = 0.24
J₃ = 0.053
J₄ = 0.017
J₅ = 0.24
set_exchange!(sys, J₁, Bond(3, 2, [1,1,0]))
set_exchange!(sys, J₄, Bond(1, 1, [0,0,1]))
set_exchange!(sys, J₂, Bond(1, 3, [0,0,0]))

ϵD = -1
ϵH = +1
ϵT = ϵD * ϵH

if ϵT == -1
    set_exchange!(sys, J₃, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₅, Bond(3, 2, [1,1,1]))
elseif ϵT == 1
    set_exchange!(sys, J₅, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₃, Bond(3, 2, [1,1,1]))
else
    throw("Provide a valid chirality")
end

k = [0, 0, 1/7]
axis = [0,0,1]
set_spiral_order_on_sublattice!(sys, 1; k, axis, S0=[1, 0, 0])
set_spiral_order_on_sublattice!(sys, 2; k, axis, S0=[-1/2, -sqrt(3)/2, 0])
set_spiral_order_on_sublattice!(sys, 3; k, axis, S0=[-1/2, +sqrt(3)/2, 0])

plot_spins(sys; color=[s[1] for s in sys.dipoles])

qs = [[0, 1, -1], [0, 1, -1+1], [0, 1, -1+2], [0, 1, -1+3]]
path = q_space_path(cryst, qs, 600)

measure = ssf_perp(sys)
swt = SpinWaveTheory(sys; measure)
energies = range(0, 6, 400)  # 0 < ω < 6 (meV)
res = intensities(swt, path; energies, kernel=gaussian2(fwhm=0.25))
axisopts = (; title=L"$ϵ_T=-1$, $ϵ_Δ=-1$, $ϵ_H=+1$", titlesize=20)
plot_intensities(res; units, axisopts, saturation=0.7, colormap=:jet)

const q_i = cryst.recipvecs * [0, 0, 1]
measure = ssf_custom(sys) do q, ssf
    (_, y, z) = blume_maleev(q_i, q)
    imag(y'*ssf*z - z'*ssf*y)
end

swt = SpinWaveTheory(sys; measure)
res = intensities(swt, path; energies, kernel=gaussian2(fwhm=0.25))
axisopts = (; title=L"$ϵ_T=-1$, $ϵ_Δ=-1$, $ϵ_H=+1$", titlesize=20)
plot_intensities(res; units, axisopts, saturation=0.7, colormap=:bwr)
