using Sunny, GLMakie

units = Units(:meV, :angstrom)
a = b = 8.539 # (Å)
c = 5.2414
latvecs = lattice_vectors(a, b, c, 90, 90, 120)
types = ["Fe", "Nb", "Ba", "Si", "O", "O", "O"]
positions = [[0.24964,0,0.5], [0,0,0], [0.56598,0,0], [2/3,1/3,0.5220],
             [2/3,1/3,0.2162], [0.5259,0.7024,0.3536], [0.7840,0.9002,0.7760]]
langasite = Crystal(latvecs, positions, 150; types)
cryst = subcrystal(langasite, "Fe")
view_crystal(cryst)

sys = System(cryst, [1 => Moment(s=5/2, g=2)], :dipole; seed=0)
J₁ = 0.85
J₂ = 0.24
J₃ = 0.053
J₄ = 0.017
J₅ = 0.24
set_exchange!(sys, J₁, Bond(3, 2, [1,1,0]))
set_exchange!(sys, J₄, Bond(1, 1, [0,0,1]))
set_exchange!(sys, J₂, Bond(1, 3, [0,0,0]))

ϵT = -1
if ϵT == -1
    set_exchange!(sys, J₃, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₅, Bond(3, 2, [1,1,1]))
elseif ϵT == 1
    set_exchange!(sys, J₅, Bond(2, 3, [-1,-1,1]))
    set_exchange!(sys, J₃, Bond(3, 2, [1,1,1]))
else
    error("Chirality must be ±1")
end

axis = [0, 0, 1]
randomize_spins!(sys)
k = minimize_spiral_energy!(sys, axis)

sys_enlarged = repeat_periodically_as_spiral(sys, (1, 1, 7); k, axis)
plot_spins(sys_enlarged; color=[S[1] for S in sys_enlarged.dipoles])

measure = ssf_perp(sys)
swt = SpinWaveTheorySpiral(sys; measure, k, axis)

qs = [[0, 1, -1], [0, 1, -1+1], [0, 1, -1+2], [0, 1, -1+3]]
path = q_space_path(cryst, qs, 400)
energies = range(0, 6, 400)
res = intensities(swt, path; energies, kernel=gaussian(fwhm=0.25))
plot_intensities(res; units, saturation=0.7, colormap=:jet, title="Scattering intensities")

measure = ssf_custom_bm(sys; u=[0, 1, 0], v=[0, 0, 1]) do q, ssf
    imag(ssf[2,3] - ssf[3,2])
end
swt = SpinWaveTheorySpiral(sys; measure, k, axis)
res = intensities(swt, path; energies, kernel=gaussian(fwhm=0.25))
plot_intensities(res; units, saturation=0.8, allpositive=false,
                 title="Im[S²³(q, ω) - S³²(q, ω)]")
