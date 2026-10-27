using Sunny, GLMakie

units = Units(:meV, :angstrom)
a = 5.2821 / sqrt(2)
b = 5.6144 / sqrt(2)
c = 7.5283
latvecs = lattice_vectors(a, b, c, 90, 90, 90)
positions = [[0, 0, 0], [0, 0, 1/2]]
types = ["V", "V"]
cryst = Crystal(latvecs, positions, 1; types)

moments = [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)]
sys = System(cryst, moments, :dipole_uncorrected; dims=(2,2,1))
Jab = 2.6
Jc  = 3.1
δ   = 0.35
K1  = 0.90
K2  = 0.97
d   = 1.15
Jc1 = [-Jc*(1+δ)+K2 0 -d; 0 -Jc*(1+δ) 0; +d 0 -Jc*(1+δ)]
Jc2 = [-Jc*(1-δ)+K2 0 +d; 0 -Jc*(1-δ) 0; -d 0 -Jc*(1-δ)]
set_exchange!(sys, Jab, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, Jab, Bond(2, 2, [1, 0, 0]))
set_exchange!(sys, Jc1, Bond(1, 2, [0, 0, 0]))
set_exchange!(sys, Jc2, Bond(2, 1, [0, 0, 1]))
set_exchange!(sys, Jab, Bond(1, 1, [0, 1, 0]))
set_exchange!(sys, Jab, Bond(2, 2, [0, 1, 0]))
set_onsite_coupling!(sys, S -> -K1*S[1]^2, 1)
set_onsite_coupling!(sys, S -> -K1*S[1]^2, 2)

view_crystal(sys)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[0.75, 0.75, 0], [0.5, 0.5, 0], [0.5, 0.5, 1]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units)
