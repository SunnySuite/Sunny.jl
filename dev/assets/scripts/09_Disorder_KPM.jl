using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.6"

latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
cryst = Crystal(latvecs, [[0, 0, 0]])
sys = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1))
set_exchange!(sys, +1.0, Bond(1, 1, [1,0,0]))

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles], ndims=2)

qs = [[0, 0, 0], [1/3, 1/3, 0], [1/2, 0, 0], [0, 0, 0]]
labels = ["Γ", "K", "M", "Γ"]
path = q_space_path(cryst, qs, 150; labels)

kernel = lorentzian(fwhm=0.4)
energies = range(0.0, 3.0, 150)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

sys_inhom = to_inhomogeneous(repeat_periodically(sys, (10, 10, 1)))

for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1,1,[1,0,0]))
    noise = randn()/3
    set_exchange_at!(sys_inhom, 1.0 + noise, site1, site2; offset)
end

minimize_energy!(sys_inhom, maxiters=5_000)
plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], ndims=2)

swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.01)
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

set_field!(sys_inhom, [0, 0, 7.5])
randomize_spins!(sys_inhom)
minimize_energy!(sys_inhom)

energies = range(0.0, 9.0, 150)
swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.1)
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

for site in eachsite(sys_inhom)
    noise = randn()/6
    sys_inhom.gs[site] = [1 0 0; 0 1 0; 0 0 1+noise]
end

swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.1)
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

set_field!(sys, [0, 0, 7.5])
randomize_spins!(sys)
minimize_energy!(sys)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, path; energies, kernel)
plot_intensities(res)
