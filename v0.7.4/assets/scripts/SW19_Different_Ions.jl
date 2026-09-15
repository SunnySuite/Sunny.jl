using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.4"

units = Units(:meV, :angstrom)
a = 3.0
b = 8.0
c = 4.0
latvecs = lattice_vectors(a, b, c, 90, 90, 90)
positions = [[0, 0, 0], [0, 1/2, 0]]
types = ["Cu2", "Fe2"]
cryst = Crystal(latvecs, positions, 1; types)
view_crystal(cryst)

J_Cu_Cu = 1.0
J_Fe_Fe = 1.0
J_Cu_Fe = -0.1
moments = [1 => Moment(s=1/2, g=2), 2 => Moment(s=2, g=2)]
sys = System(cryst, moments, :dipole; dims=(2, 1, 1))
set_exchange!(sys, J_Cu_Cu, Bond(1, 1, [-1, 0, 0]))
set_exchange!(sys, J_Fe_Fe, Bond(2, 2, [-1, 0, 0]))
set_exchange!(sys, J_Cu_Fe, Bond(2, 1, [0, 1, 0]))
set_exchange!(sys, J_Cu_Fe, Bond(1, 2, [0, 0, 0]))

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)

qs = [[0,0,0], [1,0,0]]
path = q_space_path(cryst, qs, 400)

fig = Figure(size=(768,600))

formfactors = [1 => FormFactor("Cu2"), 2 => FormFactor("Fe2")]
swt = SpinWaveTheory(sys; measure=ssf_trace(sys; formfactors))
res = intensities_bands(swt, path)
plot_intensities!(fig[1, 1], res; units, title="All correlations")

formfactors = [1 => FormFactor("Cu2"), 2 => zero(FormFactor)]
swt = SpinWaveTheory(sys; measure=ssf_trace(sys; formfactors))
res = intensities_bands(swt, path)
plot_intensities!(fig[1, 2], res; units, title="Cu-Cu correlations")

formfactors = [1 => zero(FormFactor), 2 => FormFactor("Fe2")]
swt = SpinWaveTheory(sys; measure=ssf_trace(sys; formfactors))
res = intensities_bands(swt, path)
plot_intensities!(fig[2, 2], res; units, title="Fe-Fe correlations")

fig

Sunny.magnetization_lswt_correction_dipole(swt; atol=1e-4)
