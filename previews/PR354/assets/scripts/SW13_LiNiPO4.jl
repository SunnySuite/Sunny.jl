using Sunny, GLMakie
@assert pkgversion(Sunny) >= v"0.7.6"

units = Units(:meV, :angstrom)
a = 10.02
b = 5.86
c = 4.68
latvecs = lattice_vectors(a, b, c, 90, 90, 90)
positions = [[1/4, 1/4, 0]]
types = ["Ni"]
cryst = Crystal(latvecs, positions, 62; types)
view_crystal(cryst)

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole_uncorrected)
Jbc =  1.036
Jb  =  0.6701
Jc  = -0.0469
Jac = -0.1121
Jab =  0.2977
Da  =  0.1969
Db  =  0.9097
set_exchange!(sys, Jbc, Bond(2, 3, [0, 0, 0]))
set_exchange!(sys, Jc, Bond(1, 1, [0, 0, -1]))
set_exchange!(sys, Jb, Bond(1, 1, [0, 1, 0]))
set_exchange!(sys, Jab, Bond(1, 2, [0, 0, 0]))
set_exchange!(sys, Jab, Bond(3, 4, [0, 0, 0]))
set_exchange!(sys, Jac, Bond(3, 1, [0, 0, 0]))
set_exchange!(sys, Jac, Bond(4, 2, [0, 0, 0]))
set_onsite_coupling!(sys, S -> Da*S[1]^2 + Db*S[2]^2, 1)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles])

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[0, 1, 0], [2, 1, 0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
fig = Figure(size=(768, 300))
plot_intensities!(fig[1, 1], res; units);

data_sorted = sort(res.data; dims=1, by= >(1e-12))
ax = Axis(fig[1, 2], xlabel="Momentum (r.l.u.)", ylabel="Intensity",
          xticks=res.qpts.xticks, xticklabelrotation=π/6)
lines!(ax, data_sorted[end, :]; label="Lower band")
lines!(ax, data_sorted[end-1, :]; label="Upper band")
axislegend(ax)
fig

qs =  [[0, 1, 0], [0, 1, 2]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
fig = Figure(size=(768, 300))
plot_intensities!(fig[1, 1], res; units)

data_sorted = sort(res.data; dims=1, by=x->abs(x)>1e-12)
ax = Axis(fig[1, 2], xlabel="Momentum (r.l.u.)", ylabel="Intensity",
          xticks=res.qpts.xticks, xticklabelrotation=π/6)
lines!(ax, data_sorted[end, :]; label="Lower band")
lines!(ax, data_sorted[end-1, :]; label="Upper band")
axislegend(ax)
fig
