# Reproduce Fig. 3b,c of Do et al., Nat. Commun. 12, 5331 (2021)
# (https://arxiv.org/abs/2012.05445). Calculates the dynamical structure factor
# of the easy-plane s = 1 antiferromagnet Ba₂FeSi₂O₇ along [H, 0, 1/2],
# comparing GLSWT with its ``1/s`` correction. Each uses its own fitted
# parameter set.

using Sunny, LinearAlgebra
using GLMakie
load_blas_for_threading()

# Effective s = 1 model of the Fe sublattice, Table 1 of Do et al. The
# body-centered chemical cell is reshaped into the two-site primitive cell of
# the (π, π, π) Néel order.

cryst = Crystal(lattice_vectors(8.3193, 8.3193, 5.3348, 90, 90, 90), [[0, 0, 0]], 65; types=["Fe"])

function bfso_swt(J₁, D)
    sys = System(cryst, [1 => Moment(; s=1, g=diagm([2.18, 2.18, 1.93]))], :SUN)
    set_exchange!(sys, diagm([J₁, J₁, J₁/3]), Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, 0.1 * diagm([J₁, J₁, J₁/3]), Bond(1, 1, [0, 0, 1]))
    set_onsite_coupling!(sys, S -> D * S[3]^2, 1)
    sys = reshape_supercell(sys, [1 0 1/2; 0 1 1/2; 0 0 1])
    randomize_spins!(sys)
    minimize_energy!(sys; maxiters=2000)

    # The model is U(1) symmetric. Rotate the ordered moments onto ``±x``.
    θ = atan(sys.dipoles[1][2], sys.dipoles[1][1])
    U = exp(im * θ * Matrix(spin_matrices(1)[3]))
    for site in eachsite(sys)
        set_coherent!(sys, U * sys.coherents[site], site)
    end
    measure = ssf_perp(sys; apply_g=false, formfactors=[1 => FormFactor("Fe2")])
    return SpinWaveTheory(sys; measure)
end

swtA = bfso_swt(0.245, 1.61) # Parameter set 𝒜, fitted with GLSWT
swtB = bfso_swt(0.266, 1.42) # Parameter set ℬ, fitted with 1/s corrections

# Average over the two in-plane domains, ordered along ``±x`` and ``±y``.

rotations = [([0, 0, 1], 0), ([0, 0, 1], π/2)]
weights = [1, 1]
path = q_space_path(cryst, [[H, 0, 1/2] for H in 0:3], 180; labels=string.(0:3))

η = 0.1
energies = 0:0.02:3.5
res1 = domain_average(cryst, path; rotations, weights) do path
    intensities(swtA, path; energies, kernel=lorentzian(fwhm=2η))
end
res2 = domain_average(cryst, path; rotations, weights) do path
    Sunny.corrected_intensities(swtB, path; energies, η, tol=0.003, threaded=true, verbose=true)
end

# Compare with Fig. 3b,c of Do et al.

colormap = cgrad([:white, RGBf(0.3, 0.4, 1), :darkblue, :black])
fig = Figure(size=(700, 800))
axis = (; xlabel="H on (H, 0, 1/2) (r.l.u.)", ylabel="Energy (meV)")
plot_intensities!(fig[1, 1], res1; colormap, colorrange=(0, 3), title="GLSWT (set 𝒜)", axis)
plot_intensities!(fig[2, 1], res2; colormap, colorrange=(0, 3), title="GLSWT + one-loop (set ℬ)", axis)
fig
