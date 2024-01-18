using Sunny, GLMakie, LinearAlgebra

function test_system(; dims=(4, 4, 1), J=1.0, D=0.2, h=2.0, S=1/2, mode=:dipole)
    units = Sunny.Units.theory
    crystal = Crystal(diagm([1.0, 1.0, 1.1]), [[0,0,0]])
    sys = System(crystal, dims, [SpinInfo(1; S, g=2)], mode; units)
    set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
    S = mode == :dipole_large_S ? Inf : S
    set_onsite_coupling!(sys, D*spin_matrices(S)[3]^2, 1)
    set_external_field!(sys, [0, 0, h])
    return sys, crystal
end

begin
    sys, crystal = test_system(; J=1.0, S=1/2, D=-0.1, h=0.3, mode=:dipole)

    # Find ground state
    randomize_spins!(sys)
    minimize_energy!(sys; maxiters=10_000)
    plot_spins(sys)
end

# Examine SWT dispersion
begin
    swt = SpinWaveTheory(sys)
    γ = 0.05
    formula = intensity_formula(swt, :perp; kernel=lorentzian(γ))

    path, labels = reciprocal_space_path(crystal, [[0, 0, 0], [1/2, 1/2, 0], [1, 1, 0]], 100)
    energies = range(0.0, 4.0, 200)
    is = intensities_broadened(swt, path, energies, formula)
    heatmap(1:size(is, 1), energies, is; axis=(xticks=labels,), colorrange=(0.0, 10.0))
end

# Calcualte S(q,ω) with classical dynamics
nω = 200
ωmax = 4.0 
Δt = 0.05 

sys, crystal = test_system(; J=1.0, S=1/2, D=-0.1, h=0.3, mode=:dipole, dims=(16, 16, 1))
sc = dynamical_correlations(sys; nω, ωmax, Δt, process_trajectory=:symmetrize)
randomize_spins!(sys)
minimize_energy!(sys; maxiters=10_000)

λ = 0.1
kT = 0.1
integrator = Langevin(Δt/2; λ, kT)

# Thermalize
for _ in 1:5000
    step!(sys, integrator)
end

# Sample
for _ in 1:10
    for _ in 1:2000
        step!(sys, integrator)
    end
    add_sample!(sc, sys)
end

formula = intensity_formula(sc, :perp)
is_cl = intensities_interpolated(sc, path, formula)
heatmap(is_cl; colorrange=(0.0, 0.01))