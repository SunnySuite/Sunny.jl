function optim_energy(free_dipoles, sys::System{0})
    norm_penalty(y) = 1/y^2 - 2/y
    directions = reinterpret(reshape, SVector{3, Float64}, free_dipoles)
    E = 0.0
    for site in all_sites(sys)
        polarize_spin!(sys, directions[site], site)
        E += norm_penalty(norm(directions[site]))
    end
    return E + energy(sys)
end

function optim_gradient!(buf, free_dipoles, sys::System{0}) 
    function grad_norm_penalty(dipole)
        xx = dipole ⋅ dipole 
        2*((√xx-1)/xx^2)*dipole
    end
    free_dipoles = reinterpret(reshape, SVector{3, Float64}, free_dipoles)
    Hgrad = reinterpret(reshape, SVector{3, Float64}, buf)

    # Calculate gradient of energy in original coordinates
    for site in all_sites(sys)
        polarize_spin!(sys, free_dipoles[site], site)
    end
    Sunny.set_forces!(Hgrad, sys.dipoles, sys)

    # Calculate gradient in "new coordinates" and incorporate regularizing terms
    for site in all_sites(sys)
        xx = free_dipoles[site] ⋅ free_dipoles[site]
        jac = √xx * (I - (1/xx) * (free_dipoles[site] * free_dipoles[site]'))
        Hgrad[site] = -jac * Hgrad[site] + grad_norm_penalty(free_dipoles[site]) # Note Optim expects ∇, `set_forces!` gives -∇
    end
end

function minimize_energy!(sys; method=Optim.LBFGS)
    f(spins) = optim_energy(spins, sys)
    g!(G, spins) = optim_gradient!(G, spins, sys)
    free_dipoles = Array(reinterpret(reshape, Float64, sys.dipoles))
    Optim.optimize(f, g!, free_dipoles, method()) 
end