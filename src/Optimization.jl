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
        ixx = 1/(free_dipoles[site] ⋅ free_dipoles[site])
        jac = √ixx * (I - ixx * (free_dipoles[site] * free_dipoles[site]'))
        Hgrad[site] = -jac * Hgrad[site] + grad_norm_penalty(free_dipoles[site]) # Note Optim expects ∇, `set_forces!` gives -∇
    end
end

"""
    minimize_energy!(sys; method=Optim.LBFGS, kwargs...)

Minimize the energy of a spin system using either LBFGS (`method=Optim.LBFGS`)
or Conjugate Gradient (`method=Optim.ConjugateGradient`) methods. Currently only
works for systems in dipole mode. 
"""
function minimize_energy!(sys; method=Optim.LBFGS, kwargs...)
    f(spins) = optim_energy(spins, sys)
    g!(G, spins) = optim_gradient!(G, spins, sys)
    free_dipoles = Array(reinterpret(reshape, Float64, sys.dipoles))
    Optim.optimize(f, g!, free_dipoles, method(), kwargs...) 
end

################################################################################
# Coordinate systems
################################################################################

struct StereographicPoint{Nm1}
    chart :: Int
    coord :: SVector{Nm1, Float64}
end

function vec_to_stereo(vec::SVector{N, Float64}) where N  # N parameter unnecessary here, but generalizes to ℝP^{N-1}
    remaining_indices(skip, max) = ntuple(i -> i > (skip-1) ? i + 1 : i, max-1)

    plane_idx = argmax(vec)               
    is = remaining_indices(plane_idx, N)   
    denom = 1 - vec[plane_idx]             
    StereographicPoint{N-1}(plane_idx, ntuple(i -> vec[is[i]]/denom, N-1)) 
end

function stereo_to_vec(sp::StereographicPoint{Nm1}) where Nm1
    r2 = mapreduce(x -> x^2, +, sp.coord)
    SVector{Nm1+1,Float64}(
        ntuple(Nm1 + 1) do i
            if i < sp.chart
                2sp.coord[i]/(1+r2)
            elseif i == sp.chart
                (-1 + r2)/(1+r2)
            else
                2sp.coord[i-1]/(1+r2)
            end
        end
    )
end



# Tests
# for _ in 1:20
#     N = 3
#     vec = SVector{N, Float64}(normalize(rand(N)))
#     sp = vec_to_stereo(vec)
#     println(stereo_to_vec(sp) ≈ vec ? "Good" : "Bad")
# end
# 
# begin
#     N = 3
#     vec = SVector{N, Float64}(normalize(rand(N)))
#     sp = vec_to_stereo(vec)
#     @btime vec_to_stereo($vec)
#     vec1 = stereo_to_vec(sp)
#     @btime stereo_to_vec($sp)
# end