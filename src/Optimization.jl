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
# Stereographic stuff
################################################################################
struct StereographicPoint{N}
    data :: NTuple{N, Float64}
end

function StereographicPoint(i, coords::NTuple{Nm1, Float64}) where Nm1
    StereographicPoint{Nm1+1}((coords..., Float64(i)))
end

Base.getindex(sp::StereographicPoint, i) = getindex(sp.data, i)

function Base.show(io::IO, ::MIME"text/plain", sp::StereographicPoint{Nm1}) where Nm1
    println(io, "($(sp.data[1]), $(sp.data[2])) in ⟂$(Int(sp.data[3]))-plane")
end

@inline remaining_indices(skip, max) = ntuple(i -> i > (skip-1) ? i + 1 : i, max-1)
@inline chart(sp::StereographicPoint) = Int64(sp.data[end])

function vec_to_stereo(vec::SVector{N, Float64}) where N  # N parameter unnecessary here, but generalizes to ℝP^{N-1}
    plane_idx = argmin(vec) # Project onto closest plane, i.e. the plane 
                            # perpendicular to the axis corresponding
                            # to the minimal vector component (⟂ ̂x_min)
    other_idxs = remaining_indices(plane_idx, N)   
    denom = 1 - vec[plane_idx]             
    StereographicPoint{N}(ntuple(i -> i < N ? vec[other_idxs[i]]/denom : Float64(plane_idx), N)) 
end

function stereo_to_vec(sp::StereographicPoint{N}) where N
    c = chart(sp)
    r2 = mapreduce(x -> x^2, +, sp.data[1:end-1])
    SVector{N,Float64}(
        ntuple(N) do i
            if i < c
                2sp[i]/(1+r2)
            elseif i == c
                (-1+r2)/(1+r2)
            else
                2sp[i-1]/(1+r2)
            end
        end
    )
end

@inline remap_stereo(sp) = sp |> stereo_to_vec |> vec_to_stereo

function stereo_jac!(jac, sp::StereographicPoint{N}) where N # Hardy-har
    δ(i,j) = i == j ? 1 : 0
    Nidx(i, c) = i > c ? i - 1 : i

    c = chart(sp)
    r2 = mapreduce(x -> x^2, +, sp.data[1:end-1])
    for j in 1:N, k in 1:N-1
        jac[k,j] = if j == c
            (2sp[k]/(r2+1))*(1 - (r2-1)/(r2+1))
        else
            (2/(r2+1))*(δ(Nidx(j,c),k) - (2sp[Nidx(j,c)]*sp[k])/(r2+1))
        end
    end
end

################################################################################
# Optimization with stereographic coordinates
################################################################################
function optim_energy_stereo(stereo_coords, sys::System{0})
    stereo_coords = reinterpret(reshape, StereographicPoint{3}, stereo_coords)
    for site in all_sites(sys)
        polarize_spin!(sys, stereo_to_vec(stereo_coords[site]), site)
    end
    return energy(sys)
end

function optim_gradient_stereo!(buf, stereo_coords, sys::System{0}) 
    stereo_coords = reinterpret(reshape, StereographicPoint{3}, stereo_coords)
    Hgrad = reinterpret(reshape, SVector{3, Float64}, buf)

    # Calculate gradient of energy in original coordinates
    for site in all_sites(sys)
        polarize_spin!(sys, stereo_to_vec(stereo_coords[site]), site)
    end
    Sunny.set_forces!(Hgrad, sys.dipoles, sys)

    # Calculate gradient, using Jacobian for stereographic coordinates 
    jac = zeros(3,3)  # Maybe write function to apply matrix multiply and avoid allocation
    for site in all_sites(sys)
        stereo_jac!(jac, stereo_coords[site])
        Hgrad[site] = -jac * Hgrad[site]  # Note Optim expects ∇, `set_forces!` gives -∇
    end
end



# Under construction
function minimize_energy_stereo!(sys; method=Optim.LBFGS, maxiters = 100, kwargs...)
    f(stereo_spins) = optim_energy_stereo(stereo_spins, sys)
    g!(G, stereo_spins) = optim_gradient_stereo!(G, stereo_spins, sys)
    # fg!(energy_spins, G, x) = stereo_energy_grad!(stereo_spins, G, sys)

    options = Optim.Options(iterations=10, kwargs...)

    # Quick test if any coordinates every go to infinity 
    niters = 0
    success = false
    while niters < maxiters
        niters += 1
        stereo_spins = map(vec -> vec_to_stereo(vec), sys.dipoles) 
        stereo_spins = Array(reinterpret(reshape, Float64, stereo_spins))
        # stereo_spins = map!(vec -> vec_to_stereo(vec), stereo_spins, sys.dipoles) 
        optout = Optim.optimize(f, g!, stereo_spins, method(), options)
        if optout.g_converged
            println("We got there.")
            return true
        end
        # println("Maximum: ", maximum(stereo_spins))
        # println(stereo_spins)
        if maximum(stereo_spins) > 5.0
            println("Remapping coords...")
            stereo_spins = reinterpret(reshape, StereographicPoint{3}, stereo_spins)
            @. stereo_spins = remap_stereo(stereo_spins)
        end
        println("Trying again...")
    end

    return success
end