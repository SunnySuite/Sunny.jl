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
# Coordinates for ℝP^{N-1} 
################################################################################
struct StereographicPoint{N}
    data :: NTuple{N, Float64}
end

function StereographicPoint(i, coords::NTuple{N, Float64}) where N
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
        set_coherent_state!(sys, stereo_to_vec(stereo_coords[site]), site)
    end
    return energy(sys)
end

function optim_gradient_stereo!(buf, stereo_coords, sys::System{0}) 
    stereo_coords = reinterpret(reshape, StereographicPoint{3}, stereo_coords)
    Hgrad = reinterpret(reshape, SVector{3, Float64}, buf)

    # Calculate gradient of energy in original coordinates
    for site in all_sites(sys)
        set_coherent_state!(sys, stereo_to_vec(stereo_coords[site]), site)
    end
    Sunny.set_forces!(Hgrad, sys.coherents, sys)

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

################################################################################
# Equivalent of set forces for SU(N) systems
################################################################################

# Calculates ℌz 
function set_complex_forces!(HZ, Z, sys::System{N}) where N
    B = only(get_dipole_buffers(sys, 1) )

    @. sys.dipoles = expected_spin(Z) # temporarily desyncs dipoles and coherents
    set_forces!(B, sys.dipoles, sys)

    if is_homogeneous(sys)
        ints = interactions_homog(sys)
        for site in all_sites(sys)
            Λ = ints[to_atom(site)].aniso.matrep
            HZ[site] = mul_spin_matrices(Λ, -B[site], -Z[site])  # Z negative to match sign convention of set_forces!
        end 
    else
        ints = interactions_inhomog(sys)
        for site in all_sites(sys)
            Λ = ints[site].aniso.matrep
            HZ[site] = mul_spin_matrices(Λ, -B[site], -Z[site])  # Z negative to match sign convention of set_forces!
        end 
    end
end

################################################################################
# Coordinates for ℂP^{N-1} 
################################################################################

struct CPAffineCoordinate{N}
    data :: NTuple{N, ComplexF64}
end

function CPAffineCoordinate(i, coords::NTuple{Nm1, ComplexF64}) where Nm1
    CPAffineCoordinate{Nm1+1}((coords..., ComplexF64(i)))
end

Base.getindex(sp::CPAffineCoordinate, i) = getindex(sp.data, i)

function Base.show(io::IO, ::MIME"text/plain", cp::CPAffineCoordinate{N}) where N
    println(io, "$( "(" * prod( [string(round(x; sigdigits=3))*", " for x in cp.data[1:end-2]]) * "$(string(round(cp.data[end-2]; sigdigits=3))))" ) in ⟂$(Int(cp.data[end]))-plane")
end

@inline chart(cp::CPAffineCoordinate) = Int64(cp.data[end])

function vec_to_affine(vec::SVector{N, ComplexF64}) where N  # N parameter unnecessary here, but generalizes to ℝP^{N-1}
    plane_idx = argmax(norm.(vec)) # Choose projection coordinate
    other_idxs = remaining_indices(plane_idx, N)   
    denom = vec[plane_idx] 
    CPAffineCoordinate{N}(ntuple(i -> i < N ? vec[other_idxs[i]]/denom : ComplexF64(plane_idx), N)) 
end

function affine_to_vec(cp::CPAffineCoordinate{N}) where N
    c = chart(cp)
    r2 = mapreduce(x -> x^2, +, cp.data[1:end-1])
    SVector{N,ComplexF64}(
        ntuple(N) do i
            if i < c
                cp[i]/sqrt(1+r2)
            elseif i == c
                1/sqrt(1+r2)
            else
                cp[i-1]/sqrt(1+r2)
            end
        end
    ) |> normalize
end

@inline remap_affine(cp) = cp |> affine_to_vec |> vec_to_affine

function affine_jac!(jac, cp::CPAffineCoordinate{N}) where N
    δ(i,j) = i == j ? 1 : 0
    Nidx(i, c) = i > c ? i - 1 : i

    c = chart(cp)
    r2 = mapreduce(x -> x^2, +, cp.data[1:end-1])
    for j in 1:N, k in 1:N-1
        jac[k,j] = if j == c
            -cp[k]/((1+r2)^(3/2))
        else
            (1/sqrt(1+r2)) * (δ(Nidx(j,c),k) - cp[k]/(sqrt(1+r2)))
        end
    end
end

################################################################################
# Optimization with ℂP^{N-1} affine coordinates
################################################################################
function optim_energy_affine(affine_coords, sys::System{N}) where N
    affine_coords = reinterpret(reshape, CPAffineCoordinate{N}, affine_coords)
    for site in all_sites(sys)
        set_coherent_state!(sys, affine_to_vec(affine_coords[site]), site)
    end
    return energy(sys)
end

function optim_gradient_affine!(buf, affine_coords, sys::System{N}) where N
    affine_coords = reinterpret(reshape, CPAffineCoordinate{N}, affine_coords)
    Hgrad = reinterpret(reshape, SVector{N, ComplexF64}, buf)

    # Calculate gradient of energy in original coordinates
    for site in all_sites(sys)
        set_coherent_state!(sys, affine_to_vec(affine_coords[site]), site)
    end
    Sunny.set_complex_forces!(Hgrad, sys.coherents, sys)

    # Calculate gradient, using Jacobian for affinegraphic coordinates 
    jac = zeros(ComplexF64, N, N)  # Maybe write function to apply matrix multiply and avoid allocation
    for site in all_sites(sys)
        affine_jac!(jac, affine_coords[site])
        Hgrad[site] = -jac * Hgrad[site]  # Note Optim expects ∇, `set_forces!` gives -∇
    end
end

function minimize_energy_affine!(sys::System{N}; method=Optim.LBFGS, maxiters = 100, kwargs...) where N
    f(affine_spins) = optim_energy_affine(affine_spins, sys)
    g!(G, affine_spins) = optim_gradient_affine!(G, affine_spins, sys)
    # fg!(energy_spins, G, x) = affine_energy_grad!(affine_spins, G, sys)


    # Quick test if any coordinates every go to infinity 
    # niters = 0
    # success = false
    # while niters < maxiters
    #     niters += 1
    #     affine_spins = map(vec -> vec_to_affine(vec), sys.dipoles) 
    #     affine_spins = Array(reinterpret(reshape, Float64, affine_spins))
    #     # affine_spins = map!(vec -> vec_to_affine(vec), affine_spins, sys.dipoles) 
    #     optout = Optim.optimize(f, g!, affine_spins, method(), options)
    #     if optout.g_converged
    #         println("We got there.")
    #         return true
    #     end
    #     # println("Maximum: ", maximum(affine_spins))
    #     # println(affine_spins)
    #     if maximum(affine_spins) > 5.0
    #         println("Remapping coords...")
    #         affine_spins = reinterpret(reshape, affinegraphicPoint{3}, affine_spins)
    #         @. affine_spins = remap_affine(affine_spins)
    #     end
    #     println("Trying again...")
    # end
    # return success

    affine_spins = map(vec -> vec_to_affine(vec), sys.coherents) 
    affine_spins = Array(reinterpret(reshape, ComplexF64, affine_spins))
    options = Optim.Options(iterations=1000, kwargs...)
    Optim.optimize(f, g!, affine_spins, method(), options)
end


################################################################################
# Complex projective representation of CP^N-1 vector 
################################################################################

using Sunny, StaticArrays, LinearAlgebra

function vec_to_stereo(u, z)
    # (I - z*z')*u/(1 - √(u'*z*z'*u))

    zz = z*z'
    (I - zz)*u/(1 - sqrt(u'*zz*u))

    # (I - z*z')*u/(1 - abs(z'*u))
end

function stereo_to_vec(v, z)
    v2 = v' * v
    denom = 1 + v2

    (2v + (v2-1)*z)/ denom
end

function spin_expectations(ψ::SVector{N, ComplexF64}) where N
    S = Sunny.spin_matrices(N)
    real.((ψ' * S[1] * ψ, ψ' * S[2] * ψ, ψ' * S[3] * ψ))
end

# Tests
function there_and_back(u::SVector{N, ComplexF64}, z::SVector{N, ComplexF64}; θ=nothing) where N
    # u1 = vec_to_stereo(u,z) |> v -> stereo_to_vec(v, z)
    v = vec_to_stereo(u, z) 
    u1 = stereo_to_vec(v, z)

    println("Norm of intermediate: ", norm(v))
    println("Reconstruction error: ", norm(u - u1))

    display(spin_expectations(u))
    display(spin_expectations(u1))
    u1 * exp(-im*angle(u1[1]))
end

function there_and_back2(u::SVector{N, ComplexF64}) where N
    rv = rand(SVector{N, ComplexF64})
    ov = rv - (u'*rv)*u
    z = normalize(ov)


    # u1 = vec_to_stereo(u,z) |> v -> stereo_to_vec(v, z)
    v = vec_to_stereo(u, z) 
    u1 = stereo_to_vec(v, z)

    println("Norm of intermediate: ", norm(v))
    println("Reconstruction error: ", norm(u - u1))

    display(spin_expectations(u))
    display(spin_expectations(u1))
    u1 * exp(-im*angle(u1[1]))
end
