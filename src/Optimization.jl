# Converts unnormalized representation of coherent state, α, to standard dipole
# or normalized complex vector representation.
function projective_to_conventional(α::SVector{N, T}, z::SVector{N, T}) where {N, T}
    v = (I - z*z')*α  # Automatically infers that `I` should be an SMatrix (no allocations)
    v2 = v'*v
    (2v + (1-v2)*z) / (1+v2)  # Guaranteed to be normalized
end

# Calculate du(α)/dα and apply to vec.
@inline function apply_projective_jacobian(vec, α::SVector{N, T}, z::SVector{N, T}) where {N, T}
    P = (I - z*z')
    v = P*α
    dv_dα = P 
    dv2_dα = 2*(α' - 2*(α'*z)*z')
    v2 = v' * v 

    jac = (1/(1+v2)) * ((2dv_dα - (z*dv2_dα)') - (1/(1+v2)) * dv2_dα' * ((2v + (1-v2)*z)'))

    return jac * vec
end

# Calculate H(u(α)) for SU(N) mode 
function optim_energy(αs, zs, sys::System{N}) where N
    αs = reinterpret(reshape, SVector{N, ComplexF64}, αs)
    for site in all_sites(sys)
        set_coherent_state!(sys, projective_to_conventional(αs[site], zs[site]), site)
    end
    return energy(sys)
end

# Calculate H(u(α)) for dipole mode 
function optim_energy(αs, zs, sys::System{0})
    αs = reinterpret(reshape, Vec3, αs)
    for site in all_sites(sys)
        polarize_spin!(sys, projective_to_conventional(αs[site], zs[site]), site)
    end
    return energy(sys)
end

# Calculate dH(u(α))/dα for coherent states
function optim_gradient!(buf, αs, zs, B, sys::System{N}) where N
    αs = reinterpret(reshape, SVector{N, ComplexF64}, αs)
    Hgrad = reinterpret(reshape, SVector{N, ComplexF64}, buf)

    for site in all_sites(sys)
        set_coherent_state!(sys, projective_to_conventional(αs[site], zs[site]), site)
    end
    Sunny.set_forces!(B, sys.dipoles, sys)
    Sunny.set_complex_forces!(Hgrad, B, sys.coherents, sys)

    for site in all_sites(sys)
        Hgrad[site] = apply_projective_jacobian(Hgrad[site], αs[site], zs[site])
    end
end

# Calculate dH(u(α))/dα for dipoles
function optim_gradient!(buf, αs, zs, _, sys::System{0})
    αs = reinterpret(reshape, Vec3, αs)
    Hgrad = reinterpret(reshape, Vec3, buf)

    for site in all_sites(sys)
        polarize_spin!(sys, projective_to_conventional(αs[site], zs[site]), site)
    end
    Sunny.set_forces!(Hgrad, sys.dipoles, sys)

    for site in all_sites(sys)
        Hgrad[site] = -apply_projective_jacobian(Hgrad[site], αs[site], zs[site])
    end
end

"""
    minimize_energy!(sys::System{N}; method=Optim.LBFGS, maxiters = 10000, kwargs...) where N
"""
function minimize_energy!(sys::System{N}; method=Optim.LBFGS, maxiters = 10000, kwargs...) where N
    numbertype = N == 0 ? Float64 : ComplexF64
    buffer = N == 0 ? sys.dipoles : sys.coherents
    B = N == 0 ? nothing : get_dipole_buffers(sys, 1) |> only

    zs = copy(buffer)
    αs = zeros(numbertype, (N == 0 ? 3 : N, size(buffer)...))

    f(proj_coords) = optim_energy(proj_coords, zs, sys)
    g!(G, proj_coords) = optim_gradient!(G, proj_coords, zs, B, sys)

    options = Optim.Options(iterations=maxiters, kwargs...)
    Optim.optimize(f, g!, αs, method(), options)
end