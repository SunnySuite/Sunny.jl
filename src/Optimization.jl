# The following four helper functions allow for more code resuse since the
# projective parameterization is formally the same for both dipoles and coherent
# states (the element types are just real in the first case, compelx in the
# second).
function set_forces_optim!(∇H, ∇H_dip, sys::System{N}) where {N}
    Sunny.set_energy_grad_dipoles!(∇H_dip, sys.dipoles, sys)
    Sunny.set_energy_grad_coherents!(∇H, ∇H_dip, sys.coherents, sys)
end

function set_forces_optim!(∇H, _, sys::System{0}) 
    Sunny.set_energy_grad_dipoles!(∇H, sys.dipoles, sys)
end

@inline function set_spin_optim!(sys::System{N}, α, z, site) where N
    set_coherent_state!(sys, projective_to_conventional(α, z), site)
end

@inline function set_spin_optim!(sys::System{0}, α, z, site)
    polarize_spin!(sys, projective_to_conventional(α, z), site)
end

# Converts unnormalized representation of coherent state, α, to standard dipole
# or normalized complex vector representation.
function projective_to_conventional(α, z)
    v = (I - z*z')*α
    v2 = v'*v
    return (2v + (1-v2)*z) / (1+v2)  # Guaranteed to be normalized
end

# Calculate du(α)/dα and apply to `vec`. u(α) = (2v + (1-v²)z)/(1+v²) with v =
# (1-zz†)α. Won't allocate if all inputs are StaticArrays.
@inline function apply_projective_jacobian(vec, α, z)
    P = (I - z*z')
    v = P*α
    v2 = v'*v 

    dv_dα = P 
    dv2_dα = 2*(α' - 2*(α'*z)*z')
    jac = (1/(1+v2)) * ((2dv_dα - (z*dv2_dα)') - (1/(1+v2)) * dv2_dα' * ((2v + (1-v2)*z)'))

    return jac * vec
end

# Calculate H(u(α)) 
function optim_energy(αs, zs, sys::System{N})::Float64 where N
    T = N == 0 ? Vec3 : CVec{N} 
    αs = reinterpret(reshape, T, αs)
    for site in all_sites(sys)
        set_spin_optim!(sys, αs[site], zs[site], site)
    end
    return energy(sys) # Note: `Sunny.energy` seems to allocate and is type-unstable
end

# Non-allocating check for largest unnormalized coordinate.
function maxnorm(αs)
    max = 0.0
    for α in αs
        magnitude = norm(α) 
        max = magnitude > max ? magnitude : max
    end
    return max
end

# Calculate dH(u(α))/dα 
function optim_gradient!(buf, αs, zs, B, sys::System{N}, halted, quickmode=false) where N
    T = N == 0 ? Vec3 : CVec{N} 
    αs = reinterpret(reshape, T, αs)
    Hgrad = reinterpret(reshape, T, buf)

    # Check if any coordinate has gone adrift and signal need to reset if necessary
    if !quickmode
        maxdist = maxnorm(αs)
        if maxdist > 1.5     # 1.5 found empirically, works well for both dipole and SU(3) 
            Hgrad .*= 0      # Trick Optim.jl to stop by setting gradient to 0 (triggers convergence tests) 
            halted[] = true  # Let main loop know we haven't really converged
            return
        end
    end

    # Calculate gradient of energy with respect to α
    for site in all_sites(sys)
        set_spin_optim!(sys, αs[site], zs[site], site)
    end

    set_forces_optim!(Hgrad, B, sys)

    for site in all_sites(sys)
        Hgrad[site] = apply_projective_jacobian(Hgrad[site], αs[site], zs[site])
    end
end

# Quick "touchup" optimization that assumes the system is already near the
# ground state. Never changes the parameterization of coherent states or
# dipoles. For internal use when setting up a spin wave calculation.
function minimize_energy_touchup!(sys::System{N}; method=Optim.LBFGS, maxiters = 40, kwargs...) where N
    numbertype = N == 0 ? Float64 : ComplexF64
    buffer = N == 0 ? sys.dipoles : sys.coherents
    B = N == 0 ? nothing : get_dipole_buffers(sys, 1) |> only

    zs = copy(buffer)
    αs = zeros(numbertype, (N == 0 ? 3 : N, size(buffer)...))
    halted = Ref(false)

    f(proj_coords) = optim_energy(proj_coords, zs, sys)
    g!(G, proj_coords) = optim_gradient!(G, proj_coords, zs, B, sys, halted, true) # true skips coordinate drifting test

    options = Optim.Options(iterations=maxiters, kwargs...)
    output = Optim.optimize(f, g!, αs, method(), options)
    success = Optim.converged(output)
    if !success
        @warn "Optimization failed to converge within $(output.iterations) iterations. `System` not in ground state and spin wave calculations may fail."
    end

    return success 
end

"""
    minimize_energy!(sys::System{N}; method=Optim.LBFGS, maxiters = 1000, kwargs...) where N

Optimizes the spin configuration in `sys` to minimize energy. Any method from
Optim.jl that accepts only a gradient may be used by setting the `method`
keyword. Defaults to LBFGS.
"""
function minimize_energy!(sys::System{N}; method=Optim.LBFGS, maxiters = 1000, kwargs...) where N

    # Set up type and buffer information depending on system type
    numbertype = N == 0 ? Float64 : ComplexF64
    buffer = N == 0 ? sys.dipoles : sys.coherents
    B = N == 0 ? nothing : get_dipole_buffers(sys, 1) |> only

    # Allocate buffers for optimization
    zs = copy(buffer)
    αs = zeros(numbertype, (N == 0 ? 3 : N, size(buffer)...))
    halted = Ref(false)

    # Set up energy and gradient functions using closures to pass "constants" 
    f(proj_coords) = optim_energy(proj_coords, zs, sys)
    g!(G, proj_coords) = optim_gradient!(G, proj_coords, zs, B, sys, halted, false)

    # Perform optimization, resetting parameterization of coherent states as necessary
    options = Optim.Options(iterations=maxiters, kwargs...)
    totaliters = 0
    success = false
    while totaliters < maxiters
        output = Optim.optimize(f, g!, αs, method(), options)
        if halted[]
            zs .= buffer  # Reset parameterization based on current state
            αs .*= 0      # Reset unnormalized coordinates
            halted[] = false
        else
            if Optim.converged(output)  # Convergence report only meaningful if not halted
                success = true
                break  
            end
        end
        totaliters += output.iterations
    end

    return success
end