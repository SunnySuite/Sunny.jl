"""
    Langevin(Δt::Float64; λ::Float64, kT::Float64)

Projected Heun integration scheme with noise and damping.
Use with the `step!` function to evolve a `System` forward by a time step of `Δt`:

step!(sys::System, integrator::Langevin)

If `kT > 0`, this will simulate dynamics in the presence of a thermal bath. `λ` is an
empirical parameter that determines the strength of coupling to the thermal bath and
sets a time scale for decorrelation, `1/λ`. Both keyword parameters are required.
"""
mutable struct Langevin
    Δt  :: Float64
    λ   :: Float64
    kT  :: Float64
end

Langevin(Δt; λ, kT) = Langevin(Δt, λ, kT)


"""
    ImplicitMidpoint(Δt::Float64; atol=1e-12) where N

Energy-conserving integrator for simulating dynamics without damping or noise.
Use with the `step!` function to evolve a `System` forward by a time step of `Δt`:

step!(sys::System, integrator::ImplicitMidpoint)

The above function will use the spherical midpoint integration scheme for dipole systems
and the Schrodinger midpoint integration scheme for SU(N) spin systems.
Both integration schemes are symplectic (energy-conserving) and are appropriate for
simulating dissipationless dynamics over long periods of time.
"""
mutable struct ImplicitMidpoint
    Δt   :: Float64
    atol :: Float64
end

ImplicitMidpoint(Δt; atol=1e-12) = ImplicitMidpoint(Δt, atol)



################################################################################
# Dipole integration
################################################################################

@inline rhs_dipole(s, B) = -s × B
@inline rhs_dipole(s, B, λ) = -s × (B + λ * (s × B))

"""
    step!(sys::System, dynamics)

Advance the spin configuration one dynamical time-step. The `dynamics` object
may be a continuous spin dynamics, such as [`Langevin`](@ref) or
[`ImplicitMidpoint`](@ref), or it may be a discrete Monte Carlo sampling scheme
such as [`LocalSampler`](@ref).
"""
function step! end

function step!(sys::System{0}, integrator::Langevin)
    (B, s₁, f₁, r₁, ξ) = get_dipole_buffers(sys, 5)
    (; kT, λ, Δt) = integrator
    s = sys.dipoles

    randn!(sys.rng, ξ)
    ξ .*= √(2λ*kT)

    # Euler step
    set_forces!(B, s, sys)
    @. f₁ = rhs_dipole(s, B, λ)
    @. r₁ = rhs_dipole(s, ξ)   # note absence of λ argument -- noise only appears once in rhs.
    @. s₁ = s + Δt * f₁ + √Δt * r₁

    # Corrector step
    set_forces!(B, s₁, sys)
    @. s = s + 0.5 * Δt * (f₁ + rhs_dipole(s₁, B, λ)) + 0.5 * √Δt * (r₁ + rhs_dipole(s₁, ξ))
    @. s = normalize_dipole(s, sys.κs)
    nothing
end

# The spherical midpoint method, Phys. Rev. E 89, 061301(R) (2014)
# Integrates ds/dt = s × ∂E/∂s one timestep s → s′ via implicit equations
#   s̄ = (s′ + s) / 2
#   ŝ = s̄ / |s̄|
#   (s′ - s)/Δt = 2(s̄ - s)/Δt = - ŝ × B,
# where B = -∂E/∂ŝ.
function step!(sys::System{0}, integrator::ImplicitMidpoint)
    s = sys.dipoles
    (; Δt, atol) = integrator

    (B, s̄, ŝ, s̄′) = get_dipole_buffers(sys, 4)
    
    # Initial guess for midpoint
    @. s̄ = s

    max_iters = 100
    for _ in 1:max_iters
        # Integration step for current best guess of midpoint s̄. Produces
        # improved midpoint estimator s̄′.
        @. ŝ = normalize_dipole(s̄, sys.κs)
        set_forces!(B, ŝ, sys)
        @. s̄′ = s + 0.5 * Δt * rhs_dipole(ŝ, B)

        # If converged, then we can return
        if isapprox(s̄, s̄′, atol=atol* √length(s̄))
            # Normalization here should not be necessary in principle, but it
            # could be useful in practice for finite `atol`.
            @. s = normalize_dipole(2*s̄′ - s, sys.κs)
            return
        end

        @. s̄ = s̄′
    end

    error("Spherical midpoint method failed to converge to tolerance $atol after $max_iters iterations.")
end


################################################################################
# SU(N) integration
################################################################################


@inline function proj(a::T, Z::T)  where T <: CVec
    (a - ((Z' * a) * Z))  
end


# Returns (Λ + B⋅S) Z
@generated function mul_spin_matrices(Λ, B::Sunny.Vec3, Z::Sunny.CVec{N}) where N
    S = Sunny.spin_matrices(N)
    out = map(1:N) do i
        out_i = map(1:N) do j
            terms = Any[:(Λ[$i,$j])]
            for α = 1:3
                S_αij = S[α][i,j]
                if !iszero(S_αij)
                    push!(terms, :(B[$α] * $S_αij))
                end
            end
            :(+($(terms...)) * Z[$j])
        end
        :(+($(out_i...)))
    end
    return :(CVec{$N}($(out...)))
end


function rhs_langevin!(ΔZ::Array{CVec{N}, 4}, Z::Array{CVec{N}, 4}, ξ::Array{CVec{N}, 4},
                        B::Array{Vec3, 4}, integrator::Langevin, sys::System{N}) where N
    (; kT, λ, Δt) = integrator

    @. sys.dipoles = expected_spin(Z) # temporarily desyncs dipoles and coherents
    set_forces!(B, sys.dipoles, sys)

    if is_homogeneous(sys)
        ints = interactions_homog(sys)
        for idx in all_sites(sys)
            Λ = ints[idx[4]].aniso.matrep
            HZ = mul_spin_matrices(Λ, -B[idx], Z[idx]) # HZ = (Λ - B⋅s) Z
            ΔZ′ = -im*√(2*Δt*kT*λ)*ξ[idx] - Δt*(im+λ)*HZ
            ΔZ[idx] = proj(ΔZ′, Z[idx])
        end 
    else
        ints = interactions_inhomog(sys)
        for idx in all_sites(sys)
            Λ = ints[idx].aniso.matrep
            HZ = mul_spin_matrices(Λ, -B[idx], Z[idx]) # HZ = (Λ - B⋅s) Z
            ΔZ′ = -im*√(2*Δt*kT*λ)*ξ[idx] - Δt*(im+λ)*HZ
            ΔZ[idx] = proj(ΔZ′, Z[idx])
        end 
    end
    nothing
end


function step!(sys::System{N}, integrator::Langevin) where N
    (Z′, ΔZ₁, ΔZ₂, ξ) = get_coherent_buffers(sys, 4)
    B = get_dipole_buffers(sys, 1) |> only
    Z = sys.coherents

    randn!(sys.rng, ξ)

    # Prediction
    rhs_langevin!(ΔZ₁, Z, ξ, B, integrator, sys)
    @. Z′ = normalize_ket(Z + ΔZ₁, sys.κs)

    # Correction
    rhs_langevin!(ΔZ₂, Z′, ξ, B, integrator, sys)
    @. Z = normalize_ket(Z + (ΔZ₁+ΔZ₂)/2, sys.κs)

    # Coordinate dipole data
    @. sys.dipoles = expected_spin(Z)
end


function rhs_ll!(ΔZ, Z, B, integrator, sys)
    (; Δt) = integrator

    @. sys.dipoles = expected_spin(Z) # temporarily desyncs dipoles and coherents
    set_forces!(B, sys.dipoles, sys)

    if is_homogeneous(sys)
        ints = interactions_homog(sys)
        for idx in all_sites(sys)
            Λ = ints[idx[4]].aniso.matrep
            HZ = mul_spin_matrices(Λ, -B[idx], Z[idx]) # HZ = (Λ - B⋅s) Z
            ΔZ[idx] = - Δt*im*HZ
        end 
    else
        ints = interactions_inhomog(sys)
        for idx in all_sites(sys)
            Λ = ints[idx].aniso.matrep
            HZ = mul_spin_matrices(Λ, -B[idx], Z[idx]) # HZ = (Λ - B⋅s) Z
            ΔZ[idx] = - Δt*im*HZ
        end 
    end
end

# Implicit Midpoint Method applied to the nonlinear Schrödinger dynamics, as
# proposed in Phys. Rev. B 106, 054423 (2022). Integrates dZ/dt = - i ℌ(Z) Z one
# timestep Z → Z′ via the implicit equation
#
#   (Z′-Z)/Δt = - i ℌ(Z̄) Z, where Z̄ = (Z+Z′)/2
#
function step!(sys::System{N}, integrator::ImplicitMidpoint; max_iters=100) where N
    (; atol) = integrator
    (ΔZ, Z̄, Z′, Z″) = get_coherent_buffers(sys, 4)
    B = get_dipole_buffers(sys, 1) |> only
    Z = sys.coherents

    @. Z′ = Z 
    @. Z″ = Z 

    for _ in 1:max_iters
        @. Z̄ = (Z + Z′)/2

        rhs_ll!(ΔZ, Z̄, B, integrator, sys)

        @. Z″ = Z + ΔZ

        if isapprox(Z′, Z″, atol=atol*√length(Z′))
            @. Z = normalize_ket(Z″, sys.κs)
            @. sys.dipoles = expected_spin(Z)
            return
        end

        Z′, Z″ = Z″, Z′
    end

    error("Schrödinger midpoint method failed to converge in $max_iters iterations.")
end
