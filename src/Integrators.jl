@doc raw"""
    LangevinHeunP(kT::Float64, λ::Float64, Δt::Float64)

Projected Heun integration scheme with noise and damping.
Use with the `step!` function to evolve a `SpinSystem` forward by a time step of `Δt`:

step!(sys::SpinSystem, integrator::LangevinHeunP)

If `kT > 0`, this will simulate dynamics in the presence of a thermal bath. `λ` is an
empirical parameter that determines the strength of coupling to the thermal bath and
sets a time scale for decorrelation, `1/λ`.
"""


mutable struct LangevinHeunP
    kT  :: Float64
    λ   :: Float64
    Δt  :: Float64
end


@doc raw"""
    ImplicitMidpoint(Δt::Float64; atol=1e-12) where N

Energy-conserving integrator for simulating dynamics without damping or noise.
Use with the `step!` function to evolve a `SpinSystem` forward by a time step of `Δt`:

step!(sys::SpinSystem, integrator::ImplicitMidpoint)

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
# Convenience functions

function set_expected_spins!(dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, sys::SpinSystem) where N
    @assert N > 0
    num_sites = size(dipoles)[end]
    for site in 1:num_sites
        spin_rescaling = sys.site_infos[site].spin_rescaling
        for cell in CartesianIndices(size(dipoles)[1:3]) 
            dipoles[cell,site] = spin_rescaling * expected_spin(coherents[cell,site])
        end
    end
end

set_expected_spins!(sys::SpinSystem) = set_expected_spins!(sys.dipoles, sys.coherents, sys)


################################################################################
# Dipole integration
################################################################################
# Normalize to κ value given in site_infos. For old LL dynamics only.
function normalize!(S::Array{Vec3, 4}, sys::SpinSystem)
    for site in 1:size(S, 4)
        spin_rescaling = sys.site_infos[site].spin_rescaling
        for cell in CartesianIndices(size(S)[1:3]) 
            S[cell, site] *= spin_rescaling/norm(S[cell, site])
        end
    end
end

normalize_dipoles!(sys::SpinSystem) = normalize!(sys.dipoles, sys)

@inline rhs_dipole(S, B) = -S × B
@inline rhs_dipole(S, B, λ) = -S × (B + λ * (S × B))

@doc raw"""
    step!(sys::SpinSystem, integrator)

Advance the spin system forward one step using the parameters and integration
scheme specified by `integrator`.
"""
function step!(sys::SpinSystem{0}, integrator::LangevinHeunP)
    (B, S₁, f₁, r₁, ξ) = get_dipole_buffers(sys, 5)
    (; kT, λ, Δt) = integrator
    S, ℋ = sys.dipoles, sys.hamiltonian

    randn!(sys.rng, ξ)
    ξ .*= √(2λ*kT)

    # Euler step
    field!(B, S, ℋ)
    @. f₁ = rhs_dipole(S, B, λ)
    @. r₁ = rhs_dipole(S, ξ)   # note absence of λ argument -- noise only appears once in rhs.
    @. S₁ = S + Δt * f₁ + √Δt * r₁

    # Corrector step
    field!(B, S₁, ℋ)
    @. S = S + 0.5 * Δt * (f₁ + rhs_dipole(S₁, B, λ)) + 0.5 * √Δt * (r₁ + rhs_dipole(S₁, ξ))
    normalize!(S, sys)

    nothing
end

function step!(sys::SpinSystem{0}, integrator::ImplicitMidpoint)
    S, ℋ = sys.dipoles, sys.hamiltonian
    (; Δt, atol) = integrator
    (B, S̄, Ŝ, S̄′) = get_dipole_buffers(sys, 4)
    
    # Initial guess for midpoint
    @. S̄ = S

    max_iters = 100
    for _ in 1:max_iters
        # Integration step for current best guess of midpoint S̄. Produces
        # improved midpoint estimator _S̄′.
        @. Ŝ = S̄
        normalize!(Ŝ, sys)
        field!(B, Ŝ, ℋ)
        @. S̄′ = S + 0.5 * Δt * rhs_dipole(Ŝ, B)

        # Convergence is reached if every element of _S̄ and _S̄′ agree to
        # within tolerance.
        converged = true
        for i = 1:length(S)
            converged = converged && isapprox(S̄[i], S̄′[i]; atol)
        end

        # If converged, then we can return
        if converged
            # Normalization here should not be necessary in principle, but it
            # could be useful in practice for finite `atol`.
            @. S = 2*S̄′ - S
            normalize!(S, sys)
            return
        end

        @. S̄ = S̄′
    end

    error("Spherical midpoint method failed to converge to tolerance $atol after $max_iters iterations.")
end


################################################################################
# SU(N) integration
################################################################################
# LinearAlgebra's `normalize!` doesn't normalize to 1 on complex vectors
function normalize!(Z::Array{CVec{N}, 4}) where N
    @inbounds for i in eachindex(Z)
        Z[i] = Z[i] / √(real(Z[i]' * Z[i]))
    end
end

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
                        B::Array{Vec3, 4}, integrator::LangevinHeunP, sys::SpinSystem{N}) where N
    (; kT, λ, Δt) = integrator
    (; dipoles, hamiltonian) = sys

    set_expected_spins!(dipoles, Z, sys) 
    field!(B, dipoles, hamiltonian)

    @inbounds for i in CartesianIndices(ξ)
        Λ = view(hamiltonian.sun_aniso, :, :, i[4])
        ℌZ = mul_spin_matrices(Λ, -B[i], Z[i]) # (Λ - B⋅S) Z

        # The field κ describes an effective ket rescaling, Z → Z' = √κ Z. For
        # numerical convenience, Sunny avoids explicit ket rescaling, and
        # instead performs an equivalent rescaling expectation values, ⟨A⟩ → κ
        # ⟨A⟩. In the latter approach, the noise term in Langevin dynamics must
        # also be rescaled as ξ → 1/√κ ξ. One can see the equivalence directly
        # by writing the Langevin equation in Z' and then dividing both sides √κ
        # to get an equivalent dynamics in normalized spins. The resulting
        # dynamics samples the correct Boltzmann equilibrium involving the
        # rescaled classical Hamiltonian.
        κ = sys.site_infos[i[4]].spin_rescaling
        ΔZ′ = -im*√(2*Δt*kT*λ/κ)*ξ[i] - Δt*(im+λ)*ℌZ
        ΔZ[i] = proj(ΔZ′, Z[i])
    end 
    nothing
end


function step!(sys::SpinSystem{N}, integrator::LangevinHeunP) where N
    (Z′, ΔZ₁, ΔZ₂, ξ) = get_coherent_buffers(sys, 4)
    B = get_dipole_buffers(sys, 1) |> only
    Z = sys.coherents

    randn!(sys.rng, ξ)

    # Prediction
    rhs_langevin!(ΔZ₁, Z, ξ, B, integrator, sys)
    @. Z′ = Z + ΔZ₁
    normalize!(Z′)

    # Correction
    rhs_langevin!(ΔZ₂, Z′, ξ, B, integrator, sys)
    @. Z += (ΔZ₁ + ΔZ₂)/2
    normalize!(Z)

    # Coordinate dipole data
    set_expected_spins!(sys)
end


function rhs_ll!(ΔZ, Z, B, integrator, sys)
    (; Δt) = integrator
    (; dipoles, hamiltonian) = sys

    set_expected_spins!(dipoles, Z, sys) # temporarily de-synchs dipoles and coherents
    field!(B, dipoles, hamiltonian)

    @inbounds for idx in CartesianIndices(Z)
        Λ = view(hamiltonian.sun_aniso, :, :, idx[4])
        ℌZ = mul_spin_matrices(Λ, -B[idx], Z[idx])
        ΔZ[idx] = - Δt*im*ℌZ
    end 
end


function step!(sys::SpinSystem, integrator::ImplicitMidpoint; max_iters=100)
    (; atol) = integrator
    (ΔZ, Zb, Z′, Z″) = get_coherent_buffers(sys, 4)
    B = get_dipole_buffers(sys, 1) |> only
    Z = sys.coherents

    @. Z′ = Z 
    @. Z″ = Z 

    for _ in 1:max_iters
        @. Zb = (Z + Z′)/2

        rhs_ll!(ΔZ, Zb, B, integrator, sys)

        @. Z″ = Z + ΔZ

        if isapprox(Z′, Z″; atol)
            @. Z = Z″
            set_expected_spins!(sys)
            return
        end

        Z′, Z″ = Z″, Z′
    end

    error("SchrodingerMidpoint integrator failed to converge in $max_iters iterations.")
end


################################################################################
# Langevin Sampler 
################################################################################
@doc raw"""
    LangevinSampler(integrator::LangevinHeunP, nsteps::Int)

Creates a sampler from a Langevin integrator. `nsteps` determines how many
times `step!` is called using the integrator. `nsteps` should be selected large
enough to ensure that the state of the SpinSystem after integration
is decorrelated with respect to its initial state.
"""
mutable struct LangevinSampler <: AbstractSampler
    integrator :: LangevinHeunP
    nsteps     :: Int
end

set_temp!(sampler::LangevinSampler, kT) = sampler.integrator.kT = kT
get_temp(sampler::LangevinSampler) = sampler.integrator.kT

function sample!(sys, sampler::LangevinSampler)
    (; integrator, nsteps) = sampler
    for _ ∈ 1:nsteps
        step!(sys, integrator)
    end
    return nothing
end