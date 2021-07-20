""" Structs and methods for dynamics integrators -- attempted to stay as close as possible
     to the types of David's, but didn't make it parametric over the Float type (T)
"""

""" Standard (non-damped) Landau-Lifshitz dynamics integrates:
      dS/dt = -S × H, where H = B + ∑_neigh S_neigh
"""

abstract type Integrator end

mutable struct SpinHeunP{D, L, Db} <: Integrator
    sys :: SpinSystem{D, L, Db}
    _S₁ :: Array{Vec3, Db}
    _S₂ :: Array{Vec3, Db}
    _H₁ :: Array{Vec3, Db}
    _H₂ :: Array{Vec3, Db}
    _f₁ :: Array{Vec3, Db}
end

function SpinHeunP(sys::SpinSystem)
    return SpinHeunP(
        sys, zero(sys.sites), zero(sys.sites), zero(sys.sites),
        zero(sys.sites), zero(sys.sites),
    )
end


function evolve!(integrator::SpinHeunP, Δt::Float64)
    @unpack sys, _S₁, _S₂, _H₁, _H₂, _f₁ = integrator
    # Euler step
    field!(_H₁, sys.sites, sys.interactions)
    for idx in eachindex(sys)
        _f₁[idx] = -sys[idx] × _H₁[idx]
        _S₁[idx] = sys[idx] + Δt * _f₁[idx]
    end

    # Corrector step
    field!(_H₂, _S₁, sys.interactions)
    for idx in eachindex(sys)
        f₂ = -_S₁[idx] × _H₂[idx]
        S₂ = sys[idx] + 0.5 * Δt * (_f₁[idx] + f₂)
        S₂ = S₂ / norm(S₂)
        _S₂[idx] = S₂
    end

    # Swap buffers
    sys.sites, integrator._S₂ = integrator._S₂, sys.sites
end