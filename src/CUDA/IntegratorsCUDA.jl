""" Structs and methods for dynamics integrators on the GPU using CUDA
"""

""" Standard (non-damped) Landau-Lifshitz dynamics integrates:
      dS/dt = -S × H, where H = B + ∑_neigh S_neigh
"""

abstract type IntegratorCUDA end

"""
    HeunPCUDA(sys::SpinSystem)

Integrates Landau-Lifshitz spin dynamics using the Heun method, with a final
projection step that exactly constrains |S|=1. The method is locally second
order accurate.
"""
mutable struct HeunPCUDA{D, L, Db} <: IntegratorCUDA
    sys :: SpinSystem{D, L, Db, HamiltonianCUDA{D}}
    _S₁ :: CUDA.CuArray{Vec3, Db}
    _S₂ :: CUDA.CuArray{Vec3, Db}
    _B  :: CUDA.CuArray{Vec3, Db}
    _f₁ :: CUDA.CuArray{Vec3, Db}

    function HeunPCUDA(sys::SpinSystem{D, L, Db, HamiltonianCUDA{D}}) where {D, L, Db}
        return new{D, L, Db}(
            sys, CUDA.zero(sys.sites), CUDA.zero(sys.sites),
            CUDA.zero(sys.sites), CUDA.zero(sys.sites)
        )
    end
end

"""
    evolve!(integrator, Δt)

Performs a single integrator timestep of size Δt.
"""
function evolve!(integrator::HeunPCUDA, Δt::Float64)
    # TODO
    nothing
end