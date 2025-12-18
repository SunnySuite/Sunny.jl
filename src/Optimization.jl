struct MinimizationResult
    converged :: Bool
    niters :: Int
    energy :: Float64
    data :: Optim.OptimizationResults

    function MinimizationResult(niters, res)
        return new(Optim.converged(res), niters, Optim.minimum(res), res)
    end
end

function Base.show(io::IO, ::MIME"text/plain", res::MinimizationResult)
    (; converged, niters, data) = res
    Δx = number_to_simple_string(Optim.x_abschange(data); digits=3)
    g_res = number_to_simple_string(Optim.g_residual(data); digits=3)
    if converged
        print(io, "Converged in $niters iterations")
    else
        print(io, "Non-converged after $niters iterations: |Δx|=$Δx, |∂E/∂x|=$g_res")
    end
end

function Base.show(io::IO, res::MinimizationResult)
    print(io, "MinimizationResult(")
    show(io, MIME("text/plain"), res)
    print(io, ")")
end


# Manifold where the first nspins are normalized spheres of dimension
# (nelems-1). Any remaining data is interpreted in the Euclidean metric.
struct SpinManifold <: Optim.Manifold
    nelems :: Int  # Number of scalar components in each spin
    nspins :: Int  # Number of spins to normalize
end

optim_spin_view(sm::SpinManifold, x) = view(x, 1:(sm.nelems*sm.nspins))

function Optim.retract!(sm::SpinManifold, x)
    x′ = reshape(optim_spin_view(sm, x), sm.nelems, :)
    for j in 1:sm.nspins
        xj = view(x′, :, j)
        xj ./= norm(xj)
    end
    return x
end

function Optim.project_tangent!(sm::SpinManifold, g, x)
    x = reshape(optim_spin_view(sm, x), sm.nelems, :)
    g = reshape(optim_spin_view(sm, g), sm.nelems, :)
    for j in 1:sm.nspins
        xj = view(x, :, j)
        gj = view(g, :, j)
        gj .= gj .- xj .* ((xj' * gj) / norm2(xj))
    end
    return nothing
end

function optimize_with_restarts(; calc_f, calc_g!, x, method, maxiters, options_args)
    niters = 0
    while true
        options = Optim.Options(; iterations=maxiters-niters, options_args...)
        res = Optim.optimize(calc_f, calc_g!, x, method, options)
        x = Optim.minimizer(res)
        niters += Optim.iterations(res)
        if Optim.converged(res) || niters >= maxiters
            return (; res, niters)
        end
    end
end


"""
    minimize_energy!(sys::System; maxiters=1000, jitter=1e-8, kwargs...)

Optimizes the spin configuration in `sys` to find a local minimum of the energy.
Large magnetic cells will be slower to converge; increase `maxiters` as needed.
Prior to optimization, each spin will be randomly perturbed with the
dimensionless `jitter` magnitude, which can be useful to break accidental
symmetries. Any remaining `kwargs` will be included in the `Options` constructor
of the [Optim.jl package](https://github.com/JuliaNLSolvers/Optim.jl).

Returns an object that can be inspected for optimization statistics.

!!! tip "Escaping local minima"  
    To search for the global energy minimum, a simple strategy is to repeatedly
    call [`randomize_spins!`](@ref) and then `minimize_energy!`. The example
    below keeps the minimum energy spin state found within 100 runs.

    ```julia
    tmp_sys = clone_system(sys)
    for i in 1:100
        randomize_spins!(sys)
        minimize_energy!(sys)
        if energy(sys) < energy(tmp_sys)
            copy_spins!(tmp_sys, sys)
        end
    end
    copy_spins!(sys, tmp_sys)
    ```
"""
function minimize_energy!(sys::System{N}; maxiters=1000, jitter=1e-8, kwargs...) where N
    # Small perturbation to destabilize an accidental stationary point (local
    # maximum or saddle).
    perturb_spins!(sys, jitter)

    # Optimization variables are normalized spins or coherent states. In case of
    # a vacancy, use an arbitrary representative on the sphere: [1, 1, …] / √N.
    normalize_or_fallback(x) = normalize(iszero(x) ? one.(x) : x)
    x = if iszero(N)
        collect(vec(reinterpret(Float64, normalize_or_fallback.(sys.dipoles))))
    else
        collect(vec(reinterpret(ComplexF64, normalize_or_fallback.(sys.coherents))))
    end

    # Load spins into system
    function load_spins!(x)
        if iszero(N)
            x = reinterpret(Vec3, x)
            for (i, site) in enumerate(eachsite(sys))
                set_dipole!(sys, x[i], site)
            end
        else
            x = reinterpret(CVec{N}, x)
            for (i, site) in enumerate(eachsite(sys))
                set_coherent!(sys, x[i], site)
            end
        end
    end

    # Energy and gradient callback functions
    function calc_f(x)
        load_spins!(x)
        return energy(sys)
    end
    function calc_g!(g, x)
        load_spins!(x)
        if iszero(N)
            g = reshape(reinterpret(Vec3, g), size(sys.dipoles))
            set_energy_grad_dipoles!(g, sys.dipoles, sys)
            @. g *= norm(sys.dipoles)   # Sensitivity to change in unit spin
        else
            g = reshape(reinterpret(CVec{N}, g), size(sys.coherents))
            set_energy_grad_coherents!(g, sys.coherents, sys)
            @. g *= norm(sys.coherents) # Sensitivity to change in unit ket
        end
        return nothing
    end

    # Disable check on the energy f, because we require high precision in the
    # dimensionless spin variables x. The checks x_abstol and g_abstol are in
    # the p=Inf norm (largest vector component).
    x_abstol = 1e-12
    g_abstol = 1e-12 * characteristic_energy_scale(sys)
    manifold = SpinManifold(iszero(N) ? 3 : N, length(eachsite(sys)))
    method = Optim.ConjugateGradient(; alphaguess=LineSearches.InitialHagerZhang(; αmax=10.0), manifold)
    options_args = (; g_abstol, x_abstol, x_reltol=NaN, f_reltol=NaN, f_abstol=NaN, kwargs...)
    (; res, niters) = optimize_with_restarts(; calc_f, calc_g!, x, method, maxiters, options_args)

    load_spins!(Optim.minimizer(res))
    mr = MinimizationResult(niters, res)
    if !mr.converged
        @warn repr("text/plain", mr)
    end
    return mr
end
