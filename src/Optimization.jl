struct MinimizationResult
    converged :: Bool
    iterations :: Int
    data :: Optim.OptimizationResults
end

function Base.show(io::IO, ::MIME"text/plain", res::MinimizationResult)
    (; converged, iterations, data) = res
    Δx = number_to_simple_string(Optim.x_abschange(data); digits=3)
    g_res = number_to_simple_string(Optim.g_residual(data); digits=3)
    if converged
        print(io, "Converged in $iterations iterations: |Δx|=$Δx, |∂E/∂x|=$g_res")
    else
        print(io, "Non-converged after $iterations iterations: |Δx|=$Δx, |∂E/∂x|=$g_res")
    end
end

function Base.show(io::IO, res::MinimizationResult)
    print(io, "MinimizationResult(")
    show(io, MIME("text/plain"), res)
    print(io, ")")
end


# Returns the stereographic projection u(α) = (2v + (1-v²)n)/(1+v²), which
# involves the orthographic projection v = (1-nn̄)α. The input `n` must be
# normalized. When `α=0`, the output is `u=n`, and when `|α|→ ∞` the output is
# `u=-n`. In all cases, `|u|=1`.
function stereographic_projection(α, n)
    @assert n'*n ≈ 1 || all(isnan, n)

    v = α - n*(n'*α)              # project out component parallel to `n`
    v² = real(v'*v)
    u = (2v + (1-v²)*n) / (1+v²)  # stereographic projection
    return u
end

# Calculate the vector-Jacobian-product x̄ du(α)/dα, where
#   u(v) = (2v + (1-v²)n)/(1+v²), v(α,ᾱ)=Pα, and P=1-nn̄.
#
# From the chain rule for Wirtinger derivatives,
#   du/dα = (du/dv) (dv/dα) + (du/dv̄) (dv̄/dα) = du/dv P.
#
# In the second step, we used
#   dv/dα = P
#   dv̄/dα = conj(dv/dᾱ) = 0.
# 
# The remaining Jacobian matrix is
#   du/dv = (2-2nv̄)/(1+v²) - 2(2v+(1-v²)n)/(1+v²)² v̄
#         = c - c[(1+cb)n + cv]v̄,
# where b = (1-v²)/2 and c = 2/(1+v²).
#
# Using the above definitions, return:
#   x̄ du/dα = x̄ du/dv P
#
@inline function vjp_stereographic_projection(x̄, α, n)
    all(isnan, n) && return zero(n') # No gradient when α is fixed to zero

    @assert n'*n ≈ 1

    v = α - n*(n'*α)
    v² = real(v'*v)
    b = (1-v²)/2
    c = 2/(1+v²)
    # Perform dot products first to avoid constructing outer-product
    x̄_dudv = c*x̄' - c * (x̄' * ((1+c*b)*n + c*v)) * v'
    # Apply projection P=1-nn̄ on right
    return x̄_dudv - (x̄_dudv * n) * n'
end

# Returns v such that u = (2v + (1-v²)n)/(1+v²) and v⋅n = 0
function inverse_stereographic_projection(u, n)
    all(isnan, n) && return zero(n) # NaN values denote α = v = zero

    @assert u'*u ≈ 1

    uperp = u - n*(n'*u)
    uperp² = real(uperp' * uperp)
    s = sign(u⋅n)
    if isone(s) && uperp² < 1e-5
        c = 1/2 + uperp²/8 + uperp²*uperp²/16
    else
        c = (1 - s * sqrt(max(1 - uperp², 0))) / uperp²
    end
    return c * uperp
end

function optim_set_spins!(sys::System{0}, αs, ns)
    αs = reinterpret(reshape, Vec3, αs)
    for site in eachsite(sys)
        s = stereographic_projection(αs[site], ns[site])
        set_dipole!(sys, s, site)
    end
end
function optim_set_spins!(sys::System{N}, αs, ns) where N
    αs = reinterpret(reshape, CVec{N}, αs)
    for site in eachsite(sys)
        Z = stereographic_projection(αs[site], ns[site])
        set_coherent!(sys, Z, site)
    end
end

function optim_set_gradient!(G, sys::System{0}, αs, ns)
    (αs, G) = reinterpret.(reshape, Vec3, (αs, G))
    set_energy_grad_dipoles!(G, sys.dipoles, sys)            # G = dE/dS
    @. G *= norm(sys.dipoles)                                # G = dE/dS * dS/du = dE/du
    @. G = adjoint(vjp_stereographic_projection(G, αs, ns))  # G = dE/du du/dα = dE/dα
end
function optim_set_gradient!(G, sys::System{N}, αs, ns) where N
    (αs, G) = reinterpret.(reshape, CVec{N}, (αs, G))
    set_energy_grad_coherents!(G, sys.coherents, sys)        # G = dE/dZ
    @. G *= norm(sys.coherents)                              # G = dE/dZ * dZ/du = dE/du
    @. G = adjoint(vjp_stereographic_projection(G, αs, ns))  # G = dE/du du/dα = dE/dα
end

function minimize_energy2!(sys::System{N}; maxiters=1000, method=Optim.ConjugateGradient(),
                          subiters=10, δ=1e-8, kwargs...) where N
    # Perturbation of sufficient magnitude to "almost surely" push away from an
    # unstable stationary point (e.g. local maximum or saddle).
    perturb_spins!(sys, δ)

    # Allocate buffers for optimization:
    #   - Each `ns[site]` defines a direction for stereographic projection.
    #   - Each `αs[:,site]` will be optimized in the space orthogonal to `ns[site]`.
    if iszero(N)
        ns = normalize.(sys.dipoles)
        αs = zeros(Float64, 3, size(sys.dipoles)...)
    else
        ns = normalize.(sys.coherents)
        αs = zeros(ComplexF64, N, size(sys.coherents)...)
    end

    # Functions to calculate energy and gradient for the state `αs`
    function f(αs)
        optim_set_spins!(sys, αs, ns)
        return energy(sys)
    end
    function g!(G, αs)
        optim_set_spins!(sys, αs, ns)
        optim_set_gradient!(G, sys, αs, ns)
    end

    # Repeatedly optimize using a small number (`subiters`) of steps, within
    # which the stereographic projection axes are fixed. Because we require
    # high-precision in the spin variables (x), disable check on the energy (f),
    # which is much lower precision. Disable check on x_reltol because x=[zeros]
    # is a valid configuration (all spins aligned with the stereographic
    # projection axes). Optim interprets x_abstol and g_abstol in the p=Inf norm
    # (largest vector component). Note that x is dimensionless and g=dE/dx has
    # energy units.
    x_abstol = 1e-12
    g_abstol = 1e-12 * characteristic_energy_scale(sys)
    options = Optim.Options(; iterations=subiters, g_abstol, x_abstol, x_reltol=NaN, f_reltol=NaN, f_abstol=NaN, kwargs...)
    local res
    for iter in 1 : div(maxiters, subiters, RoundUp)
        res = Optim.optimize(f, g!, αs, method, options)

        if Optim.converged(res)
            cnt = (iter-1)*subiters + res.iterations
            return MinimizationResult(true, cnt, res)
        end

        # Reset stereographic projection based on current state
        ns .= normalize.(iszero(N) ? sys.dipoles : sys.coherents)
        αs .*= 0
    end

    mr = MinimizationResult(false, maxiters, res)
    @warn repr("text/plain", mr)
    return mr
end





struct SpinManifold <: Optim.Manifold
    sys :: System
end

function rescale_to_magnitude(S, S0)
    iszero(S0) ? zero(S) : normalize(S) * S0
end
function optim_retract_aux!(sys::System{0}, x)
    x = reinterpret(reshape, Vec3, x)
    @. x = rescale_to_magnitude(x, sys.κs)
end
function optim_retract_aux!(sys::System{N}, x) where N
    x = reinterpret(reshape, CVec{N}, x)
    @. x = rescale_to_magnitude(x, sqrt(sys.κs))
end
function Optim.retract!(sm::SpinManifold, x)
    optim_retract_aux!(sm.sys, x)
    return x
end

function optim_project_tangent_aux!(sys::System{0}, g, x)
    g = reinterpret(reshape, Vec3, g)
    x = reinterpret(reshape, Vec3, x)
    @. g = proj(g, x)
end
function optim_project_tangent_aux!(sys::System{N}, g, x) where N
    g = reinterpret(reshape, CVec{N}, g)
    x = reinterpret(reshape, CVec{N}, x)
    @. g = proj(g, x)
end
function Optim.project_tangent!(sm::SpinManifold, g, x)
    optim_project_tangent_aux!(sm.sys, g, x)
    return g
end

function characteristic_vector_scale(sys::System{0})
    return sqrt(Statistics.mean(κ^2 for κ in sys.κs))
end
function characteristic_vector_scale(sys::System{N}) where N
    return sqrt(Statistics.mean(κ for κ in sys.κs))
end

function optimize_with_restarts(; calc_f, calc_g!, x, method, maxiters, options_args)
    iters = 0
    while true
        options = Optim.Options(; iterations=maxiters-iters, options_args...)
        res = Optim.optimize(calc_f, calc_g!, x, method, options)
        x = Optim.minimizer(res)
        iters += Optim.iterations(res)
        if iters >= maxiters || Optim.termination_code(res) != Optim.TerminationCode.FailedLinesearch
            return (res, iters)
        end
    end
end

"""
    minimize_energy!(sys::System; maxiters=1000, kwargs...)

Optimizes the spin configuration in `sys` to minimize energy. A total of
`maxiters` iterations will be attempted. Any remaining `kwargs` will be included
in the `Options` constructor of the [Optim.jl
package](https://github.com/JuliaNLSolvers/Optim.jl)

Convergence status is stored in the field `ret.converged` of the return value
`ret`. Additional optimization statistics are stored in the field `ret.data`.
"""
function minimize_energy!(sys::System{N}; maxiters=1000, δ=1e-8, kwargs...) where N
    # Small perturbation to destabilize an accidental stationary point (local
    # maximum or saddle).
    perturb_spins!(sys, δ)

    # Optimization variables are normalized spins or coherent states
    if iszero(N)
        x = collect(reinterpret(reshape, Float64, sys.dipoles))
    else
        x = collect(reinterpret(reshape, ComplexF64, sys.coherents))
    end

    function load_spin_state!(sys::System{0}, x)
        x = reinterpret(reshape, Vec3, x)
        for site in eachsite(sys)
            set_dipole!(sys, x[site], site)
        end
    end
    function load_spin_state!(sys::System{N}, x) where N
        x = reinterpret(reshape, CVec{N}, x)
        for site in eachsite(sys)
            set_coherent!(sys, x[site], site)
        end
    end

    # Energy and gradient callback functions
    function calc_f(x)
        load_spin_state!(sys, x)
        return energy(sys)
    end
    function calc_g!(g, x)
        load_spin_state!(sys, x)
        if iszero(N)
            g = reinterpret(reshape, Vec3, g)
            set_energy_grad_dipoles!(g, sys.dipoles, sys)
            @. g = proj(g, sys.dipoles)
        else
            g = reinterpret(reshape, CVec{N}, g)
            set_energy_grad_coherents!(g, sys.coherents, sys)
            @. g = proj(g, sys.coherents)
        end
    end

    # Disable check on the energy f, because we require high precision in the
    # dimensionless spin variables x. The checks x_abstol and g_abstol are in
    # the p=Inf norm (largest vector component).
    E0 = characteristic_energy_scale(sys)
    x0 = characteristic_vector_scale(sys)
    x_abstol = 1e-12 * x0
    g_abstol = 1e-12 * E0 / x0
    method = Optim.ConjugateGradient(; alphaguess=LineSearches.InitialHagerZhang(; αmax=10x0), manifold=SpinManifold(sys))
    options_args = (; g_abstol, x_abstol, x_reltol=NaN, f_reltol=NaN, f_abstol=NaN, kwargs...)
    (res, iters) = optimize_with_restarts(; calc_f, calc_g!, x, method, maxiters, options_args)

    load_spin_state!(sys, Optim.minimizer(res))
    mr = MinimizationResult(Optim.converged(res), iters, res)
    if !mr.converged
        @warn repr("text/plain", mr)
    end
    return mr
end
