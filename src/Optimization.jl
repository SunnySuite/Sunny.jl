
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
    @. G *= norm(sys.dipoles)                                # G = dE/dS * ds/du = dE/du
    @. G = adjoint(vjp_stereographic_projection(G, αs, ns))  # G = dE/du du/dα = dE/dα
end
function optim_set_gradient!(G, sys::System{N}, αs, ns) where N
    (αs, G) = reinterpret.(reshape, CVec{N}, (αs, G))
    set_energy_grad_coherents!(G, sys.coherents, sys)        # G = dE/dZ
    @. G *= norm(sys.coherents)                              # G = dE/dZ * dZ/du = dE/du
    @. G = adjoint(vjp_stereographic_projection(G, αs, ns))  # G = dE/du du/dα = dE/dα
end


"""
    minimize_energy!(sys::System{N}; maxiters=1000, method=Optim.ConjugateGradient(),
                     kwargs...) where N

Optimizes the spin configuration in `sys` to minimize energy. A total of
`maxiters` iterations will be attempted. The `method` parameter will be used in
the `optimize` function of the [Optim.jl
package](https://github.com/JuliaNLSolvers/Optim.jl). Any remaining `kwargs`
will be included in the `Options` constructor of Optim.jl.
"""
function minimize_energy!(sys::System{N}; maxiters=1000, subiters=10, method=Optim.ConjugateGradient(),
                          kwargs...) where N
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
    # which is much lower precision. Also disable check on the energy gradient
    # (g) because this has units of energy, and it is not obvious how to select
    # a natural energy scale. Disable check on `x_reltol` because `x = [zeros]`
    # is a valid configuration (all spins aligned with the stereographic
    # projection axes). This leaves only the check on `x_abstol`. Since Optim.jl
    # uses the default (p=2) norm, the acceptable x-tolerance scales like the
    # square root of the number of optimization parameters.
    x_abstol = 1e-14 * √length(αs)
    options = Optim.Options(; iterations=subiters, x_abstol, x_reltol=NaN, g_abstol=NaN, f_reltol=NaN, f_abstol=NaN, kwargs...)
    local res
    for iter in 1 : div(maxiters, subiters, RoundUp)
        res = Optim.optimize(f, g!, αs, method, options)

        # Exit if converged. Note that Hager-Zhang line search could report
        # failure if the step Δx vanishes, but for CG this is actually success.
        if Optim.converged(res) || Optim.termination_code(res) == Optim.TerminationCode.SmallXChange
            cnt = (iter-1)*subiters + res.iterations
            return cnt
        end

        # Reset stereographic projection based on current state
        ns .= normalize.(iszero(N) ? sys.dipoles : sys.coherents)
        αs .*= 0
    end

    f_abschange, x_abschange, g_residual = number_to_simple_string.((Optim.f_abschange(res), Optim.x_abschange(res), Optim.g_residual(res)); digits=2)
    @warn "Non-converged after $maxiters iterations: |ΔE|=$f_abschange, |Δx|=$x_abschange, |∂E/∂x|=$g_residual"
    return -1
end
