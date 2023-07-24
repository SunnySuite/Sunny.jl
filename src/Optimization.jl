
# Returns the stereographic projection u(α) = (2v + (1-v²)n)/(1+v²), which
# involves the orthographic projection v = (1-nn̄)α. The input `n` must be
# normalized. When `α=0`, the output is `u=n`, and when `|α|→ ∞` the output is
# `u=-n`. In all cases, `|u|=1`.
function stereographic_projection(α, n)
    @assert n'*n ≈ 1
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
@inline function vjp_stereographic_projection(x, α, n)
    v = α - n*(n'*α)
    v² = real(v'*v)
    b = (1-v²)/2
    c = 2/(1+v²)
    # Perform dot products first to avoid constructing outer-product
    x̄_dudv = c*x' - c * (x' * ((1+c*b)*n + c*v)) * v'
    # Apply projection P=1-nn̄ on right
    return x̄_dudv - (x̄_dudv * n) * n'
end


# function variance(αs)
#     ncomponents = length(αs) * length(first(αs))
#     return sum(real(α'*α) for α in αs) / ncomponents
# end

function optim_set_spins!(sys::System{0}, αs, ns)
    αs = reinterpret(reshape, Vec3, αs)
    for site in all_sites(sys)
        s = stereographic_projection(αs[site], ns[site])
        polarize_spin!(sys, s, site)
    end
end
function optim_set_spins!(sys::System{N}, αs, ns) where N
    αs = reinterpret(reshape, CVec{N}, αs)
    for site in all_sites(sys)
        Z = stereographic_projection(αs[site], ns[site])
        set_coherent_state!(sys, Z, site)
    end
end

function optim_set_gradient!(G, sys::System{0}, αs, ns)
    (αs, G) = reinterpret.(reshape, Vec3, (αs, G))
    set_energy_grad_dipoles!(G, sys.dipoles, sys)
    @. G = adjoint(vjp_stereographic_projection(G, αs, ns)) / norm(sys.dipoles)
end
function optim_set_gradient!(G, sys::System{N}, αs, ns) where N
    (αs, G) = reinterpret.(reshape, CVec{N}, (αs, G))
    set_energy_grad_coherents!(G, sys.coherents, sys)
    @. G = adjoint(vjp_stereographic_projection(G, αs, ns)) / norm(sys.coherents)
end


"""
    minimize_energy!(sys::System{N}; method=Optim.LBFGS, maxiters = 100, subiters=20, kwargs...) where N

Optimizes the spin configuration in `sys` to minimize energy. Any method from
Optim.jl that accepts only a gradient may be used by setting the `method`
keyword. Defaults to LBFGS.
"""
function minimize_energy!(sys::System{N}; method=Optim.LBFGS(), maxiters=100, subiters=20, kwargs...) where N
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
        return energy(sys) # TODO: `Sunny.energy` seems to allocate and is type-unstable (7/20/2023)
    end
    function g!(G, αs)
        optim_set_spins!(sys, αs, ns)
        optim_set_gradient!(G, sys, αs, ns)
    end

    # Perform optimization, resetting parameterization of coherent states as necessary
    options = Optim.Options(iterations=subiters, show_trace=true, kwargs...)

    for iter in 1 : div(maxiters, subiters, RoundUp)
        output = Optim.optimize(f, g!, αs, method, options)

        if Optim.converged(output)
            cnt = (iter-1)*subiters + output.iterations
            return cnt
        end

        # Reset parameterization based on current state
        ns .= normalize.(iszero(N) ? sys.dipoles : sys.coherents)
        αs .*= 0
    end

    @warn "Optimization failed to converge within $maxiters iterations."
    return -1
end
