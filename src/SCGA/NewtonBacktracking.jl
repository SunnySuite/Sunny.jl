function newton_with_backtracking(fgh!, x0; f_reltol=NaN, x_reltol=NaN, g_abstol=NaN, maxiters=20, armijo_c=1e-4, armijo_backoff=0.5, armijo_α_min=1e-4, show_trace=false)
    # Make a copy of the initial guess.
    x = copy(x0)

    # Preallocate buffers
    candidate_x = zero(x0)
    g = zero(x0)
    H = zeros(eltype(x0), length(x0), length(x0))
    p = zero(x0)

    # Evaluate objective function f, gradient, and Hessian.
    f = fgh!(0.0, g, H, x)
    maximum(abs, g) < g_abstol && return x

    function has_converged(x, candidate_x, f, candidate_f, g)
        return (!isnan(x_reltol) && isapprox(x, candidate_x; rtol=x_reltol)) ||
               (!isnan(f_reltol) && isapprox(f, candidate_f; rtol=f_reltol)) ||
               (!isnan(g_abstol) && maximum(abs, g) < g_abstol)
    end

    for k in 1:maxiters
        # Newton direction p = H \ g. If H is not guaranteed positive definite,
        # then use lu! instead of cholesky! decomposition.
        ldiv!(p, cholesky!(Hermitian(H)), g)

        # To be used for Armijo backtracking.
        g_dot_p = dot(g, p)

        # Start with full Newton step.
        α = 1.0

        # Candidate updates to x and f.
        @. candidate_x = x - α * p
        candidate_f = fgh!(0.0, g, H, candidate_x)

        # Backtracking line search until Armijo condition is satisfied:
        # f(candidate_x) ≤ f(x) - c * α * dot(g, p)
        while candidate_f > f - armijo_c * α * g_dot_p
            has_converged(x, candidate_x, f, candidate_f, g) && return candidate_x
            α < armijo_α_min && error("Failed to satisfy Armijo condition. Consider reducing armijo_α_min=$armijo_α_min or a tolerance parameter.")

            α *= armijo_backoff
            @. candidate_x = x - α * p
            candidate_f = fgh!(0.0, g, H, candidate_x)
        end

        if show_trace
            println("Iter $k: α=$α, |g|=$(norm(g)), f(x)=$candidate_f, x=$candidate_x")
        end
        has_converged(x, candidate_x, f, candidate_f, g) && return candidate_x

        # Accept candidate updates. Note that g and H will also be reused.
        x .= candidate_x
        f = candidate_f
    end

    error("Failed to converge in maxiters=$maxiters iterations.")
end
