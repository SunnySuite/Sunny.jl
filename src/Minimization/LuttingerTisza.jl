# Returns the Luttinger-Tisza predicted exchange energy associated with the
# propagation wavevector k between chemical cells. The LT analysis minimizes
# energy E = (1/2) Sₖ† Jₖ Sₖ, where Sₖ is some length-3Nₐ vector and Jₖ is the
# 3Nₐ×3Nₐ exchange matrix. Given a minimum energy eigenpair (ϵₖ, Sₖ) the
# LT-predicted energy is |Sₖ|² ϵₖ / 2 with normalization |Sₖ|² = ∑ᵢ|Sᵢ|², where
# index i denotes a site of the magnetic cell. If the components of Sₖ
# satisfy local spin normalization constraints, then the LT energy minimized
# over k is physically correct for the spiral ground state. In practice, the
# eigenvector Sₖ may violate local spin normalization and the LT-predicted
# energy is only a lower-bound on the exchange energy. Conditions for
# correctness are given by Xiong and Wen (2013), arXiv:1208.1512.
function luttinger_tisza_exchange(sys::System; k, η=0)
    J_k = fourier_exchange_matrix(sys; q=k)

    E = if iszero(η)
        eigmin(J_k) / 2
    else
        # Estimate the minimum eigenvalue E as a weighted average of all small
        # eigenvalues, E = tr [exp(-β J) J] / tr exp(-β J), where β = 1/η.
        # Finite but small η will regularize the potential energy surface in the
        # viscinity of degenerate eigenvalues. Imposing this smoothness may aid
        # optimization of k with gradient-based methods.
        ϵs = eigvals(J_k) / 2
        ϵmin = minimum(ϵs)
        # Scale all weights exp(-λ/η) by the factor exp(λmin/η). Scaling of
        # weights has no mathematical effect, but avoids numerical overflow.
        ws = @. exp(-(ϵs-ϵmin)/η)
        sum(ws .* ϵs) / sum(ws)
    end

    # Rescale by ∑ᵢSᵢ², which would be valid if the minimum-energy eigenvector
    # respects the local spin normalization constraints.
    return E * norm2(sys.κs)
end

# Starting from an initial guess, return the wavevector k that locally minimizes
# `luttinger_tisza_exchange`.
function minimize_luttinger_tisza_exchange(sys::System; k_guess, maxiters=1000)
    # Default "low precision" settings seem OK here. If higher precision is
    # required, see considerations documented in `minimize_energy!`.
    options = Optim.Options(; iterations=maxiters)
    # For simplicity, allow Optim.jl to estimate the gradient components using
    # finite differences. These estimates are not very reliable and cause
    # occasional test failures in `test_spiral.jl`.
    res = Optim.optimize(k_guess, Optim.ConjugateGradient(), options) do k
        luttinger_tisza_exchange(sys; k, η=1e-8)
    end

    if Optim.converged(res)
        k = Optim.minimizer(res)
        # Wrap components to [0, 1)
        k = wrap_to_unit_cell(Vec3(k); tol=1e-6)
        return k
    else
        error("Momentum optimization failed to converge within $maxiters iterations.")
    end
end
