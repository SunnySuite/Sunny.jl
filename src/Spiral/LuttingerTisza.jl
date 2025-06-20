function fourier_exchange_matrix!(Jq::Matrix{ComplexF64}, sys::System; q)
    @assert sys.mode in (:dipole, :dipole_uncorrected) "SU(N) mode not supported"
    @assert sys.dims == (1, 1, 1) "System must have only a single cell"
    q_reshaped = to_reshaped_rlu(sys, q)

    Na = natoms(sys.crystal)
    Jq .= 0
    Jq = reshape(Jq, 3, Na, 3, Na)

    for (i, int) in enumerate(sys.interactions_union)
        for coupling in int.pair
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break
            iszero(biquad) || error("Biquadratic interactions not supported")

            @assert i == bond.i
            j = bond.j

            J = exp(2π * im * dot(q_reshaped, bond.n)) * Mat3(bilin*I)
            view(Jq, :, i, :, j) .+= J
            view(Jq, :, j, :, i) .+= J'
        end
    end

    if !isnothing(sys.ewald)
        (; demag, μ0_μB²) = sys.ewald
        Aq = precompute_dipole_ewald_at_wavevector(sys.crystal, (1,1,1), demag, q_reshaped) * μ0_μB²
        Aq = reshape(Aq, Na, Na)
        for i in 1:Na, j in 1:Na
            view(Jq, :, i, :, j) .+= sys.gs[i]' * Aq[i, j] * sys.gs[j]
        end
    end

    for i in 1:Na
        onsite_coupling = sys.interactions_union[i].onsite
        (; c2, c4, c6) = onsite_coupling
        iszero(c4) && iszero(c6) || error("Single-ion anisotropy beyond quadratic order not supported")
        anisotropy = SA[c2[1]-c2[3]        c2[5] 0.5c2[2];
                              c2[5] -c2[1]-c2[3] 0.5c2[4];
                           0.5c2[2]     0.5c2[4]   2c2[3]]
        view(Jq, :, i, :, i) .+= 2 * anisotropy
    end

    Jq = reshape(Jq, 3Na, 3Na)
    @assert diffnorm2(Jq, Jq') < 1e-15
    return hermitianpart!(Jq)
end

function fourier_exchange_matrix(sys::System; q)
    Na = natoms(sys.crystal)
    Jq = zeros(ComplexF64, 3Na, 3Na)
    return fourier_exchange_matrix!(Jq, sys; q)
end

# Returns the Luttinger-Tisza predicted exchange energy associated with the
# propagation wavevector k between chemical cells. The LT analysis minimizes
# energy E = (1/2) Sₖ† Jₖ Sₖ, where Sₖ is some length-3Nₐ vector and Jₖ is the
# 3Nₐ×3Nₐ exchange matrix. Given a minimum energy eigenpair (ϵₖ, Sₖ) the
# LT-predicted energy is |Sₖ|² ϵₖ / 2 with normalization |Sₖ|² = ∑ᵢ|Sᵢ|², where
# index i denotes spin sublattice of the chemical cell. If the components of Sₖ
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
function minimize_luttinger_tisza_exchange(sys::System; k_guess, maxiters=10_000)
    options = Optim.Options(; iterations=maxiters)

    # Work around: https://github.com/JuliaNLSolvers/LineSearches.jl/issues/175
    method = Optim.LBFGS(; linesearch=Optim.LineSearches.BackTracking(order=2))
    res = Optim.optimize(k_guess, method, options) do k
        luttinger_tisza_exchange(sys; k, η=1e-8)
    end
    res = Optim.optimize(Optim.minimizer(res), Optim.ConjugateGradient(), options) do k
        luttinger_tisza_exchange(sys; k, η=1e-8)
    end

    if Optim.converged(res)
        k = Optim.minimizer(res)
        # Wrap components to [0, 1)
        k = wrap_to_unit_cell(Vec3(k); symprec=1e-6)
        return k
    else
        error("Momentum optimization failed to converge within $maxiters iterations.")
    end
end
