# Returns the minimum eigenvalue of the Luttinger-Tisza exchange matrix. It is a
# lower-bound on the exchange energy for magnetic order that has propagation
# wavevector k between Bravais lattice cells. This analysis ignores local
# normalization constraints for the spins within the cell.
function luttinger_tisza_exchange(sys::System; k, ϵ=0)
    @assert sys.mode in (:dipole, :dipole_large_S) "SU(N) mode not supported"
    @assert sys.latsize == (1, 1, 1) "System must have only a single cell"

    Na = natoms(sys.crystal)
    J_k = zeros(ComplexF64, 3, Na, 3, Na)

    for i in 1:natoms(sys.crystal)
        for coupling in sys.interactions_union[i].pair
            (; isculled, bond, bilin) = coupling
            isculled && break

            (; j, n) = bond
            J = exp(2π * im * dot(k, n)) * Mat3(bilin*I)
            J_k[:, i, :, j] += J / 2
            J_k[:, j, :, i] += J' / 2
        end
    end

    if !isnothing(sys.ewald)
        A = Sunny.precompute_dipole_ewald_at_wavevector(sys.crystal, (1,1,1), k) * sys.ewald.μ0_μB²
        A = reshape(A, Na, Na)
        for i in 1:Na, j in 1:Na
            J_k[:, i, :, j] += gs[i]' * A[i, j] * gs[j] / 2
        end
    end

    J_k = reshape(J_k, 3*Na, 3*Na)
    @assert diffnorm2(J_k, J_k') < 1e-15
    J_k = hermitianpart(J_k)

    E = if iszero(ϵ)
        eigmin(J_k)
    else
        # Estimate the minimum eigenvalue E as a weighted average of all small
        # eigenvalues, E = tr [exp(-β J) J] / tr exp(-β J), where β = 1/ϵ.
        # Finite but small ϵ will regularize the potential energy surface in the
        # viscinity of degenerate eigenvalues. Imposing this smoothness may aid
        # optimization of k with gradient-based methods.
        λs = eigvals(J_k)
        λmin = minimum(λs)
        # Scale all weights exp(-λ/ϵ) by the factor exp(λmin/ϵ). Scaling of
        # weights has no mathematical effect, but avoids numerical overflow.
        ws = @. exp(-(λs-λmin)/ϵ)
        sum(ws .* λs) / sum(ws)
    end

    # Scale minimum eigenvalue E by ∑ᵢSᵢ², which is valid under the L.T.
    # assumption that each spin is locally normalized.
    return E * norm2(sys.κs)
end

# Starting from an initial guess, return the wavevector k that locally minimizes
# `luttinger_tisza_exchange`.
function minimize_luttinger_tisza_exchange(sys::System; k_guess, maxiters=10_000)
    options = Optim.Options(; iterations=maxiters)

    # Work around: https://github.com/JuliaNLSolvers/LineSearches.jl/issues/175
    method = Optim.LBFGS(; linesearch=Optim.LineSearches.BackTracking(order=2))
    res = Optim.optimize(k_guess, method, options) do k
        luttinger_tisza_exchange(sys; k, ϵ=1e-8)
    end
    res = Optim.optimize(Optim.minimizer(res), Optim.ConjugateGradient(), options) do k
        luttinger_tisza_exchange(sys; k, ϵ=1e-8)
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
