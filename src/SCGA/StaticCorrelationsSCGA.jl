"""
    StaticCorrelationsSCGA(sys::System; measure, kT, dq)

Constructs an object to calculate [`intensities_static`](@ref) within the self
consistent gaussian approximation (SCGA). This theory assumes a classical
Boltzmann distribution with temperature `kT`. It is expected to be meangingful
above the ordering temperature, where fluctuations are approximately Gaussian.

Only `:dipole` and `:dipole_uncorrected` system modes are supported.

The theory of SCGA approximates local spin magnitude constraints with a _weaker_
global constraint condition. For each spin sublattice, the global spin sum rule
can be expressed as an integral over the unit cube ``𝐪 ∈ [0,1]^3`` for
wavevectors ``𝐪`` in reciprocal lattice units (RLU). Each such integral will be
approximated as a discrete sum over a regular grid of `floor(1/dq)^3`
wavevectors for the provided `dq` value.

If the conventional crystal cell admits a smaller primitive cell, then the SCGA
calculations can be accelerated. Construct a smaller system with
[`reshape_supercell`](@ref) and [`primitive_cell`](@ref). In this case, the
discretized ``𝐪``-point grid runs over the full Brillouin zone associated with
the primitive cell of the crystal.
"""
struct StaticCorrelationsSCGA
    sys :: System
    measure :: MeasureSpec
    β :: Float64
    λs :: Vector{Float64}

    function StaticCorrelationsSCGA(sys::System; measure::Union{Nothing, MeasureSpec}, kT::Float64, dq::Float64)
        measure = @something measure empty_measurespec(sys)
        if size(eachsite(sys)) != size(measure.observables)[2:5]
            error("Size mismatch. Check that measure is built using consistent system.")
        end

        kT > 0 || error("Temperature kT must be positive")
        β = 1 / kT

        qs = make_q_grid(sys, dq)
        Js = [fourier_exchange_matrix(sys; q) for q in qs]

        # Initial guess for the Lagrange multipliers in the physically allowed
        # space (all shifted J matrices positive definite). If `eigmin` becomes
        # a bottleneck, we can try weaker bounds using Gershgorin's circle
        # theorem: https://en.wikipedia.org/wiki/Gershgorin_circle_theorem.
        λ_init = -minimum(eigmin(J) for J in Js) + 1/β

        λs = if allequal(sys.crystal.classes)
            find_lagrange_multiplier_single(sys, Js, β, λ_init)
        else
            find_lagrange_multiplier_multi(sys, Js, β, λ_init)
        end

        return new(sys, measure, β, λs)
    end
end

function make_q_grid(sys, dq)
    wraps = [false, false, false]
    for i in 1:natoms(sys.crystal)
        for coupling in sys.interactions_union[i].pair
            (; isculled, bond) = coupling
            isculled && break
            @. wraps = wraps || !iszero(bond.n)
        end
    end

    qα = [w ? (-1/2 : dq : 1/2-dq) : [0] for w in wraps]
    return vec([to_standard_rlu(sys, Vec3(q_reshaped)) for q_reshaped in Iterators.product(qα...)])
end

# Computes the Lagrange multiplier for the standard SCGA approach with a common
# Lagrange multiplier for all sublattices.
function find_lagrange_multiplier_single(sys, Js, β, λ_init)
    evals = Iterators.flatten(eigvals(J) for J in Js)
    s² = norm2(sys.κs) * length(Js)

    function fgh!(_, gbuffer, hbuffer, λs)
        λ = λs[1]
        # λ must be large enough to shift all eigenvalues positive. Otherwise,
        # apply an infinite penalty.
        if λ + minimum(evals) <= 0
            isnothing(gbuffer) || gbuffer .= NaN
            isnothing(hbuffer) || hbuffer .= NaN
            return Inf
        end
        fbuffer = λ*s²/2 - sum(log(λ + ev) for ev in evals) / 2β
        if !isnothing(gbuffer)
            gbuffer[1] = s²/2 - sum(1 / (λ + ev) for ev in evals) / 2β
        end
        if !isnothing(hbuffer)
            hbuffer[1, 1] = sum(1 / (λ + ev)^2 for ev in evals) / 2β
        end
        return fbuffer
    end

    λs = newton_with_backtracking(fgh!, [λ_init]; x_reltol=1e-10, show_trace=false)
    return fill(λs[1], natoms(sys.crystal))
end


function find_lagrange_multiplier_multi(sys, Js, β, λ_init)
    Na = natoms(sys.crystal)
    s² = vec(sys.κs .^ 2)

    function fgh!(_, gbuffer, hbuffer, λs)
        fbuffer = 0.0
        if !isnothing(gbuffer)
            gbuffer .= 0
        end
        if !isnothing(hbuffer)
            hbuffer .= 0
        end

        Λ = Diagonal(repeat(λs, inner=3))
        A = zeros(ComplexF64, 3Na, 3Na)
        A⁻¹ = zeros(ComplexF64, 3, Na, 3, Na)

        # Determine the Lagrange multipliers λ by maximizing the "grand" free
        # energy G(λ) = log det A / 2β - ∑ᵢ λᵢ s²ᵢ / 2, where A = J + Λ.
        # Implement this numerically as minimization of the objective function f
        # = -G.
        for J in Js
            # Cholesky decomposition fails if the matrix A is not positive
            # definite. This implies unphysical λ values, which we penalize by
            # making the objective function infinite.
            @. A = J + Λ
            A_chol = cholesky!(A, RowMaximum(); check=false)
            if !issuccess(A_chol)
                isnothing(gbuffer) || gbuffer .= NaN
                isnothing(hbuffer) || hbuffer .= NaN
                return Inf
            end

            ldiv!(reshape(A⁻¹, 3Na, 3Na), A_chol, I(3Na))

            # The finite-valued objective function f
            fbuffer += (λs' * s²) / 2 - logdet(A_chol) / 2β

            # Gradient of f
            if !isnothing(gbuffer)
                for i in 1:Na
                    gbuffer[i] += s²[i] / 2 - real(tr(view(A⁻¹, :, i, :, i))) / 2β
                end
            end

            # Hessian of f
            if !isnothing(hbuffer)
                for i in 1:Na, j in 1:Na
                    hbuffer[i, j] += + norm(view(A⁻¹, :, i, :, j))^2 / 2β
                end
            end
        end

        return fbuffer
    end

    λs = fill(λ_init, Na)
    return newton_with_backtracking(fgh!, λs; x_reltol=1e-10, show_trace=false)
end


function intensities_static(scga::StaticCorrelationsSCGA, qpts)
    (; sys, measure, λs, β) = scga
    Λ = Diagonal(repeat(λs, inner=3))

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(scga.sys)

    Na = nsites(sys)
    Ncells = Na / natoms(cryst)
    Nq = length(qpts.qs)

    Nobs = num_observables(measure)
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)

    intensity = zeros(eltype(measure), Nq)

    r = sys.crystal.positions

    for (iq, q) in enumerate(qpts.qs)
        pref = zeros(ComplexF64, Nobs, Na)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        for i in 1:Na, μ in 1:3
            ff = get_swt_formfactor(measure, μ, i)
            pref[μ, i] = exp(-2π * im * dot(q_reshaped, r[i])) * compute_form_factor(ff, norm2(q_global))
        end
        Jq = fourier_exchange_matrix(sys; q)
        inverted_matrix = inv(β*Λ + β*Jq) # this is [(Iλ+J(q))^-1]^αβ_μν
        inverted_matrix = reshape(inverted_matrix, 3, Na, 3, Na)

        ssf = zero(CMat3)
        for i in 1:Na, j in 1:Na
            ssf += pref[1, i] * conj(pref[1, j]) * view(inverted_matrix, :, i, :, j)
        end

        @assert ssf ≈ ssf'
        @assert all(>=(0), real(diag(ssf)))

        map!(corrbuf, measure.corr_pairs) do (α, β)
            ssf[α, β] / Ncells
        end
        intensity[iq] = measure.combiner(q_global, corrbuf)
    end

    return StaticIntensities(sys.crystal, qpts, reshape(intensity, size(qpts.qs)))
end
