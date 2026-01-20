"""
    SCGA(sys::System; measure, kT, dq)

Constructs an object to calculate [`intensities_static`](@ref) within the
self-consistent Gaussian approximation (SCGA). This theory assumes a classical
Boltzmann distribution with temperature `kT`. It is expected to be meaningful
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
struct SCGA
    sys :: System
    measure :: MeasureSpec
    β :: Float64
    extfield :: Vec3 # Applied field in energy units, as in set_field!
    λs :: Vector{Float64}
    dipoles :: Array{Vec3, 4}

    function SCGA(sys::System; measure::Union{Nothing, MeasureSpec}, kT::Real, dq::Float64)
        sys.dims == (1, 1, 1) || error("Use system with dims=(1, 1, 1) for efficiency")

        measure = @something measure empty_measurespec(sys)
        if size(eachsite(sys)) != size(measure.observables)[2:5]
            error("Size mismatch. Check that measure is built using consistent system.")
        end

        if !(sys.mode in (:dipole, :dipole_uncorrected))
            error("SCGA requires :dipole or :dipole_uncorrected mode.")
        end

        kT > 0 || error("Temperature kT must be positive")
        β = 1 / kT

        extfield = allequal(sys.extfield) ? first(sys.extfield) : error("External field must be homogeneous")

        0 < dq < 1 || error("Select q-space resolution 0 < dq < 1.")
        qs = make_q_grid(sys, dq)
        Js = [fourier_exchange_matrix(sys; q) for q in qs]

        # Initial guess for Lagrange multipliers must ensure that all shifted J
        # matrices are positive definite.
        λ_init = -minimum(eigmin(J) for J in Js) + 1/β

        (λs, dipoles) = try
            # An external field may break the symmetry-equivalence of sites
            if allequal(sys.crystal.classes) && iszero(extfield)
                (find_lagrange_multiplier_single(sys, Js, β, λ_init), zero(sys.extfield))
            else
                find_lagrange_multiplier_multi(sys, Js, β, extfield, λ_init)
            end
        catch err
            if err isa OptimizationError
                rethrow(InstabilityError("Self-consistency failed; try raising kT or refining dq"))
            else
                rethrow(err)
            end
        end

        return new(sys, measure, β, extfield, λs, dipoles)
    end
end

function magnetic_moment(scga::SCGA, site)
    site = to_cartesian(site)
    return - scga.sys.gs[site] * scga.dipoles[site]
end

function magnetic_moment_per_site(scga::SCGA)
    Statistics.mean(magnetic_moment(scga, site) for site in eachsite(scga.sys))
end

function make_q_grid(sys, dq)
    # Round up to integer grid length
    dq = 1 / round(1 / dq, RoundUp)
    wraps = [false, false, false]
    for int in sys.interactions_union
        for coupling in int.pair
            @. wraps = wraps || !iszero(coupling.bond.n)
        end
    end

    qα = [w ? (0 : dq : 1-dq) : [0] for w in wraps]
    return vec([to_standard_rlu(sys, Vec3(q_reshaped)) for q_reshaped in Iterators.product(qα...)])
end

# If all sites are symmetry-equivalent, then solve for a single Lagrange
# multiplier λ. This has energy units and effectively shifts J(q) → J(q) + λ.
# Traditional SCGA notation, e.g. Conlon and Chalker, would use instead the
# dimensionless Lagrange multiplier λ_C = β (λ - eigmin J(q)).
function find_lagrange_multiplier_single(sys, Js, β, λ_init)
    # Calculate J(q) eigenvalues once now and then shift them by λ later
    evals = reduce(vcat, eigvals.(Js))
    Nq = length(Js)
    s² = vec(sys.κs .^ 2)
    sum_s² = sum(s²)

    function fgh!(_, gbuffer, hbuffer, λs)
        λ = λs[1]
        # λ must be large enough to shift all eigenvalues positive. Otherwise,
        # apply an infinite penalty.
        if λ + minimum(evals) <= 0
            isnothing(gbuffer) || (gbuffer .= NaN)
            isnothing(hbuffer) || (hbuffer .= NaN)
            return Inf
        end
        fbuffer = λ*sum_s²/2 - sum(log(λ + ev) for ev in evals) / (2β*Nq)
        if !isnothing(gbuffer)
            gbuffer[1] = sum_s²/2 - sum(1 / (λ + ev) for ev in evals) / (2β*Nq)
        end
        if !isnothing(hbuffer)
            hbuffer[1, 1] = sum(1 / (λ + ev)^2 for ev in evals) / (2β*Nq)
        end
        return fbuffer
    end

    g_abstol = 1e-8 * Statistics.mean(s²)
    armijo_slack = 1e-8 * sum_s² / β
    λs = newton_with_backtracking(fgh!, [λ_init]; g_abstol, armijo_slack)
    return fill(λs[1], natoms(sys.crystal))
end


function find_lagrange_multiplier_multi(sys, Js, β, extfield, λ_init)
    Na = natoms(sys.crystal)
    Nq = length(Js)
    s² = vec(sys.κs .^ 2)

    # Handle Zeeman coupling when processing J₁ = J(q=0). Elements must be
    # exactly real.
    @assert first(Js) ≈ fourier_exchange_matrix(sys; q=zero(Vec3))
    @assert iszero(imag(first(Js)))

    # Zeeman coupling enters as: - ∑ᵢ Sᵢ⋅bᵢ
    b = [- g' * extfield for g in sys.gs]
    # Expected spin dipoles ⟨Sᵢ⟩
    S = zero(b)

    A = zeros(ComplexF64, 3Na, 3Na)
    A⁻¹ = zeros(ComplexF64, 3, Na, 3, Na)

    function fgh!(_, gbuffer, hbuffer, λs)
        fbuffer = 0.0
        isnothing(gbuffer) || (gbuffer .= 0)
        isnothing(hbuffer) || (hbuffer .= 0)

        Λ = Diagonal(repeat(λs, inner=3))

        # Determine the Lagrange multipliers λ by maximizing (not minimizing!)
        # the "grand" free energy G(λ) = log det A / 2β - ∑ᵢ λᵢ s²ᵢ / 2, where A
        # = J + Λ. Implement this numerically as minimization of the objective
        # function f = -G.
        for (iq, J) in enumerate(Js)
            # Cholesky decomposition fails if the matrix A is not positive
            # definite. This implies unphysical λ values, which we penalize by
            # making the objective function infinite.
            @. A = J + Λ
            A_chol = cholesky!(A, RowMaximum(); check=false)
            if !issuccess(A_chol)
                isnothing(gbuffer) || (gbuffer .= NaN)
                isnothing(hbuffer) || (hbuffer .= NaN)
                return Inf
            end

            ldiv!(reshape(A⁻¹, 3Na, 3Na), A_chol, I(3Na))

            # The finite-valued objective function f
            fbuffer += ((λs'*s²)/2 - logdet(A_chol)/2β) / Nq

            # Gradient of f
            if !isnothing(gbuffer)
                for i in 1:Na
                    gbuffer[i] += (s²[i]/2 - real(tr(view(A⁻¹, :, i, :, i)))/2β) / Nq
                end
            end

            # Hessian of f
            if !isnothing(hbuffer)
                for i in 1:Na, j in 1:Na
                    hbuffer[i, j] += + norm2(view(A⁻¹, :, i, :, j)) / (2β*Nq)
                end
            end

            # In case of J(q=0), account for Zeeman coupling
            if iq == 1 && !iszero(extfield)
                # Solve S = A \ b
                S_vec, b_vec = vec.(reinterpret.(Float64, (S, b)))
                ldiv!(S_vec, A_chol, b_vec)

                fbuffer += real(dot(b, S))/2

                if !isnothing(gbuffer)
                    for i in 1:Na
                        gbuffer[i] += -norm2(S[i])/2
                    end
                end

                if !isnothing(hbuffer)
                    for i in 1:Na, j in 1:Na
                        hbuffer[i, j] += real(dot(S[i], Mat3(view(A⁻¹, :, i, :, j)), S[j]))
                    end
                end
            end
        end

        return fbuffer
    end

    λs = fill(λ_init, Na)
    g_abstol = 1e-8 * Statistics.mean(s²)
    armijo_slack = 1e-8 * sum(s²) / β
    λs = newton_with_backtracking(fgh!, λs; g_abstol, armijo_slack)
    # A final call to ensure that S is updated for the latest λs
    iszero(extfield) || fgh!(0.0, nothing, nothing, λs)
    return (λs, S)
end


# Returns matrix of sensitivities ∂λᵢ/∂θₖ for labeled parameters θₖ. The λ are
# defined to satisfy g(θ, λ(θ)) = 0, with g = ∂f/∂λ the gradient of the
# objective f. By the implicit function theorem, the gradient is ∂λ/∂θₖ =
# - H⁻¹ vₖ, involving the Hessian H = ∂g/∂λ and the vectors of mixed partials,
# vₖ = ∂g/∂θₖ.
function lagrange_multiplier_jacobian(sys, qs, β, λs, labels)
    Na = natoms(sys.crystal)

    J = zeros(ComplexF64, 3Na, 3Na)
    ∂J = zeros(ComplexF64, 3Na, 3Na)

    Λ = Diagonal(repeat(λs, inner=3))
    A = zeros(ComplexF64, 3Na, 3Na)
    A⁻¹ = zeros(ComplexF64, 3Na, 3Na)
    A⁻¹_block = reshape(A⁻¹, 3, Na, 3, Na)
    H = zeros(Float64, Na, Na)

    v = zeros(Float64, Na, length(labels))

    for q in qs
        fourier_exchange_matrix!(J, sys; q)
        @. A = J + Λ
        A_chol = cholesky!(A, RowMaximum())
        ldiv!(A⁻¹, A_chol, I(3Na))
        for i in 1:Na, j in 1:Na
            H[i, j] += + norm(view(A⁻¹_block, :, i, :, j))^2 / 2β
        end

        for (k, label) in enumerate(labels)
            fourier_exchange_matrix_sensitivity!(∂J, sys, label; q)
            for i in 1:Na
                y_i = reshape(view(A⁻¹_block, :, :, :, i), 3Na, 3)
                v[i, k] += real(dot(y_i, ∂J, y_i)) / 2β # TODO: remove allocation?
            end
        end
    end

    return - H \ v
end


function intensities_static(scga::SCGA, qpts; measure=nothing)
    (; sys, λs, β) = scga
    measure = @something measure scga.measure
    Λ = Diagonal(repeat(λs, inner=3))

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(scga.sys)

    Na = nsites(sys)
    Ncells = Na / natoms(cryst)

    # Temporary storage for pair correlations
    Nobs = num_observables(measure)
    Ncorr = num_correlations(measure)
    corr = zeros(ComplexF64, Ncorr)

    # Preallocation
    A = zeros(ComplexF64, 3Na, 3Na)
    O = view(measure.observables::Array{Vec3, 5}, :, 1, 1, 1, :)
    X = zeros(ComplexF64, 3Na, Nobs)
    pref = zeros(ComplexF64, 3Na, Nobs)
    pref_reshaped = reshape(pref, 3, Na, Nobs)

    data = map(qpts.qs) do q
        any(isnan, q) && return NaN * zero(eltype(measure))

        q_global = cryst.recipvecs * q

        for i in 1:Na, μ in 1:Nobs
            r_global = global_position(sys, (1, 1, 1, i)) # + offsets[μ, i]
            ff = get_swt_formfactor(measure, μ, i)
            c = exp(+ im * dot(q_global, r_global)) * compute_form_factor(ff, norm2(q_global))
            for α in 1:3
                pref_reshaped[α, i, μ] = c * O[μ, i][α]
            end
        end

        fourier_exchange_matrix!(A, sys; q)
        A .+= Λ
        A .*= β

        A_chol = cholesky!(A; check=false)
        issuccess(A_chol) || InstabilityError("Self-consistency failed at q = $(vec3_to_string(q)); try raising kT or refining dq")
        ldiv!(X, A_chol, pref)
        map!(corr, measure.corr_pairs) do (μ, ν)
            return dot(view(pref, :, μ), view(X, :, ν)) / Ncells
        end

        #=
        A⁻¹ = reshape(inv(β*Λ + β*Jq), 3, Na, 3, Na)
        map!(corr, measure.corr_pairs) do (μ, ν)
            acc = zero(ComplexF64)
            for α in 1:3, β in 1:3, i in 1:Na, j in 1:Na
                acc += conj(pref[α, i, μ]) * A⁻¹[α, i, β, j] * pref[β, j, ν]
            end
            return acc / Ncells
        end
        =#

        return measure.combiner(q_global, corr)
    end

    return StaticIntensities(cryst, qpts, data)
end


"""
    magnetic_susceptibility_per_site(scga::SCGA)

Returns the magnetic susceptibility tensor in units of inverse energy,
``\\tilde{χ} = (dμ/d𝐁) / μ_B^2`` , where ``μ`` is the magnetic dipole per site,
``𝐁`` is the physical applied field, and ``μ_B`` is the Bohr magneton. In terms
of Sunny quantities, ``\\tilde{χ}`` can be understood as the derivative of
[`magnetic_moment_per_site`](@ref) (dimensionless) with respect to the argument
to [`set_field!`](@ref) (energy units).

At zero field, fluctuation-dissipation states ``⟨M^α(q=0) M^β(q=0)⟩ / μ_B^2 N_s
= k_B T \\tilde{χ}^{αβ}``, where ``M(q=0)`` is the magnetic dipole summed over
all ``N_s`` sites in the sample. The structure factor on the left-hand side can
be calculated as [`intensities_static`](@ref) (a per-cell quantity) divided by
the number of sites in the chemical unit cell.

!!! "Conversion to physical molar susceptibility units"

    The molar susceptibility in a Gaussian unit system is ``χ = (N_A μ_0 μ_B^2 /
    4π×10^{-6}) \\tilde{χ}``. For specific [`Units`](@ref), the conversion factor to
    (emu/Oe/mol) is provided by `units.cgs_molar_susceptibility` (inverse energy).
    In inverse meV units, for example, `units.cgs_molar_susceptibility` represents

    ```math
    \\frac{\\mathrm{emu/Oe/mol}}{(N_A μ_0 μ_B^2 / 4π×10^{-6})} = 30.9331… / \\mathrm{meV}.
    ```

    Alternatively, the conversion factor to (m³/mol) in SI units is provided by
    `units.si_molar_susceptibility` (inverse energy), which uses the identity

    ```math
    1 \\mathrm{emu/Oe/mol} = 4π×10^{-6} \\mathrm{m³/mol}.
    ```

# Example

```julia
units = Units(:meV, :angstrom)

# Magnetic susceptibility per site in inverse energy (1/meV)
χ̃ = magnetic_susceptibility_per_site(scga)

# Molar susceptibility in Gaussian units (emu/Oe/mol)
χ = χ̃ / units.cgs_molar_susceptibility
```
"""
function magnetic_susceptibility_per_site(scga)
    iszero(scga.extfield) || error("Bulk susceptibility currently requires zero field")
    measure = ssf_custom((q, ssf) -> ssf, scga.sys)
    # Fluctuation dissipation: dμ/dB = ⟨δμ, δμ⟩/kT in inverse energy units
    cryst = orig_crystal(scga.sys)
    χ = intensities_static(scga, [[0, 0, 0]]; measure).data[1] * scga.β / natoms(cryst)
    @assert iszero(imag(χ))
    return real(χ)
end

### Autodiff support

function CRC.rrule(::typeof(SCGA), sys::System; measure, kT, dq)
    (; active_labels) = sys
    scga = SCGA(sys; measure, kT, dq)
    proj_scga = CRC.ProjectTo(scga)

    function pullback(Δscga)
        Δscga = proj_scga(CRC.unthunk(Δscga))
        if Δscga isa CRC.NoTangent || Δscga isa CRC.AbstractZero
            return (CRC.NoTangent(), CRC.ZeroTangent())
        else
            @assert !isempty(active_labels)
        end

        # Add the implicit λ pathway: Δθ += (∂λ/∂θ)' * Δλ
        qs = make_q_grid(sys, dq)
        Jλ = lagrange_multiplier_jacobian(sys, qs, scga.β, scga.λs, active_labels)  # Na × nlabels
        Δvals = Δscga.sys.vals + Jλ' * Δscga.λs

        return (CRC.NoTangent(), SystemTangent(Δvals))
    end

    return scga, pullback
end

# Evaluate pullback for d = combiner(q_global, corr) given cotangent Δd,
# assuming linearity in corr.
function combiner_ad(rc::CRC.RuleConfig, combiner, q_global, corr, Δd)
    Δcorr = similar(corr)
    e = zeros(Float64, length(corr))
    for i in 1:length(corr)
        fill!(e, 0.0)
        e[i] = 1.0
        Δcorr[i] = dot(combiner(q_global, e), Δd)
    end

    # Test correctness of adjoint. Use Float64 for elements of r because some
    # combiners (like in ssf_trace) require a real input.
    r = randn(Float64, length(corr))
    matches = dot(Δd, combiner(q_global, r)) ≈ dot(Δcorr, r)

    # Fall back to much slower AD if needed
    if !matches
        @warn "Combiner appears nonlinear; falling back to generic AD"
        _, pullback = CRC.rrule_via_ad(rc, combiner, q_global, corr)
        _, _, Δcorr = pullback(Δd)
    end

    return Δcorr
end

function CRC.rrule(rc::CRC.RuleConfig, ::typeof(intensities_static), scga::SCGA, qpts; measure=nothing)
    qpts = convert(AbstractQPoints, qpts)
    res = intensities_static(scga, qpts; measure)
    proj_res = CRC.ProjectTo(res)

    function pullback(Δres)
        Δres = proj_res(CRC.unthunk(Δres))
        if Δres isa CRC.NoTangent || Δres isa CRC.AbstractZero
            return CRC.NoTangent(), CRC.ZeroTangent(), CRC.NoTangent()
        end

        Δdata = Δres.data

        (; sys, λs, β) = scga
        measure = @something measure scga.measure
        cryst = orig_crystal(sys)
        Na = nsites(sys)
        Ncells = Na / natoms(cryst)

        Nobs  = num_observables(measure)
        Ncorr = num_correlations(measure)

        # Forward-work buffers
        Λ = Diagonal(repeat(λs, inner=3))
        A = zeros(ComplexF64, 3Na, 3Na)
        X = zeros(ComplexF64, 3Na, Nobs)
        pref = zeros(ComplexF64, 3Na, Nobs)
        pref_reshaped = reshape(pref, 3, Na, Nobs)
        corr = zeros(ComplexF64, Ncorr)
        O = view(measure.observables::Array{Vec3,5}, :, 1, 1, 1, :)

        # Reverse buffers
        Δvals = zeros(Float64, length(sys.active_labels))
        Δλs = zeros(Float64, length(λs))

        ΔX = zeros(ComplexF64, 3Na, Nobs)
        G = zeros(ComplexF64, 3Na, Nobs)
        ΔA = zeros(ComplexF64, 3Na, 3Na)
        ΔA_reshaped = reshape(ΔA, 3, Na, 3, Na)
        ∂J = zeros(ComplexF64, 3Na, 3Na)

        foreach(qpts.qs, res.data, Δdata) do q, d, Δd
            any(isnan, q) && return

            ### REPEAT FORWARD CALCULATION (no tape to avoid memory costs)

            q_global = cryst.recipvecs * q

            for i in 1:Na, μ in 1:Nobs
                r_global = global_position(sys, (1, 1, 1, i)) # + offsets[μ, i]
                ff = get_swt_formfactor(measure, μ, i)
                c = exp(+ im * dot(q_global, r_global)) * compute_form_factor(ff, norm2(q_global))
                for α in 1:3
                    pref_reshaped[α, i, μ] = c * O[μ, i][α]
                end
            end

            # A = β*(J+Λ)
            fourier_exchange_matrix!(A, sys; q)
            A .+= Λ
            A .*= β

            # X = A \ pref
            A_chol = cholesky!(A; check=false)
            issuccess(A_chol) || InstabilityError("Self-consistency failed at q = $(vec3_to_string(q)); try raising kT or refining dq")
            ldiv!(X, A_chol, pref)

            # corr = dot(pref_μ, X_ν) / Ncells
            map!(corr, measure.corr_pairs) do (μ, ν)
                return dot(view(pref, :, μ), view(X, :, ν)) / Ncells
            end

            # data = measure.combiner(corr)
            @assert d ≈ measure.combiner(q_global, corr)

            ### BACKWARD CALCULATION

            # Pullback on: data = measure.combiner(corr)
            Δcorr = combiner_ad(rc, measure.combiner, q_global, corr, Δd)
            Δcorr = CRC.unthunk(Δcorr)
            if Δcorr isa CRC.NoTangent || Δcorr isa CRC.AbstractZero || isnothing(Δcorr)
                return
            end

            # Pullback on: corr = dot(pref_μ, X_ν) / Ncells
            fill!(ΔX, 0)
            for (k, (μ, ν)) in enumerate(measure.corr_pairs)
                view(ΔX, :, ν) .+= Δcorr[k] .* view(pref, :, μ) ./ Ncells
            end

            # Pullback on: X = A \ pref
            ldiv!(G, A_chol', ΔX) # G = A' \ ΔX
            mul!(ΔA, G, X')       # ΔA = G * X'
            ΔA .*= -1             # ΔA = - (A' \ ΔX) * X'

            # Pullback on: A = β*(J+Λ) for J(vals)
            for (k, label) in enumerate(sys.active_labels)
                fourier_exchange_matrix_sensitivity!(∂J, sys, label; q)
                Δvals[k] += β * real(dot(vec(ΔA), vec(∂J)))
            end

            # Pullback on: A = β*(J+Λ) for Λ = Diagonal(repeat(λs, inner=3))
            for i in 1:length(λs)
                Δλs[i] += β * real(ΔA_reshaped[1, i, 1, i] +
                                   ΔA_reshaped[2, i, 2, i] +
                                   ΔA_reshaped[3, i, 3, i])
            end
        end

        Δscga = CRC.Tangent{typeof(scga)}(
            sys = SystemTangent(Δvals),
            λs  = Δλs,
        )

        return CRC.NoTangent(), Δscga, CRC.NoTangent()
    end

    return res, pullback
end
