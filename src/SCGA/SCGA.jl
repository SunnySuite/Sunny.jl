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
    λs :: Vector{Float64}

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

        0 < dq < 1 || error("Select q-space resolution 0 < dq < 1.")
        qs = make_q_grid(sys, dq)
        Js = [fourier_exchange_matrix(sys; q) for q in qs]

        # Initial guess for Lagrange multipliers must ensure that all shifted J
        # matrices are positive definite.
        λ_init = -minimum(eigmin(J) for J in Js) + 1/β

        λs = try
            if allequal(sys.crystal.classes)
                find_lagrange_multiplier_single(sys, Js, β, λ_init)
            else
                find_lagrange_multiplier_multi(sys, Js, β, λ_init)
            end
        catch err
            rethrow(InstabilityError("Self-consistency failed; try raising kT or refining dq"))
        end

        return new(sys, measure, β, λs)
    end
end

function make_q_grid(sys, dq)
    wraps = [false, false, false]
    for int in sys.interactions_union
        for coupling in int.pair
            @. wraps = wraps || !iszero(coupling.bond.n)
        end
    end

    qα = [w ? (-1/2 : dq : 1/2-dq) : [0] for w in wraps]
    return vec([to_standard_rlu(sys, Vec3(q_reshaped)) for q_reshaped in Iterators.product(qα...)])
end

# Computes the Lagrange multiplier for the standard SCGA approach with a common
# Lagrange multiplier for all sublattices.
function find_lagrange_multiplier_single(sys, Js, β, λ_init)
    evals = reduce(vcat, eigvals.(Js))
    Nq = length(Js)
    s² = vec(sys.κs .^ 2)
    sum_s² = sum(s²)

    function fgh!(_, gbuffer, hbuffer, λs)
        λ = λs[1]
        # λ must be large enough to shift all eigenvalues positive. Otherwise,
        # apply an infinite penalty.
        if λ + minimum(evals) <= 0
            isnothing(gbuffer) || gbuffer .= NaN
            isnothing(hbuffer) || hbuffer .= NaN
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


function find_lagrange_multiplier_multi(sys, Js, β, λ_init)
    Na = natoms(sys.crystal)
    Nq = length(Js)
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

        # Determine the Lagrange multipliers λ by maximizing (not minimizing!)
        # the "grand" free energy G(λ) = log det A / 2β - ∑ᵢ λᵢ s²ᵢ / 2, where A
        # = J + Λ. Implement this numerically as minimization of the objective
        # function f = -G.
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
                    hbuffer[i, j] += + norm(view(A⁻¹, :, i, :, j))^2 / (2β*Nq)
                end
            end
        end

        return fbuffer
    end

    λs = fill(λ_init, Na)
    g_abstol = 1e-8 * Statistics.mean(s²)
    armijo_slack = 1e-8 * sum(s²) / β
    return newton_with_backtracking(fgh!, λs; g_abstol, armijo_slack)
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


function intensities_static(scga::SCGA, qpts)
    (; sys, measure, λs, β) = scga
    Λ = Diagonal(repeat(λs, inner=3))

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(scga.sys)

    Na = nsites(sys)
    Ncells = Na / natoms(cryst)
    Nq = length(qpts.qs)

    # Temporary storage for pair correlations
    Nobs = num_observables(measure)
    Ncorr = num_correlations(measure)
    corrbuf = zeros(ComplexF64, Ncorr)

    # Preallocation
    A = zeros(ComplexF64, 3Na, 3Na)
    O = view(measure.observables::Array{Vec3, 5}, :, 1, 1, 1, :)
    X = zeros(ComplexF64, 3Na, Nobs)
    pref = zeros(ComplexF64, 3Na, Nobs)
    pref_reshaped = reshape(pref, 3, Na, Nobs)
    intensity = zeros(eltype(measure), Nq)

    for (iq, q) in enumerate(qpts.qs)
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
        map!(corrbuf, measure.corr_pairs) do (μ, ν)
            return dot(view(pref, :, μ), view(X, :, ν)) / Ncells
        end

        #=
        A⁻¹ = reshape(inv(β*Λ + β*Jq), 3, Na, 3, Na)
        map!(corrbuf, measure.corr_pairs) do (μ, ν)
            acc = zero(ComplexF64)
            for α in 1:3, β in 1:3, i in 1:Na, j in 1:Na
                acc += conj(pref[α, i, μ]) * A⁻¹[α, i, β, j] * pref[β, j, ν]
            end
            return acc / Ncells
        end
        =#

        intensity[iq] = measure.combiner(q_global, corrbuf)
    end

    return StaticIntensities(cryst, qpts, reshape(intensity, size(qpts.qs)))
end



########################## MOVE LATER

CRC.@non_differentiable q_space_path(::Any, ::Any, ::Any)

function with_params(sys::System, labels::Vector{Symbol}, vals::Vector{<: Real})
    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    return sys
end

# Custom System tangent that stores param sensitivities

struct SystemTangent <: CRC.AbstractTangent
    vals::Vector{Float64}
end

Base.:+(a::SystemTangent, b::SystemTangent) = SystemTangent(a.vals .+ b.vals)
CRC.zero_tangent(t::SystemTangent) = SystemTangent(zero(t.vals))
CRC.unthunk(t::SystemTangent) = t

function CRC.ProjectTo(sys::System)
    n = length(sys.active_labels)
    function project(Δsys)
        Δsys = CRC.unthunk(Δsys)
        if Δsys isa CRC.NoTangent || Δsys isa CRC.AbstractZero
            return SystemTangent(zeros(n))
        elseif Δsys isa SystemTangent
            return Δsys
        else
            # if hasproperty(Δsys, :vals)
            #     return SystemTangent(getproperty(Δsys, :vals))
            # end
            error("Unsupported cotangent for System: $(typeof(Δsys))")
        end
    end
    return project
end


function CRC.rrule(::typeof(with_params), sys::System, labels, vals)
    sys2 = with_params(sys, labels, vals)
    proj_sys = CRC.ProjectTo(sys2)
    proj_vals = CRC.ProjectTo(vals)

    function pullback(Δsys2)
        Δsys2 = proj_sys(CRC.unthunk(Δsys2))
        return (CRC.NoTangent(), CRC.NoTangent(), CRC.NoTangent(), proj_vals(Δsys2.vals))
    end

    return sys2, pullback
end


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


# Evaluate pullback for data = combiner(q_global, corrbuff) given cotangent
# Δdata, assuming linearity in corrbuff.
function combiner_ad(rc::CRC.RuleConfig, combiner, q_global, corrbuf, Δdata)
    Δcorrbuf = similar(corrbuf)
    e = zeros(Float64, length(corrbuf))
    for i in 1:length(corrbuf)
        fill!(e, 0.0)
        e[i] = 1.0
        Δcorrbuf[i] = dot(combiner(q_global, e), Δdata)
    end

    # Test correctness of adjoint. Use Float64 for elements of r because some
    # combiners (like in ssf_trace) require a real input.
    r = randn(Float64, length(corrbuf))
    matches = dot(Δdata, combiner(q_global, r)) ≈ dot(Δcorrbuf, r)

    # Fall back to much slower AD if needed
    if !matches
        @warn "Combiner appears nonlinear; falling back to generic AD"
        _, pullback = CRC.rrule_via_ad(rc, combiner, q_global, corrbuf)
        _, _, Δcorrbuf = pullback(Δdata)
    end

    return Δcorrbuf
end

function CRC.rrule(rc::CRC.RuleConfig, ::typeof(intensities_static), scga::SCGA, qpts)
    qpts = convert(AbstractQPoints, qpts)
    res = intensities_static(scga, qpts)
    proj_res = CRC.ProjectTo(res)

    function pullback(Δres)
        Δres = proj_res(CRC.unthunk(Δres))
        if Δres isa CRC.NoTangent || Δres isa CRC.AbstractZero
            return CRC.NoTangent(), CRC.ZeroTangent(), CRC.NoTangent()
        end

        Δdata = Δres.data

        (; sys, measure, λs, β) = scga
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
        corrbuf = zeros(ComplexF64, Ncorr)
        O = view(measure.observables::Array{Vec3,5}, :, 1, 1, 1, :)

        # Reverse buffers
        Δvals = zeros(Float64, length(sys.active_labels))
        Δλs = zeros(Float64, length(λs))

        ΔX = zeros(ComplexF64, 3Na, Nobs)
        G = zeros(ComplexF64, 3Na, Nobs)
        ΔA = zeros(ComplexF64, 3Na, 3Na)
        ΔA_reshaped = reshape(ΔA, 3, Na, 3, Na)
        ∂J = zeros(ComplexF64, 3Na, 3Na)

        for (iq, q) in enumerate(qpts.qs)

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

            # corrbuf = dot(pref_μ, X_ν) / Ncells
            map!(corrbuf, measure.corr_pairs) do (μ, ν)
                return dot(view(pref, :, μ), view(X, :, ν)) / Ncells
            end

            # data = measure.combiner(corrbuf)
            @assert res.data[iq] ≈ measure.combiner(q_global, corrbuf)

            ### BACKWARD CALCULATION

            # Pullback on: data = measure.combiner(corrbuf)
            Δcorrbuf = combiner_ad(rc, measure.combiner, q_global, corrbuf, Δdata[iq])
            Δcorrbuf = CRC.unthunk(Δcorrbuf)
            if Δcorrbuf isa CRC.NoTangent || Δcorrbuf isa CRC.AbstractZero || isnothing(Δcorrbuf)
                continue
            end

            # Pullback on: corrbuf = dot(pref_μ, X_ν) / Ncells
            fill!(ΔX, 0)
            for (k, (μ, ν)) in enumerate(measure.corr_pairs)
                view(ΔX, :, ν) .+= Δcorrbuf[k] .* view(pref, :, μ) ./ Ncells
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


###########


struct FittingLoss{F}
    f :: F
    sys :: System
    labels :: Vector{Symbol}
end

function fitting_loss(f, sys, labels)
    return FittingLoss(f, sys, labels)
end

function (fl::FittingLoss)(vals)
    (; f, sys, labels) = fl

    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    try
        return f(sys)
    catch err
        (err isa InstabilityError) ? Inf : rethrow(err)
    end
end

function CRC.rrule(rc::CRC.RuleConfig, fl::FittingLoss, vals)
    (; f, sys, labels) = fl

    sys = clone_system(sys)
    set_params!(sys, labels, vals)
    sys.active_labels = labels
    (y, f_pb) = try
        CRC.rrule_via_ad(rc, f, sys)
    catch err
        (err isa InstabilityError) ? (Inf, nothing) : rethrow(err)
    end

    function pullback(Δy)
        Δvals = if isinf(y)
            fill!(similar(vals), NaN)
        else
            _, Δsys = f_pb(Δy)
            CRC.unthunk(Δsys).vals
        end
        return (CRC.NoTangent(), Δvals)
    end

    return y, pullback
end

"""
    squared_error(x, y; weights=nothing, rescale=false)

Sum of squared errors, ``L ∝ \\sum_i w_i |y_i - x_i|^2``. The nonnegative
weights ``w_i`` default to 1.

If `rescale=true` then an automated rescaling will be performed whereby ``L ∝
\\min_c \\sum_i w_i |y_i - c x_i|^2``. This can be useful when fitting to
experimental intensities of unknown scale.

In all cases, a normalization is imposed so that ``0 ≤ L ≤ 1``. This
normalization also ensures symmetry in (x, y). If any ``x_i`` or ``y_i`` is NaN,
these terms will be omitted from the sum.
"""
function squared_error(x, y; weights=nothing, rescale=false)
    ty = promote_type(eltype(x), eltype(y))
    w = @something weights fill(one(real(ty)), size(x))
    size(x) == size(y) == size(w) || error("Mismatched input sizes")
    (x, y, w) = flatten_to_vec.((x, y, w))
    all(>=(0), w) || error("Negative weights detected")

    # This functional implementation is AD-friendly
    inds = findall(i -> !isnan(x[i]) && !isnan(y[i]), eachindex(x))
    @views begin
        x² = sum(@. w[inds] * abs2(x[inds]))           # |x|² ≡ ⟨x,x⟩
        y² = sum(@. w[inds] * abs2(y[inds]))           # |y|² ≡ ⟨y,y⟩
        xy = sum(@. w[inds] * conj(x[inds]) * y[inds]) # ⟨x,y⟩
    end

    return if !rescale
        # |y - x|² / 2(|x|²+|y|²)
        1/2 - real(xy) / (x² + y²)
    else
        # |y - cx|² / |y|² where c = ⟨x,y⟩ / |x|²
        1 - abs2(xy) / (x² * y²)
    end
end
