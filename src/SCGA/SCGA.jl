"""
    SCGA(sys::System; measure, kT, dq)

Constructs an object to calculate [`intensities_static`](@ref) within the self
consistent gaussian approximation (SCGA). This theory assumes a classical
Boltzmann distribution with temperature `kT`. It is expected to be meangingful
above the ordering temperature, where fluctuations are approximately Gaussian.
If the temperature is not sufficiently high, then `intensities_static` may
report negative energies, which would indicate an instability to magnetic
ordering.

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

    function SCGA(sys::System; measure::Union{Nothing, MeasureSpec}, kT::Float64, dq::Float64)
        measure = @something measure empty_measurespec(sys)
        if size(eachsite(sys)) != size(measure.observables)[2:5]
            error("Size mismatch. Check that measure is built using consistent system.")
        end

        kT > 0 || error("Temperature kT must be positive")
        β = 1 / kT

        qs = make_q_grid(sys, dq)
        Js = [fourier_exchange_matrix(sys; q) for q in qs]

        sublattice_resolved = !allequal(sys.crystal.classes)
        if sublattice_resolved
            λs = find_lagrange_multiplier_opt_sublattice(sys, Js, β)
        else
            λ = find_lagrange_multiplier(sys, Js, β)
            λs = fill(λ, natoms(sys.crystal))
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
function find_lagrange_multiplier(sys, Js, β)
    starting_offset = 0.2
    maxiters = 500
    tol = 1e-10

    evals = Iterators.flatten(eigvals(J) for J in Js)

    Nq = length(Js)
    s² = norm2(sys.κs)

    function f(λ)
        return sum(1 / (λ + ev) for ev in evals) / (β * Nq)
    end
    function J(λ)
        return -sum(1 / (λ + ev)^2 for ev in evals) / (β * Nq)
    end

    λn = starting_offset*0.1/β - minimum(evals)
    for n in 1:maxiters
        λ = λn + (1/J(λn))*(s²-f(λn))
        if abs(λ-λn) < tol
            println("Newton's method converged to within tolerance, $tol, after $n steps.")
            return λ
        else
            λn = λ
        end
    end
end


function find_lagrange_multiplier_opt_sublattice(sys, Js, β)
    tol = 1e-6
    maxiters = 500

    Na = natoms(sys.crystal)
    evals = Iterators.flatten(eigvals(J) for J in Js)
    λ_min, λ_max = extrema(evals)
    λ_init = -λ_min + (λ_max - λ_min) / 2
    λs = λ_init*ones(Float64, Na)
    s² = vec(sys.κs .^ 2)

    function fg!(_, gbuffer, λs)
        fbuffer = 0.0
        if !isnothing(gbuffer)
            gbuffer .= 0
        end

        Λ = diagm(repeat(λs, inner=3))

        for J in Js
            A = β * (J + Λ)
            T = eigen(A)
            eig_vals = T.values
            U = T.vectors
            A⁻¹ = U * Diagonal(inv.(eig_vals)) * U'
            A⁻¹ = reshape(A⁻¹, 3, Na, 3, Na)

            if minimum(eig_vals) < 0
                F = -Inf
            else
                F = (1/2β) * sum(log.(eig_vals))
            end
            # To maximize free energy G = F - λᵢ sᵢ² we should minimize f = -G
            fbuffer += λs' * s² / 2 - F

            if !isnothing(gbuffer)
                for i in 1:Na
                    gbuffer[i] += s²[i] / 2 - real(tr(view(A⁻¹, :, i, :, i))) / 2
                end
            end
        end

        return fbuffer
    end

    # f(λs) = fg!(NaN, nothing, λs)

    options = Optim.Options(; iterations=maxiters, show_trace=false, g_tol=tol)
    result = Optim.optimize(Optim.only_fg!(fg!), λs, Optim.ConjugateGradient(), options)
    min = Optim.minimizer(result)
    return real.(min)
end


function intensities_static(scga::SCGA, qpts)
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
        real(tr(ssf)) < 0 && error("Negative intensity indicates that kT is below ordering")

        map!(corrbuf, measure.corr_pairs) do (α, β)
            ssf[α, β] / Ncells
        end
        intensity[iq] = measure.combiner(q_global, corrbuf)
    end

    return StaticIntensities(sys.crystal, qpts, reshape(intensity, size(qpts.qs)))
end
