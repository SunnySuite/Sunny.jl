# The default value `sublattice_resolved=true` is generally preferred, and imposes
# one global sum-rule constraint per spin-sublattice. Selecting instead
# `sublattice_resolved=false` will employ only a _single_ constraint that is
# averaged over sublattices. This will be faster, but risks producing incorrect
# results when the sublattices are symmetry-distinct.


"""
    SCGA(sys::System; measure, Nq, quantum_sum_rule=false)

Constructs an object to calculate [`intensities_static`](@ref) within the self
consistent gaussian approximation (SCGA). This approximation assumes a classical
Boltzmann distribution, and is expected to be meangingful above the ordering
temperature, where fluctuations are approximately Gaussian. If the temperature
is not sufficiently high, then `intensities_static` may report negative
energies, which would indicate an instability to magnetic ordering.

Currently only `:dipole` and `:dipole_uncorrected` system modes are supported.

The theory of SCGA approximates local spin magnitude constraints with a _weaker_
global constraint condition. This global constraint is implemented as a sum
rule, expressed as an integral over Fourier modes. This integral is approximated
as a discrete sum over `Nq^3` wavevectors for the provided integer `Nq`.

By default, each classical spin dipole is assumed to have magnitude ``s`` that
matches the [`Moment`](@ref) specification. Selecting `quantum_sum_rule=true`
will modify this magnitude to ``√s(s+1)``.
"""
struct SCGA
    sys :: System
    measure :: MeasureSpec
    quantum_sum_rule :: Bool
    sublattice_resolved :: Bool

    qs :: Vector{Vec3}
    Js :: Vector{HermitianC64}
    eigvals :: Vector{Vector{Float64}}
    eigvecs :: Vector{Matrix{ComplexF64}}

    function SCGA(sys::System; measure::Union{Nothing, MeasureSpec}, Nq::Int, quantum_sum_rule=false, sublattice_resolved=true)
        measure = @something measure empty_measurespec(sys)
        if size(eachsite(sys)) != size(measure.observables)[2:5]
            error("Size mismatch. Check that measure is built using consistent system.")
        end

        qs = make_q_grid(sys, Nq)
        Js = [fourier_exchange_matrix(sys; k) for k in qs]
        eigs = map(eigen, Js)
        eigvals = [e.values for e in eigs]
        eigvecs = [e.vectors for e in eigs]
        return new(sys, measure, quantum_sum_rule, sublattice_resolved, qs, Js, eigvals, eigvecs)
    end
end

function make_q_grid(sys, Nq)
    Na = natoms(sys.crystal)
    dq = 1/Nq;
    wraps = [false, false, false]
    for i in 1:Na
        for coupling in sys.interactions_union[i].pair
            (; isculled, bond) = coupling
            isculled && break
            @. wraps = wraps || !iszero(bond.n)
        end
    end
    qarrays = [w ? (-1/2 : dq : 1/2-dq) : [0] for w in wraps]
    return vec(Vec3.(Iterators.product(qarrays...)))
end

# Computes the Lagrange multiplier for the standard SCGA approach with a common
# Lagrange multiplier for all sublattices.
function find_lagrange_multiplier(scga::SCGA, β)
    starting_offset = 0.2
    maxiters = 500
    tol = 1e-10

    (; sys, quantum_sum_rule, qs, eigvals) = scga
    eigvals = Iterators.flatten(eigvals)
    Na = natoms(sys.crystal)
    Nq = length(qs)

    if quantum_sum_rule
        S_sq = sum(κ * (κ + 1) for κ in sys.κs) / Na
    else
        S_sq = norm2(sys.κs) / Na
    end

    function f(λ)
        return sum(1 / (λ + ev) for ev in eigvals) / (β * Na * Nq)
    end
    function J(λ)
        return -sum(1 / (λ + ev)^2 for ev in eigvals) / (β * Na * Nq)
    end

    λn = starting_offset*0.1/β - minimum(eigvals)
    for n in 1:maxiters
        λ = λn + (1/J(λn))*(S_sq-f(λn))
        if abs(λ-λn) < tol
            println("Newton's method converged to within tolerance, $tol, after $n steps.")
            return λ
        else
            λn = λ
        end
    end
end

function intensities_static_aux(scga::SCGA, qpts; β, Λ)
    qpts = convert(AbstractQPoints, qpts)
    (; sys, measure) = scga
    Na = natoms(sys.crystal)
    Nobs = num_observables(measure)
    intensity = zeros(eltype(measure),length(qpts.qs))
    Ncorr = length(measure.corr_pairs)
    r = sys.crystal.positions

    for (iq, q) in enumerate(qpts.qs)
        pref = zeros(ComplexF64, Nobs, Na)
        corrbuf = zeros(ComplexF64, Ncorr)
        intensitybuf = zeros(eltype(measure), 3, 3)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = sys.crystal.recipvecs * q
        for i in 1:Na, μ in 1:3
            ff = get_swt_formfactor(measure, μ, i)
            pref[μ, i] = exp(-2π * im * dot(q_reshaped, r[i])) * compute_form_factor(ff, norm2(q_global))
        end
        Jq = fourier_exchange_matrix(sys; k=q_reshaped)         # FIXME: q not q_reshaped
        inverted_matrix = inv(β*Λ + β*Jq) # this is [(Iλ+J(q))^-1]^αβ_μν
        inverted_matrix = reshape(inverted_matrix, 3, Na, 3, Na)

        for i in 1:Na, j in 1:Na
            intensitybuf += pref[1,i]*conj(pref[1,j])*view(inverted_matrix, :, i, :, j)
        end
        map!(corrbuf, measure.corr_pairs) do (α, β)
            intensitybuf[α,β]
        end
        intensity[iq] = measure.combiner(q_global, corrbuf)
    end
    if extrema(intensity)[1] < -1e-2
        error("Negative intensity indicates that kT is below ordering")
    end
    return StaticIntensities(sys.crystal, qpts, reshape(intensity,size(qpts.qs)))
end


function intensities_static(scga::SCGA, qpts; kT, λs_init=nothing)
    kT > 0 || error("Temperature kT must be positive")
    β = 1 / kT
    if scga.sublattice_resolved == true
        λs = find_lagrange_multiplier_opt_sublattice(scga, λs_init, β)
        # println("Optimized Lagrange multipliers: $λs")
        Λ = Diagonal(repeat(λs, inner=3))
    else
        Na = natoms(scga.sys.crystal)
        λ = find_lagrange_multiplier(scga, β)
        Λ = Diagonal(fill(λ, 3Na))
    end

    return intensities_static_aux(scga::SCGA, qpts; β, Λ)
end


# Computes the Lagrange multiplier for the sublattice resolved SCGA method.

function find_lagrange_multiplier_opt_sublattice(scga, λs, β)
    tol = 1e-6
    maxiters = 500

    (; sys, quantum_sum_rule, qs, Js) = scga

    if quantum_sum_rule
        S_sq = vec(sys.κs .* (sys.κs .+ 1))
    else
        S_sq = vec(sys.κs.^2)
    end
    Na = natoms(sys.crystal)
    Nq = length(qs)

    function fg!(fbuffer, gbuffer, λs)
        Λ = diagm(repeat(λs, inner=3))
        A_array = [β*J .+  β*Λ for J in Js]
        eig_vals = zeros(3Na,length(A_array))
        Us = zeros(ComplexF64,3Na,3Na,length(A_array))
        for j in eachindex(A_array)
            T = eigen(A_array[j])
            eig_vals[:,j] .= T.values
            Us[:,:,j] .= T.vectors
        end
        if !isnothing(gbuffer)
            gradF = zeros(ComplexF64,Na)
            for i in 1:Na
                gradλ = diagm(zeros(ComplexF64,3Na))
                gradλ[3i-2:3i, 3i-2:3i] = diagm([1,1,1])
                gradF[i] = 0.5sum([tr(diagm(1 ./ eig_vals[:,j]) * Us[:,:,j]'*gradλ*Us[:,:,j]) for j in eachindex(A_array)])
            end
            gradG = gradF - 0.5*Nq*S_sq
            gbuffer .= -real(gradG)
        end
        if !isnothing(fbuffer)
            if minimum(eig_vals) < 0
                F = -Inf
            else
                F = (1/2β) * sum(log.(eig_vals))
            end
            G = F - (Nq/2) * sum(λs .* S_sq)
            return -G # Target is to maximize free energy G
        end
    end

    f(λs) = fg!(NaN, nothing, λs)

    if isnothing(λs)
        println("No user provided initial guess for the Lagrange multipliers. Determining a sensible starting point from the interaction matrix.")
        A_array = Js
        eig_vals = zeros(3Na,length(A_array))
        Us = zeros(ComplexF64,3Na,3Na,length(A_array))
        for j in eachindex(A_array)
            T = eigen(A_array[j])
            eig_vals[:,j] .= T.values
            Us[:,:,j] .= T.vectors
        end
        extreme_eigvals = extrema(eig_vals)
        if extreme_eigvals[1] < 0
            pos_shift = -extreme_eigvals[1]
        else
            pos_shift = 0
        end
        λ_init = pos_shift + (extreme_eigvals[2]-extreme_eigvals[1])/2
        λs = λ_init*ones(Float64, Na)
    else
        println("Using user provided initial starting point for the optimization of the Lagrange multipliers.")
        if f(λs) > 1e7
            println("Matrix is not positive definite. Shifting Lagrange multipliers to find a better starting point.")
            A_array = Js
            eig_vals = zeros(3Na,length(A_array))
            Us = zeros(ComplexF64,3Na,3Na,length(A_array))
            for j in eachindex(A_array)
                T = eigen(A_array[j])
                eig_vals[:,j] .= T.values
                Us[:,:,j] .= T.vectors
            end
            extreme_eigvals = extrema(eig_vals)
            shift = (extreme_eigvals[2] - extreme_eigvals[1])/2
            while f(λs) > 1e7
                λs += shift*ones(Float64,Na)
            end
        end
    end

    options = Optim.Options(; iterations=maxiters, show_trace=true, g_tol=tol)
    result = Optim.optimize(Optim.only_fg!(fg!), λs, Optim.ConjugateGradient(), options)
    min = Optim.minimizer(result)
    return real.(min)
end
