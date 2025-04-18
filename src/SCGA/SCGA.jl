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
    Nq :: Int
    quantum_sum_rule :: Bool
    sublattice_resolved :: Bool

    function SCGA(sys::System; measure::Union{Nothing, MeasureSpec}, Nq::Int, quantum_sum_rule=false, sublattice_resolved=true)
        measure = @something measure empty_measurespec(sys)
        if size(eachsite(sys)) != size(measure.observables)[2:5]
            error("Size mismatch. Check that measure is built using consistent system.")
        end
        return new(sys, measure, Nq, quantum_sum_rule, sublattice_resolved)
    end
end


# Computes the Lagrange multiplier for the standard SCGA approach with a common
# Lagrange multiplier for all sublattices.
function find_lagrange_multiplier(scga::SCGA, β)
    starting_offset = 0.2
    maxiters = 500
    tol = 1e-10

    (; sys, quantum_sum_rule, Nq) = scga
    dq = 1/Nq;
    Na = natoms(sys.crystal)
    bond_counter = zeros(Float64,3)
    for i in 1:Na
        for coupling in sys.interactions_union[i].pair
            (; isculled, bond) = coupling
            isculled && break
            bond_counter += abs.(bond.n)
        end
    end
    qarray = -0.5: dq : 0.5-dq
    qarrays = []
    for i in 1:3
        if bond_counter[i] == 0
            push!(qarrays,[0])
        else
            push!(qarrays,qarray)
        end
    end
    q = [[qx, qy, qz] for qx in qarrays[1], qy in qarrays[2], qz in qarrays[3]]
    Jq_array = [fourier_exchange_matrix(sys; k=q_in) for q_in in q]
    #TODO throw error if spins different
    if quantum_sum_rule
        S_sq = sum(sys.κs .* (sys.κs.+1)) / Na
    else
        S_sq = sum(sys.κs.^2) / Na
    end

    eig_vals = zeros(3Na, length(Jq_array))
    for j in eachindex(Jq_array)
         eig_vals[:, j] .= eigvals(Jq_array[j])
    end

    function f(λ)
        sum_term = sum(1 ./ (λ .+ β * eig_vals ))
        return (1/(Na*length(q))) * sum_term
    end
    function J(λ)
        sum_term = sum((1 ./ (λ .+ β * eig_vals)).^2)
        return -(1/(Na*length(q))) * sum_term
    end

    lower = - β * minimum(eig_vals)
    λn = starting_offset*0.1 + lower # Make more robust - regularized! Check optim.jl

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

# Computes the static structure factor in the standard SCGA approach, with a
# single Lagrange multiplier, over the qpts specified.

function intensities_static_single(scga::SCGA, qpts; β)
    qpts = convert(AbstractQPoints, qpts)
    (; sys, measure) = scga
    Na = natoms(sys.crystal)
    Nobs = num_observables(measure)
    λ = find_lagrange_multiplier(scga, β)
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
        J_mat = fourier_exchange_matrix(sys; k=q_reshaped)
        inverted_matrix = inv(I(3*Na)*λ[1] + β*J_mat) # this is [(Iλ+J(q))^-1]^αβ_μν
        # inverted_matrix = diagm(1 ./ (λ[1] .+ β .* eigvals(J_mat)))
        for i in 1:Na, j in 1:Na
            intensitybuf += pref[1,i]*conj(pref[1,j])*inverted_matrix[1+3(i-1):3+3(i-1),1+3(j-1):3+3(j-1)]
            # TODO allow different form factor for each observable
        end
        map!(corrbuf, measure.corr_pairs) do (α, β)
            intensitybuf[α,β]
        end
        intensity[iq] = measure.combiner(q_global, corrbuf)
    end
    if extrema(intensity)[1] < -1e-2
        @warn "Warning: negative intensities! kT is probably below the ordering temperature."
        # TODO Throw an error. This is for diagnostic purposes.
    end
    return StaticIntensities(sys.crystal, qpts, reshape(intensity,size(qpts.qs)))
end

# Computes the static structure factor in the sublattice resolved SCGA method,
# over the qpts specified.

function intensities_static_sublattice(scga::SCGA, qpts; β, λs_init)
    qpts = convert(AbstractQPoints, qpts)
    (; sys, measure) = scga
    Na = natoms(sys.crystal)
    Nobs = num_observables(measure)
    λs = find_lagrange_multiplier_opt_sublattice(scga, λs_init, β)
    intensity = zeros(eltype(measure),length(qpts.qs))
    Ncorr = length(measure.corr_pairs)
    r = sys.crystal.positions

    for (iq, q) in enumerate(qpts.qs)
        pref = zeros(ComplexF64, Nobs, Na)
        corrbuf = zeros(ComplexF64, Ncorr)
        intensitybuf = zeros(eltype(measure),3,3)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = sys.crystal.recipvecs * q
        for i in 1:Na, μ in 1:3
            ff = get_swt_formfactor(measure, μ, i)
            pref[μ, i] = exp(-2π * im * dot(q_reshaped, r[i])) * compute_form_factor(ff, norm2(q_global))
        end
        J_mat = fourier_exchange_matrix(sys; k=q_reshaped)
        Λ =  diagm(repeat(λs, inner=3))
        inverted_matrix = (inv(β*Λ + β*J_mat)) # this is [(Iλ+J(q))^-1]^αβ_μν
        for i in 1:Na, j in 1:Na
            intensitybuf += pref[1,i]*conj(pref[1,j])*inverted_matrix[1+3(i-1):3+3(i-1),1+3(j-1):3+3(j-1)]
            # TODO allow different form factor for each observable
        end
        map!(corrbuf, measure.corr_pairs) do (α, β)
            intensitybuf[α,β]
        end
        intensity[iq] = measure.combiner(q_global, corrbuf)
    end
    if extrema(intensity)[1] < -1e-2
        @warn "Warning: negative intensities! kT is probably below the ordering temperature."
        # TODO Throw an error. This is for diagnostic purposes.
    end
    println("Optimized Lagrange multipliers: $λs")
    return StaticIntensities(sys.crystal, qpts, reshape(intensity,size(qpts.qs)))
end

function intensities_static(scga::SCGA, qpts; kT, λs_init=nothing)
    kT > 0 || error("Temperature kT must be positive")
    β = 1 / kT
    if scga.sublattice_resolved == true
        return intensities_static_sublattice(scga::SCGA, qpts; β, λs_init)
    else
        return intensities_static_single(scga::SCGA, qpts; β)
    end
end

# Computes the Lagrange multiplier for the sublattice resolved SCGA method.

function find_lagrange_multiplier_opt_sublattice(scga, λs, β)
    tol = 1e-6
    maxiters = 500

    (; sys, quantum_sum_rule, Nq) = scga

    if quantum_sum_rule
        S_sq =vec(sys.κs .* (sys.κs .+ 1))
    else
        S_sq = vec(sys.κs.^2)
    end
    Na = natoms(sys.crystal)
    dq = 1/Nq;

    bond_counter = zeros(Float64,3)
    for i in 1:Na
        for coupling in sys.interactions_union[i].pair
            (; isculled, bond) = coupling
            isculled && break
            bond_counter += abs.(bond.n)
        end
    end

    qarray = -0.5: dq : 0.5-dq
    qarrays = []
    for i in 1:3
        if bond_counter[i] == 0
            push!(qarrays, [0])
        else
            push!(qarrays, qarray)
        end
    end

    q = [[qx, qy, qz] for qx in qarrays[1], qy in qarray[2], qz in qarray[3]]
    N = length(q)

    function f(λs)
        Λ =  diagm(repeat(λs, inner=3))
        A_array = [β*Sunny.fourier_exchange_matrix(sys; k=q_in) .+  β*Λ for q_in in q]
        eig_vals = zeros(3Na,length(A_array))
        Us = zeros(ComplexF64,3Na,3Na,length(A_array))
        for j in 1:length(A_array)
            T = eigen(A_array[j])
            eig_vals[:,j] .= T.values
            Us[:,:,j] .= T.vectors
        end
        if minimum(eig_vals) < 0
            F =  -Inf
        else
            F =  (1/2β) * sum(log.(eig_vals))
        end
        G = F - (1/2β) * N * sum(λs .* S_sq)
        return -G
    end

    function fg!(fbuffer,gbuffer,λs)
        Λ =  diagm(repeat(λs, inner=3))
        A_array = [β*fourier_exchange_matrix(sys; k=q_in) .+  β*Λ for q_in in q]
        eig_vals = zeros(3Na,length(A_array))
        Us = zeros(ComplexF64,3Na,3Na,length(A_array))
        for j in 1:length(A_array)
            T = eigen(A_array[j])
            eig_vals[:,j] .= T.values
            Us[:,:,j] .= T.vectors
        end
        if gbuffer !== nothing
            gradF = zeros(ComplexF64,Na)
            for i in 1:Na
                gradλ = diagm(zeros(ComplexF64,3Na))
                gradλ[3i-2:3i, 3i-2:3i] = diagm([1,1,1])
                gradF[i] =0.5sum([tr(diagm(1 ./ eig_vals[:,j]) * Us[:,:,j]'*gradλ*Us[:,:,j]) for j in 1:length(A_array)])
            end
            gradG = gradF - 0.5*N*S_sq
            gbuffer .= -real(gradG)
        end
        if !isnothing(fbuffer)
            if minimum(eig_vals) < 0
                F = -Inf
            else
                F = (1/2β) * sum(log.(eig_vals))
            end
            G = F - (N/2) * sum(λs .* S_sq)
            fbuffer= -G
        end
    end

    if isnothing(λs)
        println("No user provided initial guess for the Lagrange multipliers. Determining a sensible starting point from the interaction matrix.")
        A_array = [fourier_exchange_matrix(sys; k=q_in) for q_in in q]
        eig_vals = zeros(3Na,length(A_array))
        Us = zeros(ComplexF64,3Na,3Na,length(A_array))
        for j in 1:length(A_array)
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
        λs = λ_init*ones(Float64,Na)
    else
        println("Using user provided initial starting point for the optimization of the Lagrange multipliers.")
        if f(λs) > 1e7
            println("Matrix is not positive definite. Shifting Lagrange multipliers to find a better starting point.")
            A_array = [fourier_exchange_matrix(sys; k=q_in)  for q_in in q]
            eig_vals = zeros(3Na,length(A_array))
            Us = zeros(ComplexF64,3Na,3Na,length(A_array))
            for j in 1:length(A_array)
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
    upper = Inf
    lower = -Inf
    options = Optim.Options(; iterations=maxiters, show_trace=true, g_tol=tol)
    result = Optim.optimize(Optim.only_fg!(fg!), λs, Optim.ConjugateGradient(), options)
    min = Optim.minimizer(result)
    return real.(min)
end
