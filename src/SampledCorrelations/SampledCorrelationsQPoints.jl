mutable struct SampledCorrelationsQPoints{Q <: AbstractQPoints}
    # 𝒮^{αβ}(q,ω) data and metadata for selected finite-grid q indices only.
    const data              :: Array{ComplexF64, 5}              # (ncorrs × natoms × natoms × nq × nω)
    const M                 :: Union{Nothing, Array{Float64, 5}}
    const recipvecs         :: Mat3
    const origin_crystal    :: Crystal
    const Δω                :: Float64
    const qpts              :: Q
    const q_indices         :: Vector{CartesianIndex{3}}
    const request_to_unique :: Vector{Int}
    const request_qabs      :: Vector{Vec3}

    # Observables
    measure         :: MeasureSpec
    const positions :: Vector{Vec3}

    # Trajectory specs
    const integrator :: AbstractIntegrator
    const measperiod :: Int
    nsamples         :: Int64

    # Buffers and precomputed data.
    const samplebuf  :: Array{ComplexF64, 6}
    const corrbuf    :: Matrix{ComplexF64}
    const space_fft! :: FFTW.AbstractFFTs.Plan
    const time_fft!  :: FFTW.AbstractFFTs.Plan
    const corr_fft!  :: FFTW.AbstractFFTs.Plan
    const corr_ifft! :: FFTW.AbstractFFTs.Plan
end

function Base.show(io::IO, ::SampledCorrelationsQPoints)
    print(io, "SampledCorrelationsQPoints")
end

function Base.show(io::IO, ::MIME"text/plain", sc::SampledCorrelationsQPoints)
    nω = round(Int, size(sc.data, 5) / 2)
    printstyled(io, "SampledCorrelationsQPoints"; bold=true, color=:underline)
    println(io, " ($(Base.format_bytes(Base.summarysize(sc))))")
    print(io, "[")
    printstyled(io, "S(q,ω)"; bold=true)
    print(io, " | nω = $nω, Δω = $(round(sc.Δω, digits=4))")
    println(io, " | $(length(sc.qpts.qs)) requested q, $(length(sc.q_indices)) sampled q")
    println(io, " | $(sc.nsamples) $(sc.nsamples > 1 ? "samples" : "sample")]")
end

function Base.setproperty!(sc::SampledCorrelationsQPoints, sym::Symbol, val)
    if sym == :measure
        sc.measure.observables ≈ val.observables || error("Measure observables cannot change")
        validate_formfactors(val)
        sc.measure.offsets ≈ val.offsets || error("Measure offsets cannot change")
        sc.measure.corr_pairs == val.corr_pairs || error("Measure corr_pairs cannot change")
    end
    setfield!(sc, sym, val)
end

function as_qpoints(qpts)
    return qpts isa AbstractQPoints ? qpts : convert(AbstractQPoints, qpts)
end

function selected_wave_vector_info(recipvecs, origin_crystal, dims, qpts::AbstractQPoints)
    qs_reshaped = [recipvecs \ origin_crystal.recipvecs * q for q in qpts.qs]
    unique_lookup = Dict{CartesianIndex{3}, Int}()
    q_indices = CartesianIndex{3}[]
    request_to_unique = Int[]
    request_qabs = Vec3[]

    for q in qs_reshaped
        m = round.(Int, Vec3(dims) .* q)
        idx = CartesianIndex{3}(map(i -> mod(m[i], dims[i]) + 1, (1, 2, 3)))
        source = get(unique_lookup, idx, 0)
        if source == 0
            push!(q_indices, idx)
            source = length(q_indices)
            unique_lookup[idx] = source
        end
        push!(request_to_unique, source)
        push!(request_qabs, recipvecs * (m ./ Vec3(dims)))
    end

    return (; q_indices, request_to_unique, request_qabs)
end

"""
    SampledCorrelationsQPoints(sys::System, qpts; measure, energies, dt,
                               calculate_errors=false, integrator=ImplicitMidpoint())

Accumulate dynamical correlations only for the finite-supercell q-grid points
needed by `qpts`. This is useful for path calculations on large supercells,
where full [`SampledCorrelations`](@ref) storage scales with every point in the
three-dimensional reciprocal grid.
"""
function SampledCorrelationsQPoints(sys::System, qpts; measure, energies, dt, calculate_errors=false, integrator=ImplicitMidpoint())
    recipvecs = sys.crystal.recipvecs
    origin_crystal = orig_crystal(sys)
    qpts = as_qpoints(qpts)

    if isnothing(energies)
        error("SampledCorrelationsQPoints requires a nonempty energy grid")
    else
        nω = length(energies)
        n_all_ω = 2(Int(nω) - 1)
        ωmax = energies[end]
        iszero(energies[1]) && ωmax > 0 || error("`energies` must be a range from 0 to a positive value")
        ΔEs = energies[2:end] - energies[1:end-1]
        all(≈(ΔEs[1]), ΔEs) || error("`energies` must be equally spaced.")
        dt, measperiod = adjusted_dt_and_downsampling_factor(dt, nω, ωmax)
        Δω = ωmax / (nω - 1)
    end

    isnan(integrator.dt) || error("Timestep of `integrator` must be uninitialized.")
    integrator.dt = dt

    measure = isnothing(measure) ? ssf_trace(sys) : measure
    validate_formfactors(measure)

    nunits, nparts = size(measure.offsets)
    positions = [sys.crystal.positions[u] + measure.offsets[u, p] for p in 1:nparts for u in 1:nunits]
    natoms = length(positions)

    q_info = selected_wave_vector_info(recipvecs, origin_crystal, sys.dims, qpts)
    nq = length(q_info.q_indices)
    samplebuf = zeros(ComplexF64, num_observables(measure), sys.dims..., natoms, n_all_ω)
    corrbuf = zeros(ComplexF64, nq, n_all_ω)

    data = zeros(ComplexF64, num_correlations(measure), natoms, natoms, nq, n_all_ω)
    M = calculate_errors ? zeros(Float64, size(data)...) : nothing

    space_fft! = 1 / √prod(sys.dims) * FFTW.plan_fft!(samplebuf, (2, 3, 4))
    time_fft! = FFTW.plan_fft!(samplebuf, 6)
    corr_fft! = FFTW.plan_fft!(corrbuf, 2)
    corr_ifft! = FFTW.plan_ifft!(corrbuf, 2)
    nsamples = 0

    return SampledCorrelationsQPoints(
        data,
        M,
        recipvecs,
        origin_crystal,
        Δω,
        qpts,
        q_info.q_indices,
        q_info.request_to_unique,
        q_info.request_qabs,
        measure,
        positions,
        integrator,
        measperiod,
        nsamples,
        samplebuf,
        corrbuf,
        space_fft!,
        time_fft!,
        corr_fft!,
        corr_ifft!,
    )
end

function available_energies(sc::SampledCorrelationsQPoints; negative_energies=false)
    n_all_ω = size(sc.data, 5)
    n_non_neg_ω = div(n_all_ω, 2) + 1
    ωvals = collect(FFTW.fftfreq(n_all_ω, n_all_ω * sc.Δω))
    ωvals[n_non_neg_ω] *= -1
    return negative_energies ? ωvals : ωvals[1:n_non_neg_ω]
end

function clone_correlations(sc::SampledCorrelationsQPoints)
    space_fft! = 1 / √prod(size(sc.samplebuf)[2:4]) * FFTW.plan_fft!(sc.samplebuf, (2, 3, 4))
    time_fft! = FFTW.plan_fft!(sc.samplebuf, 6)
    corr_fft! = FFTW.plan_fft!(sc.corrbuf, 2)
    corr_ifft! = FFTW.plan_ifft!(sc.corrbuf, 2)
    M = isnothing(sc.M) ? nothing : copy(sc.M)
    return SampledCorrelationsQPoints(
        copy(sc.data),
        M,
        sc.recipvecs,
        sc.origin_crystal,
        sc.Δω,
        sc.qpts,
        copy(sc.q_indices),
        copy(sc.request_to_unique),
        copy(sc.request_qabs),
        deepcopy(sc.measure),
        copy(sc.positions),
        copy(sc.integrator),
        sc.measperiod,
        sc.nsamples,
        copy(sc.samplebuf),
        copy(sc.corrbuf),
        space_fft!,
        time_fft!,
        corr_fft!,
        corr_ifft!,
    )
end

function merge_correlations(scs::Vector{<:SampledCorrelationsQPoints})
    sc_merged = clone_correlations(scs[1])
    μ = zero(sc_merged.data)
    for sc in scs[2:end]
        n = sc_merged.nsamples
        m = sc.nsamples
        @. μ = (n / (n + m)) * sc_merged.data + (m / (n + m)) * sc.data
        if !isnothing(sc_merged.M)
            @. sc_merged.M = (sc_merged.M + n * abs(μ - sc_merged.data)^2) + (sc.M + m * abs(μ - sc.data)^2)
        end
        sc_merged.data .= μ
        sc_merged.nsamples += m
    end
    return sc_merged
end

function new_sample!(sc::SampledCorrelationsQPoints, sys::System)
    (; integrator, samplebuf, measperiod) = sc
    obs = sc.measure.observables
    observables = reshape(obs, size(obs)[1:4]..., :)

    buf_size = size(samplebuf, 6)
    nsnaps = (buf_size ÷ 2) + 1
    samplebuf[:, :, :, :, :, (nsnaps+1):end] .= 0

    trajectory!(samplebuf, sys, integrator, nsnaps, observables; measperiod)
    return nothing
end

function accum_sample!(sc::SampledCorrelationsQPoints; window)
    (; data, M, samplebuf, corrbuf, space_fft!, time_fft!, corr_fft!, corr_ifft!) = sc
    (; corr_pairs) = sc.measure
    natoms = size(samplebuf, 5)
    num_time_offsets = size(samplebuf, 6)
    T = (num_time_offsets ÷ 2) + 1
    time_offsets = FFTW.fftfreq(num_time_offsets, num_time_offsets)
    n_contrib = reshape(T .- abs.(time_offsets), 1, num_time_offsets)
    n_contrib[n_contrib .== 0] .= Inf

    samplebuf .= conj.(samplebuf)
    space_fft! * samplebuf
    samplebuf .= conj.(samplebuf)
    time_fft! * samplebuf

    count = sc.nsamples += 1
    window_func = if window == :cosine
        reshape(cos.(range(0, π, length=num_time_offsets+1)[1:end-1]).^2, 1, num_time_offsets)
    elseif window == :rectangular
        nothing
    else
        error("Unsupported window: $window")
    end

    for j in 1:natoms, i in 1:natoms, (c, (α, β)) in enumerate(corr_pairs)
        for (q, idx) in enumerate(sc.q_indices)
            sample_α = @view samplebuf[α, idx.I..., i, :]
            sample_β = @view samplebuf[β, idx.I..., j, :]
            @. corrbuf[q, :] = conj(sample_α) * sample_β
        end

        corr_ifft! * corrbuf
        corrbuf ./= n_contrib
        if !isnothing(window_func)
            corrbuf .*= window_func
        end
        corr_fft! * corrbuf

        databuf = @view data[c, i, j, :, :]
        if isnothing(M)
            for k in eachindex(databuf)
                diff = corrbuf[k] - databuf[k]
                databuf[k] += diff * (1 / count)
            end
        else
            Mbuf = @view M[c, i, j, :, :]
            for k in eachindex(databuf)
                μ_old = databuf[k]
                databuf[k] += (corrbuf[k] - databuf[k]) / count
                μ = databuf[k]
                Mbuf[k] += real((corrbuf[k] - μ_old) * conj(corrbuf[k] - μ))
            end
        end
    end

    return nothing
end

function add_sample!(sc::SampledCorrelationsQPoints, sys::System; window=:cosine)
    new_sample!(sc, sys)
    accum_sample!(sc; window)
end

function intensities(sc::SampledCorrelationsQPoints; energies=:available, kernel=nothing, kT)
    if !isnothing(kT) && kT <= 0
        error("Positive `kT` required for classical-to-quantum corrections, or set `kT=nothing` to disable.")
    end
    if !isnothing(kernel)
        error("Kernel post-processing not yet available for `SampledCorrelationsQPoints`.")
    end

    (ωs, ωidcs) = if energies == :available
        ωs = available_energies(sc; negative_energies=false)
        (ωs, axes(ωs, 1))
    elseif energies == :available_with_negative
        ωs = available_energies(sc; negative_energies=true)
        (ωs, axes(ωs, 1))
    else
        rounded_energy_information(sc, energies)
    end

    ffs = vec(sc.measure.formfactors[1, :, :])
    data = zeros(eltype(sc.measure), length(ωs), length(sc.qpts.qs))
    NCorr = Val{size(sc.data, 1)}()
    NAtoms = Val{size(sc.data, 2)}()

    for (iq, qabs) in enumerate(sc.request_qabs)
        source_q = sc.request_to_unique[iq]
        prefactors = prefactors_for_phase_averaging(qabs, sc.recipvecs, sc.positions, ffs, NCorr, NAtoms)
        for (n, iω) in enumerate(ωidcs)
            elems = zero(SVector{size(sc.data, 1), ComplexF64})
            for j in 1:size(sc.data, 2), i in 1:size(sc.data, 2)
                elems += (prefactors[i] * conj(prefactors[j])) * SVector{size(sc.data, 1)}(view(sc.data, :, i, j, source_q, iω))
            end
            data[n, iq] = sc.measure.combiner(qabs, elems)
        end
    end

    data .*= det(sc.recipvecs) / det(sc.origin_crystal.recipvecs)
    n_all_ω = size(sc.samplebuf, 6)
    data ./= (n_all_ω * sc.Δω)

    if !isnothing(kT)
        c2q = [iszero(ω) ? 1 : abs((ω / kT) / (1 - exp(-ω / kT))) for ω in ωs]
        for i in axes(data, 2)
            data[:, i] .*= c2q
        end
    end

    return Intensities(sc.origin_crystal, sc.qpts, collect(ωs), data)
end
