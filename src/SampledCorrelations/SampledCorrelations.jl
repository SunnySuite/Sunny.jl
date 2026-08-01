mutable struct SampledCorrelations
    # 𝒮^{αβ}(q,ω) data and metadata
    const data           :: Array{ComplexF64, 7}              # Raw correlation data (ncorrs × natoms × natoms × qdims × nω)
    const M              :: Union{Nothing, Array{Float64, 7}} # Running estimate of (nsamples-1)*σ² (where σ² is the variance of intensities)
    const recipvecs      :: Mat3                              # Reshaped recipvecs, sys.crystal.recipvecs
    const origin_crystal :: Crystal                           # Original chemical cell, orig_crystal(sys)
    const Δω             :: Float64                           # Energy step size 

    # Observables
    measure            :: MeasureSpec                         # Mutation must only change formfactors and combiner
    const positions    :: Vector{Vec3}                        # Fractional positions in reshaped cell (natoms = nunits × nparts)

    # Trajectory specs
    const integrator   :: AbstractIntegrator                  # Integrator for calculating sample trajectories.
    const measperiod   :: Int                                 # Steps to skip between saving observables (i.e., downsampling factor for trajectories)
    nsamples           :: Int64                               # Number of accumulated samples (single number saved as array for mutability)

    # Buffers and precomputed data 
    const samplebuf    :: Array{ComplexF64, 6}                # Buffer for observables (nobservables × qdims × natoms × nsnapshots)
    const corrbuf      :: Array{ComplexF64, 4}                # Buffer for correlations (qdims × nω)
    const space_fft!   :: FFTW.AbstractFFTs.Plan              # Pre-planned lattice FFT for samplebuf
    const time_fft!    :: FFTW.AbstractFFTs.Plan              # Pre-planned time FFT for samplebuf
    const corr_fft!    :: FFTW.AbstractFFTs.Plan              # Pre-planned time FFT for corrbuf 
    const corr_ifft!   :: FFTW.AbstractFFTs.Plan              # Pre-planned time IFFT for corrbuf 
end

function Base.show(io::IO, ::SampledCorrelations)
    print(io, "SampledCorrelations")
    # TODO: Add correlation info?
end

function Base.show(io::IO, ::MIME"text/plain", sc::SampledCorrelations)
    nω = round(Int, size(sc.data)[7]/2)
    qdims = size(sc.data[4:6])
    printstyled(io, "SampledCorrelations"; bold=true, color=:underline)
    println(io," ($(Base.format_bytes(Base.summarysize(sc))))")
    print(io,"[")
    printstyled(io,"S(q,ω)"; bold=true)
    print(io," | nω = $nω, Δω = $(round(sc.Δω, digits=4))")
    println(io," | $(sc.nsamples) $(sc.nsamples > 1 ? "samples" : "sample")]")
    println(io,"Lattice: $qdims × $(size(sc.data, 2))")
end

function Base.setproperty!(sc::SampledCorrelations, sym::Symbol, val)
    if sym == :measure
        # The sc.data is fixed, so the new measure cannot change observables,
        # offsets, or corre_pairs.
        sc.measure.observables ≈ val.observables || error("Measure observables cannot change")
        validate_formfactors(val)
        sc.measure.offsets ≈ val.offsets || error("Measure offsets cannot change")
        sc.measure.corr_pairs == val.corr_pairs || error("Measure corr_pairs cannot change")
    end
    setfield!(sc, sym, val)
end

function validate_formfactors(measure)
    # Dependence on observable dependence (index 1) disallowed for efficiency.
    all(allequal, eachslice(measure.formfactors; dims=(2, 3))) ||
        error("Observable-dependent form factors disallowed for efficiency")
end

"""
    SampledCorrelationsStatic(sys::System; measure)

An object to accumulate samples of static pair correlations. It is similar to
[`SampledCorrelations`](@ref), but no time-integration will be performed on
calls to [`add_sample!`](@ref). The resulting object can be used with
[`intensities_static`](@ref) to calculate statistics from the classical
Boltzmann distribution. Dynamical [`intensities`](@ref) data, however, will be
unavailable. Similarly, classical-to-quantum corrections that rely on the
excitation spectrum cannot be performed.
"""
struct SampledCorrelationsStatic
    parent :: SampledCorrelations
end

function SampledCorrelationsStatic(sys::System; measure, calculate_errors=false)
    parent = SampledCorrelations(sys; measure, energies=nothing, dt=NaN, calculate_errors)
    return SampledCorrelationsStatic(parent)
end

function Base.show(io::IO, ::SampledCorrelationsStatic)
    print(io, "SampledCorrelationsStatic")
end

function Base.show(io::IO, ::MIME"text/plain", sc::SampledCorrelationsStatic)
    nsamples = sc.parent.nsamples
    qdims = size(sc.parent.data[4:6])
    printstyled(io, "SampledCorrelationsStatic"; bold=true, color=:underline)
    println(io," ($(Base.format_bytes(Base.summarysize(sc))))")
    print(io,"[")
    printstyled(io,"S(q)"; bold=true)
    println(io," | $nsamples $(nsamples > 1 ? "samples" : "sample")]")
    println(io,"Lattice: $qdims × $(size(sc.parent.data, 2))")
end

"""
    clone_correlations(sc::SampledCorrelations)

Create a copy of a `SampledCorrelations`.
"""
function clone_correlations(sc::SampledCorrelations)
    dims = size(sc.data)[2:4]
    # Avoid copies/deep copies of C-generated data structures
    space_fft! = 1/√prod(dims) * FFTW.plan_fft!(sc.samplebuf, (2,3,4))
    time_fft! = FFTW.plan_fft!(sc.samplebuf, 6)
    corr_fft! = FFTW.plan_fft!(sc.corrbuf, 4)
    corr_ifft! = FFTW.plan_ifft!(sc.corrbuf, 4)
    M = isnothing(sc.M) ? nothing : copy(sc.M)
    return SampledCorrelations(
        copy(sc.data), M, sc.recipvecs, sc.origin_crystal, sc.Δω, sc.measure,
        copy(sc.positions), copy(sc.integrator), sc.measperiod, sc.nsamples,
        copy(sc.samplebuf), copy(sc.corrbuf), space_fft!, time_fft!, corr_fft!, corr_ifft!
    )
end

function clone_correlations(sc::SampledCorrelationsStatic)
    return SampledCorrelationsStatic(clone_correlations(sc.parent))
end

"""
    merge_correlations(scs::Vector{SampledCorrelations)

Accumulate a list of `SampledCorrelations` into a single, summary
`SampledCorrelations`. Useful for reducing the results of parallel computations.
"""
function merge_correlations(scs::Vector{SampledCorrelations})
    sc_merged = clone_correlations(scs[1])
    μ = zero(sc_merged.data)
    for sc in scs[2:end]
        n = sc_merged.nsamples
        m = sc.nsamples
        @. μ = (n/(n+m))*sc_merged.data + (m/(n+m))*sc.data
        if !isnothing(sc_merged.M)
            @. sc_merged.M = (sc_merged.M + n*abs(μ - sc_merged.data)^2) + (sc.M + m*abs(μ - sc.data)^2)
        end
        sc_merged.data .= μ
        sc_merged.nsamples += m
    end
    sc_merged
end

function merge_correlations(scs::Vector{SampledCorrelationsStatic})
    sc_merged = merge_correlations([sc.parent for sc in scs]) 
    return SampledCorrelationsStatic(sc_merged)
end


# Determine a step size and down sampling factor that results in precise
# satisfaction of user-specified energy values.
function adjusted_dt_and_downsampling_factor(dt, nω, ωmax)
    @assert π/dt > ωmax "Desired `ωmax` not possible with specified `dt`. Choose smaller `dt` value."

    # Assume nω is the number of non-negative frequencies and determine total
    # number of frequency bins.
    n_all_ω = 2(Int64(nω) - 1)

    # Find downsampling factor for the given `dt` that yields an `ωmax` higher
    # than or equal to given `ωmax`. Then adjust `dt` down so that specified
    # `ωmax` is satisfied exactly.
    Δω = ωmax/(nω-1)
    measperiod = ceil(Int, π/(dt * ωmax))
    dt_new = 2π/(Δω*measperiod*n_all_ω)

    # Warn the user if `dt` required drastic adjustment, which will slow
    # simulations.
    # if dt_new/dt < 0.9
    #     @warn "To satisify specified energy values, the step size adjusted down by more than 10% from a value of dt=$dt to dt=$dt_new"
    # end

    return dt_new, measperiod
end


function to_reshaped_rlu(sc::SampledCorrelations, q)
    return sc.recipvecs \ sc.origin_crystal.recipvecs * q
end

"""
    SampledCorrelations(sys::System; measure, energies, dt)

An object to accumulate samples of dynamical pair correlations. The `measure`
argument specifies a pair correlation type, e.g. [`ssf_perp`](@ref). The
`energies` must be evenly-spaced and starting from 0, e.g. `energies = range(0,
3, 100)`. Select the integration time-step `dt` according to accuracy and speed
considerations. [`suggest_timestep`](@ref) can help in selecting an appropriate
value.

Dynamical correlations will be accumulated through calls to
[`add_sample!`](@ref), which expects a spin configuration in thermal
equilibrium. A classical spin dynamics trajectory will be simulated of
sufficient length to achieve the target energy resolution. The resulting data
can can then be extracted as pair-correlation [`intensities`](@ref) with
appropriate classical-to-quantum correction factors. See also
[`intensities_static`](@ref), which integrates over energy.
"""
function SampledCorrelations(sys::System; measure, energies, dt, calculate_errors=false, integrator=ImplicitMidpoint())
    recipvecs = sys.crystal.recipvecs
    origin_crystal = orig_crystal(sys)
    if isnothing(energies)
        n_all_ω = 1
        measperiod = 1
        isnan(dt) || error("Must use dt=NaN when energies=nothing")
        Δω = NaN
    else
        nω = length(energies)
        n_all_ω = 2(Int(nω) - 1)
        ωmax = energies[end]
        iszero(energies[1]) && ωmax > 0 || error("`energies` must be a range from 0 to a positive value")
        ΔEs = energies[2:end] - energies[1:end-1]
        all(≈(ΔEs[1]), ΔEs) || error("`energies` must be equally spaced.")
        dt, measperiod = adjusted_dt_and_downsampling_factor(dt, nω, ωmax)
        Δω = ωmax/(nω-1)
    end

    isnan(integrator.dt) || error("Timestep of `integrator` must be uninitialized.")
    integrator.dt = dt

    measure = isnothing(measure) ? ssf_trace(sys) : measure
    validate_formfactors(measure)

    # Flatten the "bare" observable positions of the cell. Column-major order
    # (atom index varies fastest) matches the layout of `measure.offsets`.
    nunits, nparts = size(measure.offsets)
    positions = [sys.crystal.positions[u] + measure.offsets[u, p] for p in 1:nparts for u in 1:nunits]
    natoms = length(positions)

    samplebuf = zeros(ComplexF64, num_observables(measure), sys.dims..., natoms, n_all_ω)
    corrbuf = zeros(ComplexF64, sys.dims..., n_all_ω)

    # The output data has n_all_ω many (positive and negative and zero) frequencies
    data = zeros(ComplexF64, num_correlations(measure), natoms, natoms, sys.dims..., n_all_ω)
    M = calculate_errors ? zeros(Float64, size(data)...) : nothing

    # The normalization is defined so that the prod(sys.dims)-many estimates of
    # the structure factor produced by the correlation conj(space_fft!) *
    # space_fft! are correctly averaged over. The corresponding time-average
    # can't be applied in the same way because the number of estimates varies
    # with Δt. These conventions ensure consistency with this spec:
    # https://sunnysuite.github.io/Sunny.jl/dev/structure-factor.html
    space_fft! = 1/√prod(sys.dims) * FFTW.plan_fft!(samplebuf, (2,3,4))
    time_fft!  = FFTW.plan_fft!(samplebuf, 6)
    corr_fft!  = FFTW.plan_fft!(corrbuf, 4)
    corr_ifft! = FFTW.plan_ifft!(corrbuf, 4)

    # Initialize nsamples to zero. Make an array so can update dynamically
    # without making struct mutable.
    nsamples = 0 

    sc = SampledCorrelations(data, M, recipvecs, origin_crystal, Δω, measure,
                             positions, integrator, measperiod, nsamples, samplebuf,
                             corrbuf, space_fft!, time_fft!, corr_fft!, corr_ifft!)

    return sc
end
