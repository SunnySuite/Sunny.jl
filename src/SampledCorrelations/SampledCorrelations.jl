mutable struct SampledCorrelations{N}
    # 𝒮^{αβ}(q,ω) data and metadata
    const data           :: Array{ComplexF64, 7}                 # Raw SF with sublattice indices (ncorrs × natoms × natoms × latsize × nω)
    const M              :: Union{Nothing, Array{Float64, 7}}    # Running estimate of (nsamples - 1)*σ² (where σ² is the variance of intensities)
    const crystal        :: Crystal                              # Crystal for interpretation of q indices in `data`
    const origin_crystal :: Union{Nothing,Crystal}               # Original user-specified crystal (if different from above) -- needed for FormFactor accounting
    const Δω             :: Float64                              # Energy step size (could make this a virtual property)  
    measure              :: MeasureSpec                          # Observable, correlation pairs, and combiner

    # Trajectory specs
    const measperiod   :: Int                                    # Steps to skip between saving observables (i.e., downsampling factor for trajectories)
    const dt           :: Float64                                # Step size for trajectory integration 
    nsamples           :: Int64                                  # Number of accumulated samples (single number saved as array for mutability)

    # Buffers and precomputed data 
    const samplebuf    :: Array{ComplexF64, 6}                   # Buffer for observables (nobservables × latsize × natoms × nsnapshots)
    const corrbuf      :: Array{ComplexF64, 4}                   # Buffer for correlations (latsize × nω)
    const space_fft!   :: FFTW.AbstractFFTs.Plan                 # Pre-planned lattice FFT for samplebuf
    const time_fft!    :: FFTW.AbstractFFTs.Plan                 # Pre-planned time FFT for samplebuf
    const corr_fft!    :: FFTW.AbstractFFTs.Plan                 # Pre-planned time FFT for corrbuf 
    const corr_ifft!   :: FFTW.AbstractFFTs.Plan                 # Pre-planned time IFFT for corrbuf 
end

function Base.show(io::IO, ::SampledCorrelations{N}) where N
    modename = N == 0 ? "Dipole" : "SU($(N))"
    print(io, "SampledCorrelations{$modename}")
    # TODO: Add correlation info?
end

function Base.show(io::IO, ::MIME"text/plain", sc::SampledCorrelations{N}) where N
    modename = N == 0 ? "Dipole" : "SU($(N))"
    printstyled(io, "SampledCorrelations";bold=true, color=:underline)
    print(io, "{$modename}")
    print(io," ($(Base.format_bytes(Base.summarysize(sc))))\n")
    print(io,"[")
    if size(sc.data)[7] == 1
        printstyled(io,"S(q)";bold=true)
    else
        printstyled(io,"S(q,ω)";bold=true)
        print(io," | nω = $(round(Int, size(sc.data)[7]/2)), Δω = $(round(sc.Δω, digits=4))")
    end
    print(io," | $(sc.nsamples) sample")
    (sc.nsamples > 1) && print(io,"s")
    print(io,"]\n")
    println(io,"Lattice: $(sc.latsize)×$(natoms(sc.crystal))")
    # TODO: Add correlation info?
end

Base.getproperty(sc::SampledCorrelations, sym::Symbol) = sym == :latsize ? size(sc.samplebuf)[2:4] : getfield(sc, sym)

function Base.setproperty!(sc::SampledCorrelations, sym::Symbol, val)
    if sym == :measure
        @assert sc.measure.observables ≈ val.observables "New MeasureSpec must contain identical observables."
        @assert all(x -> x == 1, sc.measure.corr_pairs .== val.corr_pairs) "New MeasureSpec must contain identical correlation pairs."
        setfield!(sc, :measure, val)
    else
        setfield!(sc, sym, val)
    end
end

function clone_correlations(sc::SampledCorrelations{N}) where N
    dims = size(sc.data)[2:4]
    # Avoid copies/deep copies of C-generated data structures
    space_fft! = 1/√prod(dims) * FFTW.plan_fft!(sc.samplebuf, (2,3,4))
    time_fft! = FFTW.plan_fft!(sc.samplebuf, 6)
    corr_fft! = FFTW.plan_fft!(sc.corrbuf, 4)
    corr_ifft! = FFTW.plan_ifft!(sc.corrbuf, 4)
    M = isnothing(sc.M) ? nothing : copy(sc.M)
    return SampledCorrelations{N}(
        copy(sc.data), M, sc.crystal, sc.origin_crystal, sc.Δω, deepcopy(sc.measure), 
        sc.measperiod, sc.dt, sc.nsamples,
        copy(sc.samplebuf), copy(sc.corrbuf), space_fft!, time_fft!, corr_fft!, corr_ifft!
    )
end

"""
    merge_correlations(scs::Vector{SampledCorrelations)

Accumulate a list of `SampledCorrelations` into a single, summary
`SampledCorrelations`. Useful for reducing the results of parallel computations.
"""
function merge_correlations(scs::Vector{SampledCorrelations{N}}) where N
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
    if dt_new/dt < 0.9
        @warn "To satisify specified energy values, the step size adjusted down by more than 10% from a value of dt=$dt to dt=$dt_new"
    end

    return dt_new, measperiod
end


"""
    SampledCorrelations(sys::System{N}; measure, energies, dt=NaN, calculate_errors=false) where N

Create a `SampledCorrelations` for accumulating samples of spin-spin
correlations. Requires a measurement to determine which correlation pairs to
calculate, e.g. `measure=ssf_perp(sys)`.

The stored correlations may either be static (instantaneous), generated from
sampled spin configurations, or dynamic, generated from time-evolved
trajectories. To configure a `SampledCorrelations` for static correlations, set
set `energies=nothing`. To configure a `SampledCorrelations` for dynamic
correlations, provide an evenly-spaced range of energies starting with 0, e.g.
`energies=range(0, 3.0, 100)`. Dynamic correlations also require a time step,
`dt`. See [suggest_timestep](@ref) for help selecting an appropriate value.

"""
function SampledCorrelations(sys::System{N}; measure, energies, dt=NaN, calculate_errors=false) where N

    if isnothing(energies)
        n_all_ω = 1
        measperiod = 1
        dt = NaN
        Δω = NaN
    else
        nω = length(energies)
        n_all_ω = 2(Int(nω) - 1)
        ωmax = energies[end]
        @assert iszero(energies[1]) && ωmax > 0 "`energies` must be a range from 0 to a positive value."
        ΔEs = energies[2:end] - energies[1:end-1]
        @assert all(x -> x ≈ ΔEs[1], ΔEs) "`energies` must be equally spaced."
        dt, measperiod = adjusted_dt_and_downsampling_factor(dt, nω, ωmax)
        Δω = ωmax/(nω-1)
    end

    # Preallocation
    na = natoms(sys.crystal)

    # The sample buffer holds n_non_neg_ω measurements, and the rest is a zero buffer
    measure = isnothing(measure) ? ssf_trace(sys) : measure
    println(typeof(measure))
    num_observables(measure)
    samplebuf = zeros(ComplexF64, num_observables(measure), sys.latsize..., na, n_all_ω)
    corrbuf = zeros(ComplexF64, sys.latsize..., n_all_ω)

    # The output data has n_all_ω many (positive and negative and zero) frequencies
    data = zeros(ComplexF64, num_correlations(measure), na, na, sys.latsize..., n_all_ω)
    M = calculate_errors ? zeros(Float64, size(data)...) : nothing

    # The normalization is defined so that the prod(sys.latsize)-many estimates
    # of the structure factor produced by the correlation conj(space_fft!) * space_fft!
    # are correctly averaged over. The corresponding time-average can't be applied in
    # the same way because the number of estimates varies with Δt. These conventions
    # ensure consistency with this spec:
    # https://sunnysuite.github.io/Sunny.jl/dev/structure-factor.html
    space_fft! = 1/√prod(sys.latsize) * FFTW.plan_fft!(samplebuf, (2,3,4))
    time_fft!  = FFTW.plan_fft!(samplebuf, 6)
    corr_fft!  = FFTW.plan_fft!(corrbuf, 4)
    corr_ifft! = FFTW.plan_ifft!(corrbuf, 4)

    # Initialize nsamples to zero. Make an array so can update dynamically
    # without making struct mutable.
    nsamples = 0 

    # Make Structure factor and add an initial sample
    origin_crystal = isnothing(sys.origin) ? nothing : sys.origin.crystal
    sc = SampledCorrelations{N}(data, M, sys.crystal, origin_crystal, Δω, measure, measperiod, dt, nsamples,
                                samplebuf, corrbuf, space_fft!, time_fft!, corr_fft!, corr_ifft!)

    return sc
end