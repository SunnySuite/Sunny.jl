"""
    SampledCorrelations

Basic data type for storing sampled correlation data. A `SampleCorrelations` is
initialized by calling either [`dynamic_correlations`](@ref) or
[`instant_correlations`](@ref).
"""
struct SampledCorrelations{N}
    # 𝒮^{αβ}(q,ω) data and metadata
    data           :: Array{ComplexF64, 7}                 # Raw SF with sublattice indices (ncorrs × natoms × natoms × latsize × nω)
    M              :: Union{Nothing, Array{Float64, 7}}    # Running estimate of (nsamples - 1)*σ² (where σ² is the variance of intensities)
    crystal        :: Crystal                              # Crystal for interpretation of q indices in `data`
    origin_crystal :: Union{Nothing,Crystal}               # Original user-specified crystal (if different from above) -- needed for FormFactor accounting
    Δω             :: Float64                              # Energy step size (could make this a virtual property)  
    measure        :: MeasureSpec                          # Observable, correlation pairs, and combiner

    # Trajectory specs
    measperiod   :: Int                                    # Steps to skip between saving observables (i.e., downsampling factor for trajectories)
    dt           :: Float64                                # Step size for trajectory integration 
    nsamples     :: Array{Int64, 1}                        # Number of accumulated samples (single number saved as array for mutability)

    # Buffers and precomputed data 
    samplebuf    :: Array{ComplexF64, 6}                   # Buffer for observables (nobservables × latsize × natoms × nsnapshots)
    corrbuf      :: Array{ComplexF64, 4}                   # Buffer for correlations (latsize × nω)
    space_fft!   :: FFTW.AbstractFFTs.Plan                 # Pre-planned lattice FFT for samplebuf
    time_fft!    :: FFTW.AbstractFFTs.Plan                 # Pre-planned time FFT for samplebuf
    corr_fft!    :: FFTW.AbstractFFTs.Plan                 # Pre-planned time FFT for corrbuf 
    corr_ifft!   :: FFTW.AbstractFFTs.Plan                 # Pre-planned time IFFT for corrbuf 
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
    print(io," | $(sc.nsamples[1]) sample")
    (sc.nsamples[1] > 1) && print(io,"s")
    print(io,"]\n")
    println(io,"Lattice: $(sc.latsize)×$(natoms(sc.crystal))")
    # TODO: Add correlation info?
end

Base.getproperty(sc::SampledCorrelations, sym::Symbol) = sym == :latsize ? size(sc.samplebuf)[2:4] : getfield(sc, sym)

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
        sc.measperiod, sc.dt, copy(sc.nsamples),
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
        n = sc_merged.nsamples[1] 
        m = sc.nsamples[1]
        @. μ = (n/(n+m))*sc_merged.data + (m/(n+m))*sc.data
        if !isnothing(sc_merged.M)
            @. sc_merged.M = (sc_merged.M + n*abs(μ - sc_merged.data)^2) + (sc.M + m*abs(μ - sc.data)^2)
        end
        sc_merged.data .= μ
        sc_merged.nsamples[1] += m
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
    dynamic_correlations(sys::System; dt, nω, ωmax, 
                           observables=nothing, correlations=nothing) 

Creates an empty `SampledCorrelations` object for calculating and storing
dynamical structure factor intensities ``𝒮(𝐪,ω)``. Call [`add_sample!`](@ref)
to accumulate data for the given configuration of a spin system. Internally,
this will run a dynamical trajectory and measure time correlations. The
``𝒮(𝐪,ω)`` data can be retrieved by calling
[`intensities_interpolated`](@ref). Alternatively,
[`instant_intensities_interpolated`](@ref) will integrate out ``ω`` to obtain
``𝒮(𝐪)``, optionally applying classical-to-quantum correction factors.

Three keywords are required to specify the dynamics used for the trajectory
calculation.

- `dt`: The time step used for calculating the trajectory from which dynamic
    spin-spin correlations are calculated. The trajectories are calculated with
    an [`ImplicitMidpoint`](@ref) integrator.
- `ωmax`: The maximum energy, ``ω``, that will be resolved. Note that allowed
    values of `ωmax` are constrained by the given `dt`, so Sunny will choose the
    smallest possible value that is no smaller than the specified `ωmax`.
- `nω`: The number of energy bins to calculated between 0 and `ωmax`.

Additional keyword options are the following:
- `observables`: Allows the user to specify custom observables. The
    `observables` must be given as a list of complex `N×N` matrices or
    `LinearMap`s. It's recommended to name each observable, for example:
    `observables = [:A => a_observable_matrix, :B => b_map, ...]`. By default,
    Sunny uses the 3 components of the dipole, `:Sx`, `:Sy` and `:Sz`.
- `correlations`: Specify which correlation functions are calculated, i.e. which
    matrix elements ``αβ`` of ``𝒮^{αβ}(q,ω)`` are calculated and stored.
    Specified with a vector of tuples. By default Sunny records all auto- and
    cross-correlations generated by all `observables`. To retain only the xx and
    xy correlations, one would set `correlations=[(:Sx,:Sx), (:Sx,:Sy)]` or
    `correlations=[(1,1),(1,2)]`.
"""
function instant_correlations(sys::System; kwargs...)
    dynamic_correlations(sys; dt=NaN, nω=NaN, ωmax=NaN, kwargs...)
end

function SampledCorrelations(sys::System{N}; dt, energies, measure=nothing, calculate_errors=false) where N

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
        @assert length(unique(energies[2:end] - energies[1:end-1])) == 1 "`energies` must be uniformly spaced"
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
    nsamples = Int64[0]

    # Make Structure factor and add an initial sample
    origin_crystal = isnothing(sys.origin) ? nothing : sys.origin.crystal
    sc = SampledCorrelations{N}(data, M, sys.crystal, origin_crystal, Δω, measure, measperiod, dt, nsamples,
                                samplebuf, corrbuf, space_fft!, time_fft!, corr_fft!, corr_ifft!)

    return sc
end



"""
    instant_correlations(sys::System; process_trajectory=:none, observables=nothing, correlations=nothing) 

Creates an empty `SampledCorrelations` object for calculating and storing
instantaneous structure factor intensities ``𝒮(𝐪)``. Call
[`add_sample!`](@ref) to accumulate data for the given configuration of a spin
system. Call [`instant_intensities_interpolated`](@ref) to retrieve averaged
``𝒮(𝐪)`` data.

_Important note_: When dealing with continuous (non-Ising) spins, consider
creating using [`dynamic_correlations`](@ref) instead of
`instant_correlations`. The former will provide full ``𝒮(𝐪,ω)`` data, from
which ``𝒮(𝐪)`` can be obtained by integrating out ``ω``. During this
integration step, Sunny can incorporate temperature- and ``ω``-dependent
classical-to-quantum correction factors to produce more accurate ``𝒮(𝐪)``
estimates. See [`instant_intensities_interpolated`](@ref) for more information.

The following optional keywords are available:

- `observables`: Allows the user to specify custom observables. The
    `observables` must be given as a list of complex `N×N` matrices or
    `LinearMap`s. It's recommended to name each observable, for example:
    `observables = [:A => a_observable_matrix, :B => b_map, ...]`. By default,
    Sunny uses the 3 components of the dipole, `:Sx`, `:Sy` and `:Sz`.
- `correlations`: Specify which correlation functions are calculated, i.e. which
    matrix elements ``αβ`` of ``𝒮^{αβ}(q,ω)`` are calculated and stored.
    Specified with a vector of tuples. By default Sunny records all auto- and
    cross-correlations generated by all `observables`. To retain only the xx and
    xy correlations, one would set `correlations=[(:Sx,:Sx), (:Sx,:Sy)]` or
    `correlations=[(1,1),(1,2)]`.
"""
