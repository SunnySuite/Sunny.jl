"""
    SampledCorrelations

Basic data type for storing sampled correlation data. A `SampleCorrelations` is
initialized by calling either [`dynamical_correlations`](@ref) or
[`instant_correlations`](@ref).
"""
struct SampledCorrelations{N}
    # 𝒮^{αβ}(q,ω) data and metadata
    data           :: Array{ComplexF64, 7}                 # Raw SF data for 1st BZ (numcorrelations × natoms × natoms × latsize × energy)
    M              :: Union{Nothing, Array{Float64, 7}}    # Running estimate of (nsamples - 1)*σ² (where σ² is the variance of intensities)
    crystal        :: Crystal                              # Crystal for interpretation of q indices in `data`
    origin_crystal :: Union{Nothing,Crystal}               # Original user-specified crystal (if different from above) -- needed for FormFactor accounting
    Δω             :: Float64                              # Energy step size (could make this a virtual property)  
    observables    :: ObservableInfo

    # Specs for sample generation and accumulation
    samplebuf    :: Array{ComplexF64, 6}   # New sample buffer
    fft!         :: FFTW.AbstractFFTs.Plan # Pre-planned FFT
    measperiod   :: Int                    # Steps to skip between saving observables (downsampling for dynamical calcs)
    apply_g      :: Bool                   # Whether to apply the g-factor
    dt           :: Float64                # Step size for trajectory integration 
    nsamples     :: Array{Int64, 1}        # Number of accumulated samples (array so mutable)
    processtraj! :: Function               # Function to perform post-processing on sample trajectories
end

function Base.show(io::IO, sc::SampledCorrelations{N}) where N
    modename = N == 0 ? "Dipole" : "SU($(N))"
    print(io,"SampledCorrelations{$modename}")
    print(io,all_observable_names(sc.observables))
end

function Base.show(io::IO, ::MIME"text/plain", sc::SampledCorrelations{N}) where N
    printstyled(io, "SampledCorrelations";bold=true, color=:underline)
    modename = N == 0 ? "Dipole" : "SU($(N))"
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
    print(io,"$(num_correlations(sc.observables)) correlations in $modename mode:\n")

    show(io,"text/plain",sc.observables)
end

Base.getproperty(sc::SampledCorrelations, sym::Symbol) = sym == :latsize ? size(sc.samplebuf)[2:4] : getfield(sc,sym)

function clone_correlations(sc::SampledCorrelations{N}) where N
    dims = size(sc.data)[2:4]
    nω = size(sc.data, 7)
    normalizationFactor = 1/(√(nω * prod(dims)))
    fft! = normalizationFactor * FFTW.plan_fft!(sc.samplebuf, (2,3,4,6)) # Avoid copies/deep copies of C-generated data structures
    M = isnothing(sc.M) ? nothing : copy(sc.M)
    return SampledCorrelations{N}(copy(sc.data), M, sc.crystal, sc.origin_crystal, sc.Δω,
        deepcopy(sc.observables), copy(sc.samplebuf), fft!, sc.measperiod, sc.apply_g, sc.dt,
        copy(sc.nsamples), sc.processtraj!)
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

"""
    dynamical_correlations(sys::System; dt, nω, ωmax, 
        process_trajectory=:none, observables=nothing, correlations=nothing) 

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
function dynamical_correlations(sys::System{N}; dt=nothing, Δt=nothing, nω, ωmax,
                                apply_g=true, observables=nothing, correlations=nothing,
                                calculate_errors=false, process_trajectory=:none) where N
    if !isnothing(Δt)
        @warn "`Δt` argument is deprecated! Use `dt` instead."
        dt = @something dt Δt
    end
    isnothing(dt) && error("`dt` parameter required")

    observables = parse_observables(N; observables, correlations, g = apply_g ? sys.gs : nothing)

    # Determine trajectory measurement parameters
    nω = Int64(nω)
    if nω != 1
        @assert π/dt > ωmax "Desired `ωmax` not possible with specified `dt`. Choose smaller `dt` value."
        measperiod = floor(Int, π/(dt * ωmax))
        nω = 2nω-1  # Ensure there are nω _non-negative_ energies.
        Δω = 2π / (dt*measperiod*nω)
    else
        measperiod = 1
        dt = Δω = NaN
    end

    # Set up trajectory processing function (e.g., symmetrize)
    processtraj! = if process_trajectory == :none 
        no_processing
    elseif process_trajectory == :symmetrize
        symmetrize!
    elseif process_trajectory == :subtract_mean
        subtract_mean!
    else
        error("Unknown argument for `process_trajectory`")
    end

    # Preallocation
    na = natoms(sys.crystal)
    samplebuf = zeros(ComplexF64, num_observables(observables), sys.latsize..., na, nω) 
    data = zeros(ComplexF64, num_correlations(observables), na, na, sys.latsize..., nω)
    M = calculate_errors ? zeros(Float64, size(data)...) : nothing

    # Specially Normalized FFT.
    # This is designed so that when it enters ifft(fft * fft) later (in squared fashion)
    # it will result in the 1/N factor needed to average over the
    # N-many independent estimates of the correlation.
    normalizationFactor = 1/(√(nω * prod(sys.latsize)))
    fft! = normalizationFactor * FFTW.plan_fft!(samplebuf, (2,3,4,6))

    # Other initialization
    nsamples = Int64[0]

    # Make Structure factor and add an initial sample
    origin_crystal = isnothing(sys.origin) ? nothing : sys.origin.crystal
    sc = SampledCorrelations{N}(data, M, sys.crystal, origin_crystal, Δω, observables,
                                samplebuf, fft!, measperiod, apply_g, dt, nsamples, processtraj!)

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
creating using [`dynamical_correlations`](@ref) instead of
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
function instant_correlations(sys::System; kwargs...)
    dynamical_correlations(sys; dt=NaN, nω=1, ωmax=NaN, kwargs...)
end
