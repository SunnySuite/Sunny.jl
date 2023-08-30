"""
    SampledCorrelations

Basic data type for storing sampled correlation data. A `SampleCorrelations` is
initialized by calling either [`dynamical_correlations`](@ref) or
[`instant_correlations`](@ref).
"""
struct SampledCorrelations{N}
    # 𝒮^{αβ}(q,ω) data and metadata
    data           :: Array{ComplexF64, 7}                 # Raw SF data for 1st BZ (numcorrelations × natoms × natoms × latsize × energy)
    variance       :: Union{Nothing, Array{Float64, 7}}    # Running variance calculation for Welford's algorithm 
    crystal        :: Crystal                              # Crystal for interpretation of q indices in `data`
    origin_crystal :: Union{Nothing,Crystal}               # Original user-specified crystal (if different from above) -- needed for FormFactor accounting
    Δω             :: Float64                              # Energy step size (could make this a virtual property)  

    # Correlation info (αβ indices of 𝒮^{αβ}(q,ω))
    observables    :: Vector{LinearMap}  # Operators corresponding to observables
    observable_ixs :: Dict{Symbol,Int64} # User-defined observable names
    correlations   :: SortedDict{CartesianIndex{2}, Int64}  # (α, β) to save from 𝒮^{αβ}(q, ω)

    # Specs for sample generation and accumulation
    samplebuf    :: Array{ComplexF64, 6}   # New sample buffer
    fft!         :: FFTW.AbstractFFTs.Plan # Pre-planned FFT
    copybuf      :: Array{ComplexF64, 4}   # Copy cache for accumulating samples
    measperiod   :: Int                    # Steps to skip between saving observables (downsampling for dynamical calcs)
    apply_g      :: Bool                   # Whether to apply the g-factor
    Δt           :: Float64                # Step size for trajectory integration 
    nsamples     :: Array{Int64, 1}        # Number of accumulated samples (array so mutable)
    processtraj! :: Function               # Function to perform post-processing on sample trajectories
end

function Base.show(io::IO, sc::SampledCorrelations{N}) where N
    modename = N == 0 ? "Dipole" : "SU($(N))"
    print(io,"SampledCorrelations{$modename}")
    print(io,all_observable_names(sc))
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
    print(io,"$(size(sc.data)[1]) correlations in $modename mode:\n")

    # Reverse the dictionary
    observable_names = Dict(value => key for (key, value) in sc.observable_ixs)

    for i = 1:length(sc.observables)
        print(io,i == 1 ? "╔ " : i == length(sc.observables) ? "╚ " : "║ ")
        for j = 1:length(sc.observables)
            if i > j
                print(io,"⋅ ")
            elseif haskey(sc.correlations,CartesianIndex(i,j))
                print(io,"⬤ ")
            else
                print(io,"• ")
            end
        end
        print(io,observable_names[i])
        println(io)
    end
    printstyled(io,"")
end

Base.getproperty(sc::SampledCorrelations, sym::Symbol) = sym == :latsize ? size(sc.samplebuf)[2:4] : getfield(sc,sym)

"""
    merge!(sc::SampledCorrelations, others...)

Accumulate the samples in `others` (one or more `SampledCorrelations`) into `sc`.
"""
function merge!(sc::SampledCorrelations, others...)
    for scnew in others
        nnew = scnew.nsamples[1]
        ntotal = sc.nsamples[1] + nnew
        @. sc.data = sc.data + (scnew.data - sc.data) * (nnew/ntotal)
        sc.nsamples[1] = ntotal
    end
end

# Finds the linear index according to sc.correlations of each correlation in corrs, where
# corrs is of the form [(:A,:B),(:B,:C),...] where :A,:B,:C are observable names.
function lookup_correlations(sc::SampledCorrelations,corrs; err_msg = αβ -> "Missing correlation $(αβ)")
    indices = Vector{Int64}(undef,length(corrs))
    for (i,(α,β)) in enumerate(corrs)
        αi = sc.observable_ixs[α]
        βi = sc.observable_ixs[β]
        # Make sure we're looking up the correlation with its properly sorted name
        αi,βi = minmax(αi,βi)
        idx = CartesianIndex(αi,βi)

        # Get index or fail with an error
        indices[i] = get!(() -> error(err_msg(αβ)),sc.correlations,idx)
    end
    indices
end

function all_observable_names(sc::SampledCorrelations)
    observable_names = Dict(value => key for (key, value) in sc.observable_ixs)
    [observable_names[i] for i in 1:length(observable_names)]
end

"""
    dynamical_correlations(sys::System; Δt, nω, ωmax, 
        process_trajectory=:none, observables=nothing, correlations=nothing) 

Creates a `SampledCorrelations` for calculating and storing ``𝒮(𝐪,ω)`` data.
This information will be obtained by running dynamical spin simulations on
equilibrium snapshots and measuring pair-correlations. The ``𝒮(𝐪,ω)`` data can
be retrieved by calling [`intensities_interpolated`](@ref). Alternatively,
[`instant_intensities_interpolated`](@ref) will integrate out ``ω`` to obtain
``𝒮(𝐪)``, optionally applying classical-to-quantum correction factors.
        
The `SampleCorrelations` that is returned will contain no correlation data.
Samples are generated and accumulated by calling [`add_sample!`](@ref)`(sc,
sys)` where `sc` is a `SampleCorrelations` and `sys` is an appropriately
equilibrated `System`. Note that the `sys` should be thermalized before each
call of `add_sample!` such that the spin configuration in the system represents
a new (fully decorrelated) sample.

Three keywords are required to specify the dynamics used for the trajectory
calculation.

- `Δt`: The time step used for calculating the trajectory from which dynamic
    spin-spin correlations are calculated. The trajectories are calculated with
    an [`ImplicitMidpoint`](@ref) integrator.
- `ωmax`: The maximum energy, ``ω``, that will be resolved.
- `nω`: The number of energy bins to calculated between 0 and `ωmax`.

Additional keyword options are the following:
- `process_trajectory`: Specifies a function that will be applied to the sample
    trajectory before correlation analysis. Current options are `:none` and
    `:symmetrize`. The latter will symmetrize the trajectory in time, which can
    be useful for removing Fourier artifacts that arise when calculating the
    correlations.
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
function dynamical_correlations(sys::System{N}; Δt, nω, ωmax,
                                apply_g = true, observables = nothing, correlations = nothing,
                                calculate_errors = false, process_trajectory = :none) where {N}

    # Set up correlation functions (which matrix elements αβ to save from 𝒮^{αβ})
    if isnothing(observables)
        # Default observables are spin x,y,z
        # projections (SU(N) mode) or components (dipole mode)
        observable_ixs = Dict(:Sx => 1,:Sy => 2,:Sz => 3)
        if N == 0
            dipole_component(i) = FunctionMap{Float64}(s -> s[i],1,3)
            observables = dipole_component.([1,2,3])
        else
            # SQTODO: Make this use the more optimized expected_spin function
            # Doing this will also, by necessity, allow users to make the same
            # type of optimization for their vector-valued observables.
            observables = LinearMap{ComplexF64}.(spin_matrices(;N))
        end
    else
        # If it was given as a list, preserve the user's preferred
        # ordering of observables
        if observables isa AbstractVector
            # If they are pairs (:A => [...]), use the names
            # and otherwise use alphabetical names
            if !isempty(observables) && observables[1] isa Pair
                observables = OrderedDict(observables)
            else
                dict = OrderedDict{Symbol,LinearMap}()
                for i = 1:length(observables)
                    dict[Symbol('A' + i - 1)] = observables[i]
                end
                observables = dict
            end
        end

        # If observables were provided as (:name => matrix) pairs,
        # reformat them to (:name => idx) and matrices[idx]
        observable_ixs = Dict{Symbol,Int64}()
        matrices = Vector{LinearMap}(undef,length(observables))
        for (i,name) in enumerate(keys(observables))
            next_available_ix = length(observable_ixs) + 1
            if haskey(observable_ixs,name)
                error("Repeated observable name $name not allowed.")
            end
            observable_ixs[name] = next_available_ix

            # Convert dense matrices to LinearMap
            if observables[name] isa Matrix
                matrices[i] = LinearMap(observables[name])
            else
                matrices[i] = observables[name]
            end
        end
        observables = matrices
    end

    # By default, include all correlations
    if isnothing(correlations)
        correlations = []
        for oi in keys(observable_ixs), oj in keys(observable_ixs)
            push!(correlations, (oi, oj))
        end
    elseif correlations isa AbstractVector{Tuple{Int64,Int64}}
        # If the user used numeric indices to describe the correlations,
        # we need to convert it to the names, so need to temporarily reverse
        # the dictionary.
        observable_names = Dict(value => key for (key, value) in observable_ixs)
        correlations = [(observable_names[i],observable_names[j]) for (i,j) in correlations]
    end

    # Construct look-up table for correlation matrix elements
    idxinfo = SortedDict{CartesianIndex{2},Int64}() # CartesianIndex's sort to fastest order
    for (α,β) in correlations
        αi = observable_ixs[α]
        βi = observable_ixs[β]
        # Because correlation matrix is symmetric, only save diagonal and upper triangular
        # by ensuring that all pairs are in sorted order
        αi,βi = minmax(αi,βi)
        idx = CartesianIndex(αi,βi)

        # Add this correlation to the list if it's not already listed
        get!(() -> length(idxinfo) + 1,idxinfo,idx)
    end
    correlations = idxinfo

    # Determine trajectory measurement parameters
    nω = Int64(nω)
    if nω != 1
        @assert π/Δt > ωmax "Desired `ωmax` not possible with specified `Δt`. Choose smaller `Δt` value."
        measperiod = floor(Int, π/(Δt * ωmax))
        nω = 2nω-1  # Ensure there are nω _non-negative_ energies
        Δω = 2π / (Δt*measperiod*nω)
    else
        measperiod = 1
        Δt = Δω = NaN
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
    ncorr = length(correlations)
    samplebuf = zeros(ComplexF64, length(observables), sys.latsize..., na, nω) 
    copybuf = zeros(ComplexF64, sys.latsize..., nω) 
    data = zeros(ComplexF64, ncorr, na, na, sys.latsize..., nω)
    variance = calculate_errors ? zeros(Float64, size(data)...) : nothing

    # Normalize FFT according to physical convention
    normalizationFactor = 1/(nω * √(prod(sys.latsize)))
    fft! = normalizationFactor * FFTW.plan_fft!(samplebuf, (2,3,4,6))

    # Other initialization
    nsamples = Int64[0]

    # Make Structure factor and add an initial sample
    origin_crystal = isnothing(sys.origin) ? nothing : sys.origin.crystal
    sc = SampledCorrelations{N}(data, variance, sys.crystal, origin_crystal, Δω, observables, observable_ixs, correlations,
                                samplebuf, fft!, copybuf, measperiod, apply_g, Δt, nsamples, processtraj!)

    return sc
end



"""
    instant_correlations(sys::System; process_trajectory=:none, observables=nothing, correlations=nothing) 

Creates a `SampledCorrelations` object for calculating and storing instantaneous
structure factor intensities ``𝒮(𝐪)``. This data will be calculated from the
spin-spin correlations of equilibrium snapshots, absent any dynamical
information. ``𝒮(𝐪)`` data can be retrieved by calling
[`instant_intensities_interpolated`](@ref).

_Important note_: When dealing with continuous (non-Ising) spins, consider
creating using [`dynamical_correlations`](@ref) instead of 
`instant_correlations`. The former will provide full ``𝒮(𝐪,ω)`` data, from
which ``𝒮(𝐪)`` can be obtained by integrating out ``ω``. During this
integration step, Sunny can incorporate temperature- and ``ω``-dependent
classical-to-quantum correction factors to produce more accurate ``𝒮(𝐪)``
estimates. See [`instant_intensities_interpolated`](@ref) for more information.

Prior to calling `instant_correlations`, ensure that `sys` represents a good
equilibrium sample. Additional sample data may be accumulated by calling
[`add_sample!`](@ref)`(sc, sys)` with newly equilibrated `sys` configurations.

The following optional keywords are available:

- `process_trajectory`: Specifies a function that will be applied to the sample
    trajectory before correlation analysis. Current options are `:none` and
    `:symmetrize`. The latter will symmetrize the trajectory in time, which can
    be useful for removing Fourier artifacts that arise when calculating the
    correlations.
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
    dynamical_correlations(sys; Δt=NaN, nω=1, ωmax=NaN, kwargs...)
end
