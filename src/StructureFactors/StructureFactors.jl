struct StructureFactor{N}
    # 𝒮^{αβ}(q,ω) data and metadata
    data           :: Array{ComplexF64, 7}   # Raw SF data for 1st BZ (numcorrelations × natoms × natoms × latsize × energy)
    crystal        :: Crystal                # Crystal for interpretation of q indices in `data`
    origin_crystal :: Union{Nothing,Crystal} # Original user-specified crystal (if different from above)
    Δω             :: Float64                # Energy step size

    # Correlation info (αβ indices of 𝒮^{αβ}(q,ω))
    dipole_corrs :: Bool                                  # Whether using all correlations from dipoles 
    observables  :: Array{ComplexF64, 3}                  # Operators corresponding to observables
    idxinfo      :: SortedDict{CartesianIndex{2}, Int64}  # (α, β) to save from 𝒮^{αβ}(q, ω)

    # Specs for sample generation and accumulation
    samplebuf    :: Array{ComplexF64, 6}  # New sample buffer
    fft!         :: FFTW.AbstractFFTs.Plan # Pre-planned FFT
    copybuf      :: Array{ComplexF64, 4}  # Copy cache for accumulating samples
    measperiod   :: Int                   # Steps to skip between saving observables (downsampling for dynamical calcs)
    apply_g      :: Bool                  # Whether to apply the g-factor
    integrator   :: ImplicitMidpoint      # Integrator for dissipationless trajectories (will likely move to add_sample!)
    nsamples     :: Array{Int64, 1}       # Number of accumulated samples (array so mutable)
    processtraj! :: Function              # Function to perform post-processing on sample trajectories
end

function Base.show(io::IO, ::MIME"text/plain", sf::StructureFactor)
    modename = sf.dipole_corrs ? "Dipole correlations" : "Custom correlations"
    printstyled(io, "StructureFactor [$modename]\n"; bold=true, color=:underline)
end

"""
    merge!(sf::StructureFactor, others...)

Accumulate the samples in `others` (one or more `StructureFactors`) into `sf`.
"""
function merge!(sf::StructureFactor, others...)
    for sfnew in others
        nnew = sfnew.nsamples[1]
        ntotal = sf.nsamples[1] + nnew
        @. sf.data = sf.data + (sfnew.data - sf.data) * (nnew/ntotal)
        sf.nsamples[1] = ntotal
    end
end


"""
    StructureFactor

An object holding ``𝒮(𝐪,ω)`` or ``𝒮(𝐪)`` data. Construct a `StructureFactor`
using [`DynamicStructureFactor`](@ref) or [`InstantStructureFactor`](@ref),
respectively.
"""
function StructureFactor(sys::System{N}; Δt, nω, measperiod,
                            apply_g = true, observables = nothing, correlations = nothing,
                            process_trajectory = :none) where {N}

    # Set up correlation functions (which matrix elements αβ to save from 𝒮^{αβ})
    default_observables = false
    default_correlations = false
    if isnothing(observables)
        observables = zeros(ComplexF64, 0, 0, 3)  # observables are empty in this case
        default_observables = true
    else
        (N == 0) && error("Structure Factor Error: Cannot provide matrices for observables when using dipolar `System`")
    end
    nops = size(observables, 3)
    if isnothing(correlations)
        correlations = []
        for i in 1:nops, j in i:nops
            push!(correlations, (i, j))
        end
        default_correlations = true
    end
    dipole_corrs = default_observables && default_correlations

    # Construct look-up table for matrix elements
    count = 1
    pairs = []
    for αβ in correlations
        α, β = αβ
        α, β = α < β ? (α, β) : (β, α)  # Because SF is symmetric, only save diagonal and upper triangular
        push!(pairs, (α, β) => count)
        count += 1
    end
    pairs = map(i -> CartesianIndex(i.first) => i.second, pairs) # Convert to CartesianIndices
    idxinfo = SortedDict{CartesianIndex{2}, Int64}(pairs) # CartesianIndices sort to fastest order

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
    samplebuf = zeros(ComplexF64, nops, sys.latsize..., na, nω) 
    copybuf = zeros(ComplexF64, sys.latsize..., nω) 
    data = zeros(ComplexF64, ncorr, na, na, sys.latsize..., nω)

    # Normalize FFT according to physical convention
    normalizationFactor = 1/(nω * √(prod(sys.latsize)))
    fft! = normalizationFactor * FFTW.plan_fft!(samplebuf, (2,3,4,6))

    # Other initialization
    nsamples = Int64[0]
    integrator = ImplicitMidpoint(Δt)
    Δω = nω == 1 ? 0.0 : 2π / (Δt*measperiod*nω)
    origin_crystal = !isnothing(sys.origin) ? sys.origin.crystal : nothing

    # Make Structure factor and add an initial sample
    sf = StructureFactor{N}(data, sys.crystal, origin_crystal, Δω, dipole_corrs,
                            observables, idxinfo, samplebuf, fft!, copybuf, measperiod, apply_g, integrator,
                            nsamples, processtraj!)
    add_sample!(sf, sys; processtraj!)

    return sf
end


"""
    DynamicStructureFactor(sys::System; Δt, nω, ωmax, 
        process_trajectory=:none, observables=nothing, correlations=nothing) 

Creates a `StructureFactor` for calculating and storing ``𝒮(𝐪,ω)`` data. This
information will be obtained by running dynamical spin simulations on
equilibrium snapshots, and measuring pair-correlations. The ``𝒮(𝐪,ω)`` data
can be retrieved by calling [`intensities`](@ref). Alternatively,
[`instant_intensities`](@ref) will integrate out ``ω`` to obtain ``𝒮(𝐪)``,
optionally applying classical-to-quantum correction factors.
        
Prior to calling `DynamicStructureFactor`, ensure that `sys` represents a good
equilibrium sample. Additional sample data may be accumulated by calling
[`add_sample!`](@ref)`(sf, sys)` with newly equilibrated `sys` configurations.

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
- `observables`: Enables an advanced feature for SU(_N_) mode, allowing the user
    to specify custom observables other than the three components of the dipole.
    To use this features, `observables` must be given an `N×N×numops` array,
    where the final index is used to retrieve each `N×N` operator.
- `correlations`: Specify which correlation functions are calculated, i.e. which
    matrix elements ``αβ`` of ``𝒮^{αβ}(q,ω)`` are calculated and stored.
    Specified with a vector of tuples. By default Sunny records all auto- and
    cross-correlations generated by the x, y, and z dipolar components (1, 2,
    and 3 respectively). To retain only the xx and xy correlations, one would
    set `correlations=[(1,1), (1,2)]`. If custom observables (`observables`) are
    given, the indices are ordered in the same manner as the final index of
    `ops`.
"""
function DynamicStructureFactor(sys::System; Δt, nω, ωmax, kwargs...) 
    nω = Int64(nω)
    @assert π/Δt > ωmax "Desired `ωmax` not possible with specified `Δt`. Choose smaller `Δt` value."
    measperiod = floor(Int, π/(Δt * ωmax))
    nω = 2nω-1  # Ensure there are nω _non-negative_ energies
    StructureFactor(sys; Δt, nω, measperiod, kwargs...)
end


"""
    InstantStructureFactor(sys::System; process_trajectory=:none,
                            observables=nothing, correlations=nothing) 

Creates a `StructureFactor` object for calculating and storing instantaneous
structure factor intensities ``𝒮(𝐪)``. This data will be calculated from the
spin-spin correlations of equilibrium snapshots, absent any dynamical
information. ``𝒮(𝐪)`` data can be retrieved by calling
[`instant_intensities`](@ref).

_Important note_: When dealing with continuous (non-Ising) spins, consider
creating a full [`DynamicStructureFactor`](@ref) object instead of an
`InstantStructureFactor`. The former will provide full ``𝒮(𝐪,ω)`` data, from
which ``𝒮(𝐪)`` can be obtained by integrating out ``ω``. During this
integration step, Sunny can incorporate temperature- and ``ω``-dependent
classical-to-quantum correction factors to produce more accurate ``𝒮(𝐪)``
estimates. See [`instant_intensities`](@ref) for more information.

Prior to calling `InstantStructureFactor`, ensure that `sys` represents a good
equilibrium sample. Additional sample data may be accumulated by calling
[`add_sample!`](@ref)`(sf, sys)` with newly equilibrated `sys` configurations.

The following optional keywords are available:

- `process_trajectory`: Specifies a function that will be applied to the sample
    trajectory before correlation analysis. Current options are `:none` and
    `:symmetrize`. The latter will symmetrize the trajectory in time, which can
    be useful for removing Fourier artifacts that arise when calculating the
    correlations.
- `observables`: Enables an advanced feature for SU(_N_) mode, allowing the user
    to specify custom observables other than the three components of the dipole.
    To use this features, `observables` must be given an `N×N×numops` array,
    where the final index is used to retrieve each `N×N` operator.
- `correlations`: Specify which correlation functions are calculated, i.e. which
    matrix elements ``αβ`` of ``𝒮^{αβ}(q,ω)`` are calculated and stored.
    Specified with a vector of tuples. By default Sunny records all auto- and
    cross-correlations generated by the x, y, and z dipolar components (1, 2,
    and 3 respectively). To retain only the xx and xy correlations, one would
    set `correlations=[(1,1), (1,2)]`. If custom observables (`observables`) are
    given, the indices are ordered in the same manner as the final index of
    `observables`.
"""
function InstantStructureFactor(sys::System; kwargs...)
    StructureFactor(sys; Δt=0.1, nω=1, measperiod=1, kwargs...)
end
