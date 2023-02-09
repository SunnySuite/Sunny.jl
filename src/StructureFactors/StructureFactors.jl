struct StructureFactor{N, NumCorr, NBasis}
    # 𝒮^{αβ}(q,ω) data and metadata
    data           :: Array{ComplexF64, 7}   # Raw SF data for 1st BZ (numcorrelations x nbasis x nbasis x latsize x energy)
    crystal        :: Crystal                # Crystal for interpretation of q indices in `data`
    origin_crystal :: Union{Nothing,Crystal} # Original user-specified crystal (if different from above)
    Δω             :: Float64                # Energy step size

    # Correlation info (αβ indices of 𝒮^{αβ}(q,ω))
    dipole_corrs :: Bool                                  # Whether using all correlations from dipoles 
    observables  :: Array{ComplexF64, 3}                  # Operators corresponding to observables
    idxinfo      :: SortedDict{CartesianIndex{2}, Int64}  # (α, β) to save from 𝒮^{αβ}(q, ω)

    # Specs for sample generation and accumulation
    samplebuf    :: Array{ComplexF64, 6}  # New sample buffer
    measperiod   :: Int                   # Steps to skip between saving observables (downsampling for dynamical calcs)
    apply_g      :: Bool                  # Whether to apply the g-factor
    integrator   :: ImplicitMidpoint      # Integrator for dissipationless trajectories (will likely move to add_sample!)
    nsamples     :: Array{Int64, 1}       # Number of accumulated samples (array so mutable)
    processtraj! :: Function              # Function to perform post-processing on sample trajectories
end

"""
    StructureFactor(sys::System; Δt, nω, measperiod, apply_g = true, observables = nothing,
                        correlations = nothing, process_trajectory = :none)

`StructureFactor` is the basic type for calculating ``𝒮(q,ω)`` or ``𝒮(q)``
data, storing the results, and retrieving intensity information. 

Instead of creating `StructureFactor` directly, one should call either
either [`DynamicStructureFactor`](@ref) or [`InstantStructureFactor`](@ref).

Data may be retrieved from a `StructureFactor` by calling [`intensities`](@ref) 
or [`instant_intensities`](@ref). 
"""
function StructureFactor(sys::System{N}; Δt, nω, measperiod,
                            apply_g = true, observables = nothing, correlations = nothing,
                            process_trajectory = :none) where N

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
    nb = nbasis(sys.crystal)
    ncorr = length(pairs)
    samplebuf = zeros(ComplexF64, nops, sys.latsize..., nb, nω) 
    data = zeros(ComplexF64, length(correlations), nb, nb, sys.latsize..., nω)

    # Other initialization
    nsamples = Int64[0]
    integrator = ImplicitMidpoint(Δt)
    Δω = nω == 1 ? 0.0 : 2π / (Δt*measperiod*nω)
    origin_crystal = !isnothing(sys.origin) ? sys.origin.crystal : nothing

    # Make Structure factor and add an initial sample
    sf = StructureFactor{N, ncorr, nb}(data, sys.crystal, origin_crystal, Δω, dipole_corrs,
                                       observables, idxinfo, samplebuf, measperiod, apply_g, integrator,
                                       nsamples, processtraj!)
    add_sample!(sf, sys; processtraj!)

    return sf
end


"""
    DynamicStructureFactor(sys::System; Δt, nω, ωmax, 
        apply_g=true, process_trajectory=:none, observables=nothing, correlations=nothing) 

Creates a `StructureFactor` for calculating and storing ``𝒮(q,ω)`` data. When
calculating a dynamic structure factor from classical dynamics, it is necessary
to generate spin trajetories which are used to calculate correlations. The
initial conditions for these trajectories must be sample spin configurations
drawn from the equilibrium distribution at the desired temperature. One such
trajectory is calculated immediately when initializing a
`DynamicStructureFactor`, so the spins in the `sys` must represent a good sample
before calling this function. Additional sample trajectories are created and
accumulated into the `DynamicStructureFactor` by calling
[`add_sample!`](@ref)`(sf, sys)`. The spins in the `sys` should be set to new
sample configurations before each call to `add_sample!`. This can be achieved,
for example, with the [`Langevin`](@ref) dynamics.

Three keywords are required to specify the dynamics used for the trajectory
calculation.

- `Δt`: The time step used for calculating the trajectory from which dynamic
    spin-spin correlations are calculated. The trajectories are calculated with
    an [`ImplicitMidpoint`](@ref) integrator.
- `ωmax`: The maximum energy, ``ω``, that will be resolved.
- `nω`: The number of energy bins to calculated between 0 and `ωmax`.

Additional keyword options are the following:
- `apply_g`: Determines whether to apply the g-factor when calculating
    trajectories.
- `process_trajectory`: Specifies a function that will be applied to the sample
    trajectory before correlation analysis. Current options are `:none` and
    `:symmetrize`. The latter will symmetrize the trajectory in time, which can
    be useful for removing Fourier artifacts that arise when calculating the
    correlations.
- `observables`: Enables an advanced feature for SU(_N_) mode, allowing the user to
    specify custom observables other than the three components of the dipole. To
    use this features, `observables` must be given an `N×N×numops` array, where the
    final index is used to retrieve each `N×N` operator.
- `correlations`: Specify which correlation functions are calculated, i.e. which
    matrix elements ``αβ`` of ``𝒮^{αβ}(q,ω)`` are calculated and stored.
    Specified with a vector of tuples. By default Sunny records all auto- and
    cross-correlations generated by the x, y, and z dipolar components (1, 2,
    and 3 respectively). To retain only the xx and xy correlations, one would
    set `correlations=[(1,1), (1,2)]`. If custom observables (`observables`) are given,
    the indices are ordered in the same manner as the final index of `ops`.

The ``𝒮(q,ω)`` data can be retrieved by calling [`intensities`](@ref) or
[`instant_intensities`](@ref). 
"""
function DynamicStructureFactor(sys::System; Δt, nω, ωmax, kwargs...) 
    nω = Int64(nω)
    @assert π/Δt > ωmax "Desired `ωmax` not possible with specified `Δt`. Choose smaller `Δt` value."
    measperiod = floor(Int, π/(Δt * ωmax))
    nω = 2nω-1  # Ensure there are nω _non-negative_ energies
    StructureFactor(sys; Δt, nω, measperiod, kwargs...)
end


"""
    InstantStructureFactor(sys::System; apply_g=true, process_trajectory=:none,
                            observables=nothing, correlations=nothing) 

Creates a `StructureFactor` for calculating and storing ``𝒮(q)`` data, i.e.,
spin-spin correlation data calculated at single time steps. An initial sample is
generated from the spins in `sys` when calling `InstantStructureFactor`, so the
spins should represent a good equilibrium sample before this function is called.
Additional samples may be generated by calling [`add_sample!`](@ref)`(sf, sys)`.
The spins in the `sys` should be resampled before each call to `add_sample!`.

The the following optional keywords are available:

- `apply_g`: Determines whether to apply the g-factor when calculating
    trajectories.
- `process_trajectory`: Specifies a function that will be applied to the sample
    trajectory before correlation analysis. Current options are `:none` and
    `:symmetrize`. The latter will symmetrize the trajectory in time, which can
    be useful for removing Fourier artifacts that arise when calculating the
    correlations.
- `observables`: Enables an advanced feature for SU(_N_) mode, allowing the user to
    specify custom observables other than the three components of the dipole. To
    use this features, `observables` must be given an `N×N×numops` array, where the
    final index is used to retrieve each `N×N` operator.
- `correlations`: Specify which correlation functions are calculated, i.e. which
    matrix elements ``αβ`` of ``𝒮^{αβ}(q,ω)`` are calculated and stored.
    Specified with a vector of tuples. By default Sunny records all auto- and
    cross-correlations generated by the x, y, and z dipolar components (1, 2,
    and 3 respectively). To retain only the xx and xy correlations, one would
    set `correlations=[(1,1), (1,2)]`. If custom observables (`observables`) are given,
    the indices are ordered in the same manner as the final index of `observables`.

``𝒮(𝐪)`` data can be retrieved by calling [`instant_intensities`](@ref).

NOTE: It is often advisable to generate a instantaneous structure factor,
``𝒮(𝐪)``, from a dynamic structure factor, ``𝒮(𝐪,ω)``, by integrating out
``ω``, rather than calculating ``𝒮(𝐪)`` directly from spin-spin correlations
at single instances of time. This approach makes it possible to apply a
temperature- and ``ω``-dependent classical-to-quantum intensity rescaling to the
results. This can be done in Sunny by calculating a `DynamicStructureFactor` and
retrieving ``𝒮(𝐪)`` data with [`instant_intensities`](@ref), taking care to set
the `kT` keyword to the appropriate value. `instant_intensities` will then
integrate the ``ω`` information out after applying intensity corrections.
"""
function InstantStructureFactor(sys::System; kwargs...)
    StructureFactor(sys; Δt=0.1, nω=1, measperiod=1, kwargs...)
end