################################################################################
# Basic functions for retrieving 𝒮(q, ω) values
################################################################################

# Function for getting a single 𝒮(q, ω) intensity -- primarily internal
function calc_intensity(sf::StructureFactor{N, NumCorr}, k, latidx, ω, iω, contractor, temp, ffdata) where {N, NumCorr}
    (; crystal, data) = sf.sfdata
    elems = phase_averaged_elements(view(data,:,:,:,latidx,iω), k, crystal, ffdata, Val(NumCorr))
    intensity = contract(elems, k, contractor)
    if !isnothing(temp)
        intensity *= classical_to_quantum(ω, temp)
    end
    return intensity
end

# Note that requests for intensities often come in lists of nearby q values. Since the data
# is inherently discretized, this often results in repeated calls for values at the same 
# discrete points. Since basis reduction is done for each of this calls, this results in 
# a large amount of repeated calculation. This function analyzes repetitions in advance
# and prunes out repetitions. 
# This is ugly, but the speedup when tested on a few simple, realistic examples was 3-5x.
function pruned_stencil_info(sf::StructureFactor, qs, interp::InterpolationScheme{N}) where N
    # Count the number of contiguous regions with unchanging values.
    # If all values are unique, returns the length of q_info.
    # Note comparison is on m values rather than index values and the m values are the first
    # element of the a tuple, that is, we're checking x[1] == y[1] in the map.
    m_info = map(q -> stencil_points(sf, q, interp), qs)
    numregions = sum(map((x,y) -> x[1] == y[1] ? 0 : 1, m_info[1:end-1], m_info[2:end])) + 1
    
    # Remove repeated stencil points and count number of instances of each
    ms_ref, idcs_ref = stencil_points(sf, qs[1], interp)
    ms_all  = fill(ntuple(x->zero(Vec3), N), numregions)
    ms_all[1] = ms_ref 
    idcs_all = fill(ntuple(x->CartesianIndex((-1,-1,-1)), N), numregions)
    idcs_all[1] = idcs_ref 
    counts = zeros(Int64, numregions)
    c = counts[1] = 1
    for q in qs[2:end] 
        ms, idcs = stencil_points(sf, q, interp)
        if ms != ms_ref
            ms_ref = ms 
            c += 1
            ms_all[c] =  ms
            idcs_all[c] = idcs 
        end
        counts[c] += 1
    end
    @assert sum(counts) == length(m_info)

    # Calculate corresponding q (RLU) and k (global) vectors
    (; crystal, latsize) = sf.sftraj.sys
    recip_vecs = 2π*inv(crystal.lat_vecs)

    qs_all = map(ms_all) do ms
       map(m -> m ./ latsize, ms) 
    end

    ks_all = map(qs_all) do qs
        map(q -> recip_vecs * q, qs)
    end
    
    return (; qs_all, ks_all, idcs_all, counts)
end


"""
    get_intensities(sf::StructureFactor, qs::Array;
                       contraction = :trace, interpolation = nothing,
                       kT = nothing, newbasis = nothing, 
                       formfactors = nothing, negative_energies = false)

The basic function for retrieving ``𝒮(q,ω)`` information from a
`StructureFactor`. Takes an array of wave vectors of any dimension, `qs`, and
returns an array of intensities of the same dimension plus an added final index
for energy values. The energy values corresponding the final index can be
retrieved by calling [`ωvals`](@ref). The wave vectors should be specified as
3-vectors or 3-tuples in terms of reciprocal lattice units (RLUs).  

- `contraction`: Determines the operation performed on the matrix elements α and
    β of ``𝒮^{αβ}(q,ω)`` when requesting an intensity. By default, Sunny
    returns the trace: ``∑_α 𝒮^{αα}(q,ω)``. Polarization correction can be
    achieved by setting `contraction=:depolarize`. To retrieve a single matrix
    element, pass `contraction` a tuple of `Int`s. For example, to retrieve the
    `xy` correlation, set `contraction=(1,2)`.
- `interpolation`: Since ``𝒮(q,ω)`` is calculated on a finite lattice, data is
    only available at discrete wave vectors. By default, Sunny will round a
    requested `q` to the nearest available wave vector. Linear interpolation can
    be applied by setting `interpolation=:linear`.
- `kT`: If a temperature is provided, the intensities will be rescaled by a
    temperature- and ω-dependent classical-to-quantum factor. `kT` should be
    specified when making comparisons with spin wave calculations or
    experimental data.
- `newbasis`: If the `qs` are given in terms of a basis other than the
    reciprocal lattice vectors, `newbasis` should be given a 3x3 matrix, the
    columns of which determine the new basis by specifying linear combinations
    of the reciprocal vectors. 
- `formfactors`: To apply form factor corrections, provide this keyword with a
    vector of `FormFactor`s, one for each unique site in the unit cell. Sunny
    will symmetry propagate the results to all equivalent sites.
- `negative_energies`: If set to `true`, Sunny will return the periodic
    extension of the energy axis. Most users will not want this.
"""
function get_intensities(sf::StructureFactor, qs::Array;
    contraction = :trace,
    interpolation = nothing,
    kT = nothing,
    newbasis = nothing,
    formfactors = nothing,
    negative_energies = false,
)

    # Set up interpolation scheme
    interp = if isnothing(interpolation) 
        NoInterp()
    elseif interpolation == :linear
        LinearInterp()
    end

    # Set up element contraction
    contract = if contraction == :trace
        Trace(sf)
    elseif contraction == :depolarize
        Depolarize(sf)
    elseif typeof(contraction) <: Tuple{Int, Int}
        Element(sf, contraction)
    end

    # Apply basis transformation (if newbasis provided)
    q_targets = Vec3.(qs) 
    if !isnothing(newbasis)
        newbasis = Mat3(newbasis)
        q_targets = map(q -> newbasis*q, q_targets)
    end

    # Propagate form factor information (if any)
    cryst = sf.sfdata.crystal
    if isnothing(formfactors)
        formfactors = [FormFactor{EMPTY_FF}(; atom) for atom in unique(cryst.classes)]
    end
    ffdata = propagate_form_factors(cryst, formfactors)

    # Precompute index information and preallocate
    ωs = negative_energies ? ωvals_all(sf) : ωvals(sf)
    nω = length(ωs) 
    intensities = zeros(contract, size(q_targets)..., nω)
    stencil_info = pruned_stencil_info(sf, q_targets, interp) 
    
    # Call type stable version of the function
    get_intensities!(intensities, sf, q_targets, ωs, interp, contract, kT, ffdata, stencil_info) #ddtodo: Track down allocations

    # ddtodo: See if worth it to apply classical-to-quantum rescaling here instead of inside loop

    return intensities
end

Base.zeros(::Contraction{T}, dims...) where T = zeros(T, dims...)

# Type stable version
function get_intensities!(intensities, sf::StructureFactor, q_targets::Array, ωs, interp::InterpolationScheme, contraction::Contraction{T}, temp, ffdata, stencil_info) where {T}
    li_intensities = LinearIndices(intensities)
    ci_qs = CartesianIndices(q_targets)
    (; qs_all, ks_all, idcs_all, counts) = stencil_info 
    for (iω, ω) in enumerate(ωs)
        iq = 0
        for (qs, ks, idcs, numrepeats) in zip(qs_all, ks_all, idcs_all, counts)
            local_intensities = stencil_intensities(sf, ks, idcs, ω, iω, interp, contraction, temp, ffdata) 
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(ci_qs[iq], iω)]
                intensities[idx] = interpolated_intensity(sf, q_targets[iq], qs, local_intensities, interp) 
            end
        end
    end
    return intensities
end


function get_intensities(sf::StructureFactor, q::NTuple{3,Number}; kwargs...)
    return get_intensities(sf, [Vec3(q)]; kwargs...)[1,:]
end

function get_intensities(sf::StructureFactor, q::AbstractArray{Number,1}; kwargs...)
    if length(q) != 3
        error("Q point should have three components.")
    end
    return get_intensities(sf, [Vec3(q)]; kwargs...)[1,:]
end


"""
    get_static_intensities(sf::StructureFactor, qs::Array; kwargs...)

Return the static structure factor intensities at wave vectors `qs`. The
functionality is very similar to [`get_intensities`](@ref), except the returned
array has dimensions identical to `qs`. The energy axis has been summed out.
"""
function get_static_intensities(sf::StructureFactor, qs::Array; kwargs...)
    datadims = size(qs)
    ndims = length(datadims)
    intensities = get_intensities(sf, qs; kwargs...)
    static_intensities = sum(intensities, dims=(ndims+1,))
    return reshape(static_intensities, datadims)
end

function get_static_intensities(sf::StructureFactor, q::NTuple{3, T}; kwargs...) where T <: Number
    intensities = get_intensity(sf, [Vec3(q)]; kwargs...)
    return sum(intensities)
end

function get_static_intensities(sf::StructureFactor, q::AbstractArray{T,1}; kwargs...) where T <: Number
    if length(q) != 3
        error("Q point should have three components.")
    end
    intensities = get_intensity(sf, [Vec3(q)]; kwargs...)
    return sum(intensities)
end



"""
    intensity_grid(sf::StructureFactor;
                       bzsize=(1,1,1), negative_energies = false, index_labels = false, kwargs...)

Returns intensities at discrete wave vectors for which there is exact
information. Shares all keywords with [`get_intensities`](@ref), and provides
two additional options:

- `bzsize`: Specifies the number of Brillouin zones to return, given as a
  3-tuple of integers.
- `index_labels`: If set to `true`, will return axis label information for the
    data, which may be upacked as: `(; intensities, qpoints, ωs)`.
"""
function intensity_grid(sf::StructureFactor;
                            bzsize=(1,1,1), negative_energies = false, index_labels = false, kwargs...)
    qpoints = qgrid(sf; bzsize)
    intensities = get_intensities(sf, qpoints; negative_energies, kwargs...)

    if index_labels
        ωs =  negative_energies ? ωvals_all(sf) : ωvals(sf)
        return (; intensities, qpoints, ωs)
    end
    return intensities
end


function path_points(points::Vector, density)
    legs = []
    for i ∈ 1:length(points)-1
        leg = []
        p1, p2 = points[i], points[i+1]
        dist = norm(p2 - p1)
        numpoints = dist*density
        for n in 1:numpoints
            push!(leg, Vec3((1 - (n-1)/numpoints)*p1 + (n-1)*p2/numpoints))
        end
        push!(legs, leg)
    end
    push!(legs[end], Vec3(points[end]))
    return vcat(legs...)
end

"""
    path(sf::StructureFactor, points::Vector; density = 10, index_labels=false, kwargs...)

Takes a list of ordered wave vectors, `points`, and extracts a path along a
lines connecting these wave vectors. The number of wave vectors sampled between
the specified points is determined by the `density` keyword, which determines
the number of points per inverse angstrom. If `index_labels` is set to `true`,
will return index label information which may be unpacked as `(; intensities,
qs, ωs)`. `intensities` will be a two dimension array, with the first index
corresponding to wave vectors (`q`) and the second to energy (`ω`).

All other keywords are shared with [`get_intensities`](@ref).
"""
function path(sf::StructureFactor, points::Vector; density = 10, index_labels=false, kwargs...)
    qs = path_points(Vec3.(points), density)
    intensities = Sunny.get_intensities(sf, qs; kwargs...) 
    if index_labels
        ωs = ωvals(sf)
        qs = map(q -> q.data, qs) # Return Tuples instead of StaticArrays
        return (; intensities, qs, ωs)
    end
    return intensities
end