################################################################################
# Basic functions for retrieving 𝒮(q, ω) values
################################################################################

# Internal function for getting a single 𝒮(q, ω) intensity
function calc_intensity(sf::StructureFactor, k, latidx, ω, iω, contractor, kT, ffdata)
    elems = phase_averaged_elements(view(sf.data,:,:,:,latidx,iω), k, sf, ffdata)
    intensity = contract(elems, k, contractor)
    return intensity * classical_to_quantum(ω, kT)
end

classical_to_quantum(ω, kT::Float64) = iszero(ω) ? 0.0 : ω/(kT*(1 - exp(-ω/kT)))
classical_to_quantum(ω, ::Nothing) = 1.0


# Note that requests for intensities often come in lists of nearby q values.
# Since the data is inherently discretized, this often results in repeated calls
# for values at the same discrete points. Since basis reduction is done for each
# of this calls, this results in a large amount of repeated calculation. This
# function analyzes repetitions in advance and prunes out repetitions. This is
# ugly, but the speedup when tested on a few simple, realistic examples was
# 3-5x.
function pruned_stencil_info(sf::StructureFactor, qs, interp::InterpolationScheme{N}) where N
    # Count the number of contiguous regions with unchanging values. If all
    # values are unique, returns the length of q_info. Note comparison is on m
    # values rather than index values and the m values are the first element of
    # the a tuple, that is, we're checking x[1] == y[1] in the map.
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
    recip_vecs = 2π*inv(sf.crystal.lat_vecs)'
    latsize = size(sf.samplebuf)[2:4]
    qs_all = map(ms_all) do ms
       map(m -> m ./ latsize, ms) 
    end

    ks_all = map(qs_all) do qs
        map(q -> recip_vecs * q, qs)
    end
    
    return (; qs_all, ks_all, idcs_all, counts)
end


"""
    intensities(sf::StructureFactor, qs, mode; interpolation = nothing,
                    kT = nothing, formfactors = nothing, negative_energies = false)

The basic function for retrieving ``𝒮(𝐪,ω)`` information from a
`StructureFactor`. Maps an array of wave vectors `qs` to an array of structure
factor intensities, including an additional energy index. The values of ``ω``
associated with the energy index can be retrieved by calling [`ωs`](@ref).
The three coordinates of each wave vector are measured in reciprocal lattice
units, i.e., multiples of the reciprocal lattice vectors.

- `mode`: Should be one of `:trace`, `:perp`, or `:full`. Determines an optional
    contraction on the indices ``α`` and ``β`` of ``𝒮^{αβ}(q,ω)``. Setting
    `trace` yields ``∑_α 𝒮^{αα}(q,ω)``. Setting `perp` will employ a
    polarization correction on the traced value. Setting `full` will return all
    elements ``𝒮^{αβ}(𝐪,ω)`` with contraction.
- `interpolation`: Since ``𝒮(𝐪, ω)`` is calculated on a finite lattice, data is
    only available at discrete wave vectors. By default, Sunny will round a
    requested `q` to the nearest available wave vector. Linear interpolation can
    be applied by setting `interpolation=:linear`.
- `kT`: If a temperature is provided, the intensities will be rescaled by a
    temperature- and ω-dependent classical-to-quantum factor. `kT` should be
    specified when making comparisons with spin wave calculations or
    experimental data.
- `formfactors`: To apply form factor corrections, provide this keyword with a
    vector of `FormFactor`s, one for each unique site in the unit cell. Sunny
    will symmetry propagate the results to all equivalent sites.
- `negative_energies`: If set to `true`, Sunny will return the periodic
    extension of the energy axis. Most users will not want this.
"""
function intensities(sf::StructureFactor, qs, mode;
    interpolation = :none,
    kT = nothing,
    formfactors = nothing,
    negative_energies = false,
    static_warn = true
)
    qs = Vec3.(qs)

    # Make sure it's a dynamical structure factor 
    if static_warn && size(sf.data, 7) == 1
        error("`intensities` given a static structure factor. Call `static_intensities` to retrieve static structure factor data.")
    end

    # Set up interpolation scheme
    interp = if interpolation == :none
        NoInterp()
    elseif interpolation == :linear
        LinearInterp()
    end

    # Set up element contraction
    contractor = if mode == :trace
        Trace(sf)
    elseif mode == :perp
        DipoleFactor(sf)
    elseif mode == :full
        FullTensor(sf)
    elseif typeof(mode) <: Tuple{Int, Int}
        Element(sf, mode)
    end

    # Propagate form factor information (if any)
    if isnothing(formfactors)
        formfactors = [FormFactor{EMPTY_FF}(; atom) for atom in unique(sf.crystal.classes)]
    end
    ffdata = propagate_form_factors(sf.crystal, formfactors)

    # Precompute index information and preallocate
    ωvals = ωs(sf; negative_energies)
    nω = length(ωvals) 
    intensities = zeros(contractor, size(qs)..., nω)
    stencil_info = pruned_stencil_info(sf, qs, interp) 
    
    # Call type stable version of the function
    intensities!(intensities, sf, qs, ωvals, interp, contractor, kT, ffdata, stencil_info) #ddtodo: Track down allocations

    # ddtodo: See if worth it to apply classical-to-quantum rescaling here instead of inside loop (removes branching)

    return intensities
end


# Type stable version
function intensities!(intensities, sf::StructureFactor, q_targets::Array, ωvals, interp::InterpolationScheme, contraction::Contraction{T}, temp, ffdata, stencil_info) where {T}
    li_intensities = LinearIndices(intensities)
    ci_qs = CartesianIndices(q_targets)
    (; qs_all, ks_all, idcs_all, counts) = stencil_info 
    for (iω, ω) in enumerate(ωvals)
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



"""
    static_intensities(sf::StructureFactor, qs, mode; kwargs...)

Return the static structure factor intensities at wave vectors `qs`. The
functionality is very similar to [`intensities`](@ref), except the returned
array has dimensions identical to `qs`. The energy axis has been integrated out.
"""
function static_intensities(sf::StructureFactor, qs, mode; kwargs...)
    datadims = size(qs)
    ndims = length(datadims)
    vals = intensities(sf, qs, mode; static_warn=false, kwargs...)
    static_vals = sum(vals, dims=(ndims+1,))
    return reshape(static_vals, datadims)
end


"""
    connected_path(qs::Vector, density)

Takes a list of wave vectors, `qs`, and builds an expanded list of wave vectors
that traces a path through the provided points. Also returned is a list of
marker indices corresponding to the intput points. The `density` parameter
controls the frequency of sampling.
"""
function connected_path(qs::Vector, density)
    @assert length(qs) >= 2 "The list `qs` should include at least two wavevectors."

    qs = Vec3.(qs)
    path = Vec3[]
    markers = Int[]
    for i in 1:length(qs)-1
        push!(markers, length(path)+1)
        q1, q2 = qs[i], qs[i+1]
        dist = norm(q2 - q1)
        npoints = round(Int, dist*density)
        for n in 1:npoints
            push!(path, (1 - (n-1)/npoints)*q1 + (n-1)*q2/npoints)
        end
    end
    push!(markers, length(path)+1)
    push!(path, qs[end])
    return (path, markers)
end


"""
    lorentzian(x, η) 

Returns ``η^2/(x^2 + η^2)``.
"""
lorentzian(x, η) = η^2/(x^2 + η^2)

"""
    broaden_energy(sf::StructureFactor, vals, kernel::Function; negative_energies=false)

Performs a real-space convolution along the energy axis of an array of
intensities. Assumes the format of the intensities array corresponds to what
would be returned by [`intensities`](@ref). `kernel` must be a function that
takes two numbers: `kernel(ω, ω₀)`, where `ω` is a frequency, and `ω₀` is the
center frequency of the kernel. Sunny provides [`lorentzian`](@ref)
for the most common use case:

```
newvals = broaden_energy(sf, vals, (ω, ω₀) -> lorentzian(ω-ω₀, 0.2))
```
"""
function broaden_energy(sf::StructureFactor, is, kernel::Function; negative_energies=false)
    dims = size(is)
    ωvals = ωs(sf; negative_energies)
    out = zero(is)
    for (ω₀i, ω₀) in enumerate(ωvals)
        for (ωi, ω) in enumerate(ωvals)
            for qi in CartesianIndices(dims[1:end-1])
                out[qi,ωi] += is[qi,ω₀i]*kernel(ω, ω₀)
            end
        end
    end
    return out
end