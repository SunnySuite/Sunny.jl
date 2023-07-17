abstract type InterpolationScheme{NumInterp} end
struct NoInterp <: InterpolationScheme{1} end
struct LinearInterp <: InterpolationScheme{8} end

#=
ddtodo: Explanation of interpolation "API"
=# 

function interpolated_intensity(::StructureFactor, _, _, stencil_intensities, ::NoInterp) 
    return only(stencil_intensities)
end

function interpolated_intensity(::StructureFactor, q_target, qs, stencil_intensities, ::LinearInterp) 
    q000,    _,    _,    _,    _,    _,    _, q111 = qs 
    c000, c100, c010, c110, c001, c101, c011, c111 = stencil_intensities

    x, y, z = q_target
    x0, y0, z0 = q000 
    x1, y1, z1 = q111 

    xd = (x-x0)/(x1-x0)
    yd = (y-y0)/(y1-y0)
    zd = (z-z0)/(z1-z0)

    c00 = c000*(1-xd) + c100*xd
    c01 = c001*(1-xd) + c101*xd
    c10 = c010*(1-xd) + c110*xd
    c11 = c011*(1-xd) + c111*xd

    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd

    return c0*(1-zd) + c1*zd
end


function stencil_points(sf::StructureFactor, q, ::NoInterp)

    # Each of the following lines causes a 32 byte allocation
    Ls = size(sf.samplebuf)[2:4] 
    m = round.(Int, Ls .* q)
    im = map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex{3}

    ## The following lines cause no allocations, but don't seem to be any faster.
    #     _, L1, L2, L3, _, _ = size(sf.samplebuf)
    #     m = (round(Int, L1*q[1]), round(Int, L2*q[2]), round(Int, L3*q[3]))
    #     im = CartesianIndex{3}(mod(m[1], L1)+1, mod(m[2], L2)+1, mod(m[3], L3)+1)

    return (m,), (im,)
end


function stencil_points(sf::StructureFactor, q, ::LinearInterp)
    Ls = size(sf.samplebuf)[2:4] 
    base = map(x -> floor(Int64, x), Ls .* q) 
    offsets = (
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 0), 
        (1, 1, 0), 
        (0, 0, 1),
        (1, 0, 1), 
        (0, 1, 1),
        (1, 1, 1)
    )
    ms = map(x -> x .+ base, offsets) 
    ims = map(ms) do m
        map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex{3}
    end
    return ms, ims
end

# Note that requests for intensities often come in lists of nearby q values.
# Since the data is inherently discretized, this often results in repeated calls
# for values at the same discrete points. Since basis reduction is done for each
# of this calls, this results in a large amount of repeated calculation. This
# function analyzes repetitions in advance and prunes them out. This is
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
    recip_vecs = 2π*inv(sf.crystal.latvecs)'  # Note, qs will be in terms of sf.crystal by this point, not origin_crystal
    qs_all = map(ms_all) do ms
       map(m -> m ./ sf.latsize, ms) 
    end

    ks_all = map(qs_all) do qs
        map(q -> recip_vecs * q, qs)
    end
    
    return (; qs_all, ks_all, idcs_all, counts)
end



"""
    intensities_interpolated(sf::StructureFactor, qs; interpolation = nothing, formula = intensity_formula(sf,:perp), negative_energies = false)

The basic function for retrieving ``𝒮(𝐪,ω)`` information from a
`StructureFactor`. Maps an array of wave vectors `qs` to an array of structure
factor intensities, including an additional energy index. The values of ``ω``
associated with the energy index can be retrieved by calling [`ωs`](@ref). The
three coordinates of each wave vector are measured in reciprocal lattice units,
i.e., multiples of the reciprocal lattice vectors.

- `interpolation`: Since ``𝒮(𝐪, ω)`` is calculated on a finite lattice, data
    is only available at discrete wave vectors. By default, Sunny will round a
    requested `q` to the nearest available wave vector. Linear interpolation can
    be applied by setting `interpolation=:linear`.
- `negative_energies`: If set to `true`, Sunny will return the periodic
    extension of the energy axis. Most users will not want this.
"""
function intensities_interpolated(sf::StructureFactor, qs;
    formula = intensity_formula(sf,:perp) :: ClassicalIntensityFormula,
    interpolation = :round,
    negative_energies = false,
    static_warn = true
)
    qs = Vec3.(qs)

    # If working on reshaped system, assume qs given as coordinates in terms of
    # reciprocal vectors of original crystal and convert them to qs in terms of
    # the reciprocal vectors of the reshaped crystal.
    if !isnothing(sf.origin_crystal)
        rvecs_reshaped = inv(sf.crystal.latvecs)'       # Note, leading 2π will cancel
        rvecs_origin = inv(sf.origin_crystal.latvecs)'
        qs = map(q -> rvecs_reshaped \ rvecs_origin * q, qs)
    end

    # Make sure it's a dynamical structure factor 
    if static_warn && size(sf.data, 7) == 1
        error("`intensities_interpolated` given a StructureFactor with no dynamical information. Call `instant_intensities_interpolated` to retrieve instantaneous (static) structure factor data.")
    end

    # Set up interpolation scheme
    interp = if interpolation == :round
        NoInterp()
    elseif interpolation == :linear
        LinearInterp()
    end

    # Precompute index information and preallocate
    ωvals = ωs(sf; negative_energies)
    nω = length(ωvals) 
    stencil_info = pruned_stencil_info(sf, qs, interp) 
    intensities = zeros(formula.calc_intensity.return_type, size(qs)..., nω)
    
    # Call type stable version of the function
    intensities_interpolated!(intensities, sf, qs, ωvals, interp, formula, stencil_info, formula.calc_intensity.return_type)

    return intensities
end


# Actual intensity calculation
function intensities_interpolated!(intensities, sf::StructureFactor, q_targets::Array, ωvals, interp::InterpolationScheme{NInterp}, formula, stencil_info, T) where {NInterp}
    li_intensities = LinearIndices(intensities)
    ci_qs = CartesianIndices(q_targets)
    (; qs_all, ks_all, idcs_all, counts) = stencil_info 
    for (iω, ω) in enumerate(ωvals)
        iq = 0
        for (qs, ks, idcs, numrepeats) in zip(qs_all, ks_all, idcs_all, counts)
            local_intensities = SVector{NInterp, T}(formula.calc_intensity(sf, ks[n], idcs[n], iω) for n in 1:NInterp)
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
    instant_intensities_interpolated(sf::StructureFactor, qs; kwargs...)

Return ``𝒮(𝐪)`` intensities at wave vectors `qs`. The functionality is very
similar to [`intensities_interpolated`](@ref), except the returned array has dimensions
identical to `qs`. If called on a `StructureFactor` with dynamical information,
i.e., ``𝒮(𝐪,ω)``, the ``ω`` information is integrated out.
"""
function instant_intensities_interpolated(sf::StructureFactor, qs; kwargs...)
    datadims = size(qs)
    ndims = length(datadims)
    vals = intensities_interpolated(sf, qs; static_warn=false, kwargs...)
    static_vals = sum(vals, dims=(ndims+1,))
    return reshape(static_vals, datadims)
end


"""
    connected_path(recip_vecs, qs::Vector, density)

Takes a list of wave vectors, `qs`, and builds an expanded list of wave vectors
that traces a path through the provided points. Also returned is a list of
marker indices corresponding to the input points. The `density` parameter is
given in samples per inverse Å.

Instead of `recip_vecs`, the first argument may be either a `StructureFactor` or
a `SpinWaveTheory`.
"""
function connected_path(recip_vecs, qs::Vector, density)
    @assert length(qs) >= 2 "The list `qs` should include at least two wavevectors."
    qs = Vec3.(qs)

    path = Vec3[]
    markers = Int[]
    for i in 1:length(qs)-1
        push!(markers, length(path)+1)
        q1, q2 = qs[i], qs[i+1]
        dist = norm(recip_vecs*(q1 - q2))
        npoints = round(Int, dist*density)
        for n in 1:npoints
            push!(path, (1 - (n-1)/npoints)*q1 + (n-1)*q2/npoints)
        end
    end
    push!(markers, length(path)+1)
    push!(path, qs[end])

    return (path, markers)
end
connected_path(sf::StructureFactor, qs::Vector, density) = connected_path(2π*inv(sf.crystal.latvecs)', qs, density)
connected_path(sw::SpinWaveTheory, qs::Vector, density) = connected_path(sw.recipvecs_chem, qs, density)

