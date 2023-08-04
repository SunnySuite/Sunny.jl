abstract type InterpolationScheme{NumInterp} end
struct NoInterp <: InterpolationScheme{1} end
struct LinearInterp <: InterpolationScheme{8} end

#=
ddtodo: Explanation of interpolation "API"
=# 

function interpolated_intensity(::SampledCorrelations, _, _, stencil_intensities, ::NoInterp) 
    return only(stencil_intensities)
end

function interpolated_intensity(::SampledCorrelations, q_target, qs, stencil_intensities, ::LinearInterp) 
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


function stencil_points(sc::SampledCorrelations, q, ::NoInterp)

    # Each of the following lines causes a 32 byte allocation
    Ls = size(sc.samplebuf)[2:4] 
    m = round.(Int, Ls .* q)
    im = map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex{3}

    ## The following lines cause no allocations, but don't seem to be any faster.
    #     _, L1, L2, L3, _, _ = size(sc.samplebuf)
    #     m = (round(Int, L1*q[1]), round(Int, L2*q[2]), round(Int, L3*q[3]))
    #     im = CartesianIndex{3}(mod(m[1], L1)+1, mod(m[2], L2)+1, mod(m[3], L3)+1)

    return (m,), (im,)
end


function stencil_points(sc::SampledCorrelations, q, ::LinearInterp)
    Ls = size(sc.samplebuf)[2:4] 
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
function pruned_stencil_info(sc::SampledCorrelations, qs, interp::InterpolationScheme{N}) where N
    # Count the number of contiguous regions with unchanging values. If all
    # values are unique, returns the length of q_info. Note comparison is on m
    # values rather than index values and the m values are the first element of
    # the a tuple, that is, we're checking x[1] == y[1] in the map.
    m_info = map(q -> stencil_points(sc, q, interp), qs)
    numregions = sum(map((x,y) -> x[1] == y[1] ? 0 : 1, m_info[1:end-1], m_info[2:end])) + 1
    
    # Remove repeated stencil points and count number of instances of each
    ms_ref, idcs_ref = stencil_points(sc, qs[1], interp)
    ms_all  = fill(ntuple(x->zero(Vec3), N), numregions)
    ms_all[1] = ms_ref 
    idcs_all = fill(ntuple(x->CartesianIndex((-1,-1,-1)), N), numregions)
    idcs_all[1] = idcs_ref 
    counts = zeros(Int64, numregions)
    c = counts[1] = 1
    for q in qs[2:end] 
        ms, idcs = stencil_points(sc, q, interp)
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
    qs_all = map(ms_all) do ms
       map(m -> m ./ sc.latsize, ms) 
    end

    ks_all = map(qs_all) do qs
        map(q -> sc.crystal.recipvecs * q, qs) # FIXME: should q still by in sc.crystal, or in absolute units?
    end
    
    return (; qs_all, ks_all, idcs_all, counts)
end



"""
    intensities_interpolated(sc::SampledCorrelations, qs; interpolation = nothing, formula = intensity_formula(sc,:perp), negative_energies = false)

The basic function for retrieving ``𝒮(𝐪,ω)`` information from a
`SampledCorrelations`. Maps an array of wave vectors `qs` to an array of structure
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
function intensities_interpolated(sc::SampledCorrelations, qs;
    formula = intensity_formula(sc,:perp) :: ClassicalIntensityFormula,
    interpolation = :round,
    negative_energies = false,
    static_warn = true
)
    qs = Vec3.(qs)

    # If working on reshaped system, assume qs given as coordinates in terms of
    # reciprocal vectors of original crystal and convert them to qs in terms of
    # the reciprocal vectors of the reshaped crystal.
    if !isnothing(sc.origin_crystal)
        rvecs_reshaped = inv(sc.crystal.latvecs)'       # Note, leading 2π will cancel
        rvecs_origin = inv(sc.origin_crystal.latvecs)'
        qs = map(q -> rvecs_reshaped \ rvecs_origin * q, qs)
    end

    # Make sure it's a dynamical structure factor 
    if static_warn && size(sc.data, 7) == 1
        error("`intensities_interpolated` given a SampledCorrelations with no dynamical information. Call `instant_intensities_interpolated` to retrieve instantaneous (static) structure factor data.")
    end

    # Set up interpolation scheme
    interp = if interpolation == :round
        NoInterp()
    elseif interpolation == :linear
        LinearInterp()
    end

    # Precompute index information and preallocate
    ωvals = ωs(sc; negative_energies)
    nω = length(ωvals) 
    stencil_info = pruned_stencil_info(sc, qs, interp) 

    # Pre-allocate the output array based on the return type of the formula
    return_type = typeof(formula).parameters[1]
    intensities = zeros(return_type, size(qs)..., nω)
    
    # Call type stable version of the function
    intensities_interpolated!(intensities, sc, qs, ωvals, interp, formula, stencil_info, return_type)

    return intensities
end


# Actual intensity calculation
function intensities_interpolated!(intensities, sc::SampledCorrelations, q_targets::Array, ωvals, interp::InterpolationScheme{NInterp}, formula, stencil_info, T) where {NInterp}
    li_intensities = LinearIndices(intensities)
    ci_qs = CartesianIndices(q_targets)
    (; qs_all, ks_all, idcs_all, counts) = stencil_info 
    for (iω, ω) in enumerate(ωvals)
        iq = 0
        for (qs, ks, idcs, numrepeats) in zip(qs_all, ks_all, idcs_all, counts)
            local_intensities = SVector{NInterp, T}(formula.calc_intensity(sc, ks[n], idcs[n], iω) for n in 1:NInterp)
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(ci_qs[iq], iω)]
                intensities[idx] = interpolated_intensity(sc, q_targets[iq], qs, local_intensities, interp) 
            end
        end
    end
    return intensities
end

"""
    instant_intensities_interpolated(sc::SampledCorrelations, qs; kwargs...)

Return ``𝒮(𝐪)`` intensities at wave vectors `qs`. The functionality is very
similar to [`intensities_interpolated`](@ref), except the returned array has dimensions
identical to `qs`. If called on a `SampledCorrelations` with dynamical information,
i.e., ``𝒮(𝐪,ω)``, the ``ω`` information is integrated out.
"""
function instant_intensities_interpolated(sc::SampledCorrelations, qs; kwargs...)
    datadims = size(qs)
    ndims = length(datadims)
    vals = intensities_interpolated(sc, qs; static_warn=false, kwargs...)
    static_vals = sum(vals, dims=(ndims+1,))
    return reshape(static_vals, datadims)
end


"""
    connected_path_from_rlu(cryst, qs_rlu, density)

Returns a pair `(path, xticks)`. The first return value is a path in reciprocal
space that samples linearly between the wavevectors in `qs`. The elements in
`qs` are defined in reciprocal lattice units (RLU) associated with the
[`reciprocal_lattice_vectors``](@ref) for `cryst`. The `density` parameter has
units of inverse length, and controls the number of samples between elements of
`qs`.

The second return value `xticks` can be used for plotting. The `xticks` object
is itself a pair `(numbers, labels)`, which give the locations of the
interpolating ``q``-points and a human-readable string.
"""
function connected_path_from_rlu(cryst::Crystal, qs_rlu::Vector, density)
    @assert length(qs_rlu) >= 2 "The list `qs` should include at least two wavevectors."
    qs_rlu = Vec3.(qs_rlu)

    qs = Ref(cryst.recipvecs) .* qs_rlu

    path = Vec3[]
    markers = Int[]
    for i in 1:length(qs)-1
        push!(markers, length(path)+1)
        q1, q2 = qs[i], qs[i+1]
        dist = norm(q1 - q2)
        npoints = round(Int, dist*density)
        for n in 1:npoints
            push!(path, (1 - (n-1)/npoints)*q1 + (n-1)*q2/npoints)
        end
    end
    push!(markers, length(path)+1)
    push!(path, qs[end])

    labels = map(qs_rlu) do q_rlu
        "[" * join(number_to_math_string.(q_rlu), ",") * "]"
    end
    xticks = (markers, labels)

    return (path, xticks)
end
