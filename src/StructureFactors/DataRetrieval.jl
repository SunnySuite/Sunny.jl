################################################################################
# Basic functions for retrieving 𝒮(q, ω) values
################################################################################

# Function for getting a single 𝒮(q, ω) intensity -- primarily internal -- really wants k, not m
function calc_intensity(sf::StructureFactor{N, NumCorr}, m, im, ω, iω, contractor, temp, ffdata) where {N, NumCorr}
    (; crystal, data) = sf.sfdata
    k = 2π*inv(crystal.lat_vecs)' * (m ./ sf.sftraj.sys.latsize)
    elems = phase_averaged_elements(view(data,:,:,:,im,iω), k, crystal, ffdata, Val(NumCorr))
    intensity = contract(elems, k, contractor)
    if !isnothing(temp)
        intensity *= classical_to_quantum(ω, temp)
    end
    return intensity
end


"""
Propagates from the form factor information, provided when requesting intensities,
to all symmetry equivalent sites. Works on the back of `all_symmetry_related_couplings`,
which also takes g-factor information, which is irrelevent here. We could perhaps write a simpler
symmetry propagation function to make this code nicer, but it would only ever be used here.
"""
function propagate_form_factors(sf::StructureFactor, form_factors)
    sys = sf.sftraj.sys
    natoms = size(sys.dipoles, 4)
    all_form_factors = Vector{Union{FormFactor, Nothing}}(nothing, natoms)
    if !isnothing(form_factors)
        specified_atoms = Int[]
        for form_factor in form_factors 
            atom = form_factor.atom
            # Using all_symmetry_related_couplings for convenience -- don't need transformed gs
            (sym_bonds, sym_gs) = all_symmetry_related_couplings(
                sys.crystal,
                Bond(atom, atom, [0,0,0]),
                2*I(3)  # This is irrelevant -- see note above this function
            )
            for (sym_bond, _) in zip(sym_bonds, sym_gs)
                sym_atom = sym_bond.i
                if sym_atom in specified_atoms
                    error("Provided `FormFactor` information for two symmetry equivalent sites!")
                else
                    push!(specified_atoms, sym_atom)
                end
                all_form_factors[sym_atom] = form_factor
            end
        end
    end
    return all_form_factors
end


"""
Note that requests for intensities often come in lists of nearby q values. Since the data
is inherently discretized, this often results in repeated calls for values at the same 
discrete points. Since basis reduction is done for each of this calls, this results in 
a large amount of repeated calculation. This analyzes in advance whenever the raw data 
that must undergo basis reduction changes as one iterates through a list of q values. This
information enables us to avoid some repeated identical calculations.

This is ugly, but the speedup when tested on a few simple, realistic examples was 3-5x.
"""
function prune_stencil_points(sfd, qs, interp::InterpolationScheme{N}) where N
    # Count the number of contiguous regions with unchanging values.
    # If all values are unique, returns the length of q_info.
    m_info = map(q -> stencil_points(sfd, q, interp), qs)
    numregions = sum(map((x,y) -> x[1] == y[1] ? 0 : 1, m_info[1:end-1], m_info[2:end])) + 1

    ms_ref, ims_ref = stencil_points(sfd, qs[1], interp)
    ms_all  = fill(ntuple(x->zero(Vec3), N), numregions)
    ms_all[1] = ms_ref 
    ims_all = fill(ntuple(x->CartesianIndex((-1,-1,-1)), N), numregions)
    ims_all[1] = ims_ref 
    counts = zeros(Int64, numregions)
    c = counts[1] = 1

    for q in qs[2:end] 
        ms, ims = stencil_points(sfd, q, interp)
        if ms != ms_ref
            ms_ref = ms 
            c += 1
            ms_all[c] =  ms
            ims_all[c] = ims 
        end
        counts[c] += 1
    end
    
    @assert sum(counts) == length(m_info)
    ms, ims = ms_all, ims_all

    return (; ms, ims, counts)
end

function precompute_form_factors(sf::StructureFactor, ffdata, qs)
    cryst = sf.sftraj.sys.crystal
    nb = nbasis(cryst)
    @assert length(ffdata) == nb
    ffdata_static = fill(ntuple(x->0.0, nb), size(qs))
    for (i, q) in enumerate(qs)
        k = 2π*inv(cryst.lat_vecs)' * q
        ffdata_static[i] = ntuple(x->isnothing(ffdata[i]) ? 1.0 : compute_form(norm(k), ffdata[i]), nb)
    end
    return ffdata_static
end

# Type stable version
function get_intensities(sf::StructureFactor, q_targets::Array, interp::InterpolationScheme, contraction::Contraction{T};
    temp = nothing, formfactors = nothing, negative_energies = false, newbasis = nothing,
) where {T}
    ωs = negative_energies ? ωvals_all(sf) : ωvals(sf)
    nω = length(ωs) 
    ffdata = propagate_form_factors(sf, formfactors)
    q_targets = Vec3.(q_targets) 
    if !isnothing(newbasis)
        newbasis = Mat3(newbasis)
        q_targets = map(q -> newbasis*q, q_targets)
    end

    intensities = zeros(T, size(q_targets)..., nω)
    li_intensities = LinearIndices(intensities)
    ci_qs = CartesianIndices(q_targets)
    (; counts, ims, ms) = prune_stencil_points(sf.sfdata, q_targets, interp) # needs (and gets) user qs

    for iω in 1:nω
        iq = 0
        for (c, numrepeats) in enumerate(counts)
            m, im = ms[c], ims[c]
            local_intensities = stencil_intensities(sf, m, im, ωs[iω], iω, interp, contraction, temp, ffdata) # really needs k, not q or m, for both basis reduction and contraction
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(ci_qs[iq], iω)]
                intensities[idx] = interpolated_intensity(sf, q_targets[iq], m, local_intensities, interp) # really needs q (divide by L inside)
            end
        end
        # @assert iq == length(q_targets)
    end

    return intensities
end

# This is the version of the function that will typically
# be called by the user. It sets up the appropriate data types
# and calls the type stable internal version of the function.
function get_intensities(sf::StructureFactor, qs::Array;
    contraction = :trace, interpolation = nothing, kwargs...)

    interp = if isnothing(interpolation) 
        NoInterp()
    elseif interpolation == :linear
        LinearInterp()
    end

    contract = if contraction == :trace
        Trace(sf)
    elseif contraction == :depolarize
        Depolarize(sf)
    elseif typeof(contraction) <: Tuple{Int, Int}
        Element(sf, contraction)
    end

    # Call type stable version of the function
    return get_intensities(sf, qs, interp, contract; kwargs...)
end



function get_intensity(sf::StructureFactor, q; kwargs...) 
    if length(q) != 3
        error("Q point should have three components.")
    end
    return get_intensities(sf, [Vec3(q)]; kwargs...)'
end

function get_static_intensity(sf::StructureFactor, q; kwargs...)
    intensities = get_intensity(sf, q; kwargs...)
    return sum(intensities)
end

function get_static_intensities(sf::StructureFactor, q_targets::Array; kwargs...)
    dims = size(q_targets)
    if sum(dims) < 2  
        error("To call get_static_intensities, must provide at least 2 Q values. For a single point, call `get_static_intensity`.")
    end
    ndims = length(dims)
    intensities = get_intensities(sf, q_targets; kwargs...)
    static_intensities = sum(intensities, dims=(ndims+1,))

    return reshape(static_intensities, dims)
end


################################################################################
# Functions for pulling out large amounts of data at once (paths, 2D slices) 
# NOTE: 2D slice functionality for testing only. This needs to be developed in
# a general way, probably in conjunction with SF plotting tools.
################################################################################
function intensity_grid(sf::StructureFactor;
    bzsize=(1,1,1), negative_energies = false, index_labels = false, kwargs...
)
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


function path(sf::StructureFactor, points::Vector; 
    density = 10, index_labels=false, kwargs...
)
    qs = path_points(Vec3.(points), density)
    intensities = Sunny.get_intensities(sf, qs; kwargs...) 
    if index_labels
        ωs = ωvals(sf)
        qs = map(q -> q.data, qs)
        return (; intensities, qs, ωs)
    end
    return intensities
end


function static_slice_points(p1, p2, z, density)
    dx, dy = p2[1] - p1[1], p2[2] - p1[2] 
    nx, ny = round(Int, dx*density), round(Int, dy*density) 
    points = zeros(Vec3, nx, ny)
    for i in CartesianIndices(points)
        x, y = i.I
        points[i] = Vec3(
            (1-((x-1)/nx))*p1[1] + ((x-1)/nx)*p2[1],
            (1-((y-1)/ny))*p1[2] + ((y-1)/ny)*p2[2],
            z
        )
    end
    return points
end

function static_slice(sf::StructureFactor, p1, p2, z = 0.0; 
    density = 10, index_labels=false, kwargs...
)
    qpoints = static_slice_points(p1, p2, z, density)
    intensities = Sunny.get_static_intensities(sf, qpoints; kwargs...)

    if index_labels
        return (; intensities, qpoints)
    end
    return intensities
end
