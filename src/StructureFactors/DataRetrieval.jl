################################################################################
# Basic functions for retrieving 𝒮(q, ω) values
################################################################################

# Function for getting a single 𝒮(q, ω) intensity -- primarily internal
function calc_intensity(sf::StructureFactor, q, iq, ω, iω, contractor, temp, ffdata)
    (; crystal, data) = sf.sfdata
    elems = phase_averaged_elements(@view(data[:,:,:,iq, iω]), q, crystal, ffdata)
    intensity = contract(elems, q, contractor)
    if !isnothing(temp)
        intensity *= classical_to_quantum(ω, temp)
    end
    return intensity
end

function Base.zeros(::Contraction{T}, args...) where T
    zeros(T, args...)
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
                2*I(3)  # This is irrelevant -- see note above call to this function
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
If the user provides a custom basis (newbasis), it is assumed they have provided q-values
in terms of this basis. This function converts the user-provided q-values into coordinates
with respect to the reciprocal lattice vectors.
"""
function change_basis(qs, newbasis)
    return map(qs) do q
        sum(newbasis .* q) # Make newbasis a matrix instead of list of Vec3s?
    end
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
function prune_stencil_qs(sfd, q_targets, interp::InterpolationScheme{N}) where N
    q_info = map(q -> stencil_qs(sfd, q, interp), q_targets)
    # Count the number of contiguous regions with unchanging values.
    # Note: if all values are unique, returns the length of q_info.
    numregions = sum(map((x,y) -> x == y ? 0 : 1, q_info[1:end-1], q_info[2:end])) + 1

    counts = zeros(Int64, numregions)
    qis_all = Array{NTuple{N, CartesianIndex{3}}}(undef, numregions)
    qs_all = Array{NTuple{N, Vec3}}(undef, numregions) 

    qs, qis = stencil_qs(sfd, q_targets[1], interp)

    qis_all[1] = qis_ref = CartesianIndex.(qis)
    qs_all[1] = qs
    c = counts[1] = 1
    for q in q_targets[2:end]
        qs, qis = stencil_qs(sfd, q, interp)
        if qis != qis_ref
            qis_ref = qis
            c += 1
            qis_all[c] = qis
            qs_all[c] = qs
        end
        counts[c] += 1
    end
    # @assert sum(counts) == length(q_info)
    return (; counts, qis_all, qs_all)
end

"""
"""
function get_intensities(sf::StructureFactor, q_targets::Array;
    interp = NoInterp(), contraction = Trace(), temp = nothing,
    formfactors = nothing, negative_energies = false, newbasis = nothing,
) 
    nq = length(q_targets)
    ωs = negative_energies ? ωvals_all(sf) : ωvals(sf)
    nω = length(ωs) 
    contractor = contraction(sf)
    ffdata = propagate_form_factors(sf, formfactors)
    if !isnothing(newbasis)
        q_targets = change_basis(q_targets, newbasis)
    end
    q_targets = convert.(Vec3, q_targets)

    intensities = zeros(contractor, size(q_targets)..., nω)
    li_intensities = LinearIndices(intensities)
    ci_qs = CartesianIndices(q_targets)

    (; counts, qis_all, qs_all) = prune_stencil_qs(sf.sfdata, q_targets, interp)
    @time for iω in 1:nω
        iq = 0
        for (c, numrepeats) in enumerate(counts)
            qs, qis = qs_all[c], qis_all[c]
            local_intensities = stencil_intensities(sf, qs, qis, ωs[iω], iω, interp, contractor, temp, ffdata)
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(ci_qs[iq], iω)]
                intensities[idx] = interpolated_intensity(sf, q_targets[iq], qs, local_intensities, interp)
            end
        end
        # @assert iq == length(q_targets)
    end

    return nq == 1 ? reshape(intensities, nω) : intensities
end


function get_intensity(sf::StructureFactor, q; kwargs...) 
    if length(q) != 3
        error("Q point should have three components.")
    end
    return get_intensities(sf, [Vec3(q...)]; kwargs...)
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
    density = 10, interp=Sunny.NoInterp(), contraction=Sunny.depolarize, temp=nothing,
    index_labels=false, kwargs...
)
    qpoints = path_points(points, density)
    intensities = Sunny.get_intensities(sf, qpoints; interp, contraction, temp, kwargs...) 

    if index_labels
        ωs = ωvals(sf)
        return (; intensities, qpoints, ωs)
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
    density = 10, interp=Sunny.NoInterp(), contraction=Sunny.depolarize, temp=nothing,
    index_labels=false
)
    qpoints = static_slice_points(p1, p2, z, density)
    intensities = Sunny.get_static_intensities(sf, qpoints; interp, contraction, temp)

    if index_labels
        return (; intensities, qpoints)
    end
    return intensities
end
