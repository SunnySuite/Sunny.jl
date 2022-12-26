################################################################################
# Basic functions for retrieving 𝒮(q, ω) values
################################################################################
function get_intensity(sf::StructureFactor, q; kwargs...) 
    if length(q) != 3
        error("Q point should have three components. If ")
    end
    return get_intensities(sf, [Vec3(q...)]; kwargs...)
end

function get_static_intensity(sf::StructureFactor, q; kwargs...)
    intensities = get_intensity(sf, q; kwargs...)
    return sum(intensities)
end

function Base.zeros(::Contraction{T}, args...) where T
    zeros(T, args...)
end

function get_intensities(sf::StructureFactor, q_targets::Array;
    interp = NoInterp(), contraction = Depolarize(), temp = nothing,
    negative_energies = false
) 
    nq = length(q_targets)
    ωs = negative_energies ? ωvals_all(sf) : ωvals(sf)
    nω = length(ωs) 
    contractor = contraction(sf)

    intensities = zeros(contractor, size(q_targets)..., nω)
    for iω in 1:nω
        for iq ∈ CartesianIndices(q_targets)
            q_target = convert(Vec3, q_targets[iq])
            qs, iqs = stencil_qs(sf.sfdata, q_target, interp)
            local_intensities = stencil_intensities(sf, qs, iqs, ωs[iω], iω, interp, contractor, temp)
            intensities[iq, iω] = interpolated_intensity(sf, q_target, qs, local_intensities, interp)
        end
    end

    return nq == 1 ? reshape(intensities, nω) : intensities
end

function get_static_intensities(sf::StructureFactor, q_targets::Array; kwargs...)
    dims = size(q_targets)
    if sum(dims) < 2
        error("To call get_static_intensities, must provide at least 2 Q values")
    end
    ndims = length(dims)
    intensities = get_intensities(sf, q_targets; kwargs...)
    println(size(intensities))
    static_intensities = sum(intensities, dims=(ndims+1,))

    return reshape(static_intensities, dims)
end


# Internal function for getting a single 𝒮(q, ω) intensity
function calc_intensity(sf::StructureFactor, q, iq, ω, iω, contractor, temp)
    (; crystal, site_infos, data) = sf.sfdata

    nelems, natoms = size(data, 1), size(data, 5)
    data_point = SArray{Tuple{nelems, natoms, natoms}, ComplexF64, 3, nelems*natoms*natoms}(
        data[:,iq,:,:,iω]
    )
    elems = phase_averaged_elements(data_point, q, crystal, site_infos)
    intensity = contract(elems, q, contractor)
    if !isnothing(temp)
        intensity *= classical_to_quantum(ω, temp)
    end

    return intensity
end


################################################################################
# Bulk extraction 
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
    index_labels=false
)
    qpoints = path_points(points, density)
    intensities = Sunny.get_intensities(sf, qpoints; interp, contraction, temp) 

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
