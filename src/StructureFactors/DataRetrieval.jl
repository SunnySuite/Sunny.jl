function get_intensities(sf::StructureFactor, qs::Vector;
    contraction = depolarize,
    c2q_temp = nothing,
)
    (; crystal, idx_info, site_infos, data) = sf.sfdata

    # find q values and indices
    qs = sort(nearest_q.(Ref(sf.sfdata), qs), by=x->x[2])
    ωs = ωvals(sf)

    # Reduce data at given q values
    nelems, _, _, _, _, natoms, nω = size(data)
    nω = div(nω, 2) + 1  # Only take positive energies
    intensities = zeros(length(qs), nω)
    # data_point = zeros(ComplexF64, nelems, natoms, natoms)
    for iω in 1:nω
        for (n, qdata) in enumerate(qs)
            q, iq = qdata
            data_point = SArray{Tuple{nelems, natoms, natoms}, ComplexF64, 3, nelems*natoms*natoms}(
                data[:,iq,:,:,iω]
            )
            # @. data_point = @view sf.sfdata.data[:,iq,:,:,iω]
            elems = phase_averaged_elements(data_point, q, crystal, site_infos)
            intensities[n, iω] = contraction(elems, q, idx_info)
            if !isnothing(c2q_temp)
                intensities[n, iω] *= classical_to_quantum(ωs[iω], c2q_temp)
            end
        end
    end

    return intensities
end

function get_intensity(sf::StructureFactor, q; kwargs...)
    q = convert(Vec3, q)
    get_intensities(sf, [q]; kwargs...)
end


function get_intensity_grid(sf::StructureFactor;
    contraction = depolarize,
    c2q_temp = nothing,
)
    (; crystal, idx_info, site_infos, data) = sf.sfdata

    # Get data point labels and indices
    qs = qvals(sf)
    ωs = ωvals(sf)

    # Reduce all data 
    nelems, na, nb, nc, _, natoms, nω = size(data)
    nω = div(nω, 2) + 1  # Only take positive energies
    intensities = zeros(na, nb, nc, nω)
    data_point = zeros(ComplexF64, nelems, natoms, natoms)
    for iω in 1:nω
        for iq in CartesianIndices((na, nb, nc))
            data_point = SArray{Tuple{nelems, natoms, natoms}, ComplexF64, 3, nelems*natoms*natoms}(
                data[:,iq,:,:,iω]
            )
            # @. data_point = @view sf.sfdata.data[:,iq,:,:,iω]
            elems = phase_averaged_elements(data_point, qs[iq], crystal, site_infos)
            intensities[iq, iω] = contraction(elems, qs[iq], idx_info)
            if !isnothing(c2q_temp)
                intensities[iq, iω] *= classical_to_quantum(ωs[iω], c2q_temp)
            end
        end
    end

    return intensities
end

## get_intensities with checks

# Approach: 
#    - Construct a list of qis from qs (Remove repeated or just do repeated calculations?)
#    - Extract all of these at once in a big slice
#    - Perform reductions and contractions efficiently
#    - Then interpolate resulting intensities

