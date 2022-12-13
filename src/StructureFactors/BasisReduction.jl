# This can be optimized
function phase_averaged_elements(q_data, q::Vec3, cryst::Crystal, site_infos::Vector{SiteInfo})
    nelem = size(q_data, 1)
    elems = zero(SVector{nelem, ComplexF64})

    # Each q index ranges from 0 to 1. Then k (with appropriate periodicity)
    # lives in the first Brillouin zone.
    k = inv(cryst.lat_vecs)' * q

    for l in 1:nbasis(cryst)
        r = position(cryst, l)
        ffp1 = site_infos[l].ff_params
        ff1 = isnothing(ffp1) ? 1.0 : compute_form(norm(q), ffp1)
        for l′ in 1:nbasis(cryst)
            r′ = position(cryst, l′)
            ffp2 = site_infos[l′].ff_params
            ff2 = isnothing(ffp2) ? 1.0 : compute_form(norm(q), ffp2)
            phase =  exp(-im*(k ⋅ (r′ - r)))
            elems += phase * ff1 * ff2 * q_data[:, l, l′] 
        end
    end
    return elems
end
