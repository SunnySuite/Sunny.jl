# This can be optimized -- note, also accumulates finite complex values where there should be none (on order of 1e-17)
function phase_averaged_elements(q_data, q::Vec3, cryst::Crystal, ffdata::Vector{Union{FormFactorParams, Nothing}})
    nelem = size(q_data, 1)
    elems = zero(SVector{nelem, ComplexF64})

    # If the q components range from 0 to 1, then k (with appropriate periodicity)
    # lives in the first Brillouin zone. 
    k = 2π*inv(cryst.lat_vecs)' * q

    for l in 1:nbasis(cryst)
        r = position(cryst, l)
        ffp1 = ffdata[l] 
        ff1 = isnothing(ffp1) ? 1.0 : compute_form(norm(k), ffp1)
        for l′ in 1:nbasis(cryst)
            r′ = position(cryst, l′)
            ffp2 = ffdata[l′]
            ff2 = isnothing(ffp2) ? 1.0 : compute_form(norm(k), ffp2)
            phase =  exp(-im*(k ⋅ (r′ - r)))
            elems += phase * ff1 * ff2 * @view(q_data[:, l, l′])
        end
    end
    return elems
end
