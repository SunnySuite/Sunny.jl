function phase_averaged_elements(data, k::Vec3, cryst::Crystal, ffdata::Vector{FormFactor{FFType}}, ::Val{NumElem}) where {NumElem, FFType}
    elems = zero(SVector{NumElem, ComplexF64})
    knorm = norm(k)
    for b2 in 1:nbasis(cryst)
        r2 = position(cryst, b2)
        ffp2 = ffdata[b2]
        ff2 = compute_form(knorm, ffp2)
        for b1 in 1:nbasis(cryst)
            r1 = position(cryst, b1)
            ffp1 = ffdata[b1] 
            ff1 = compute_form(knorm, ffp1)
            phase = exp(im*(k ⋅ (r2 - r1)))
            elems += phase * ff1 * ff2 * view(data, :, b1, b2)
        end
    end
    return elems
end
