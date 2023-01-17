function phase_averaged_elements(data, k::Vec3, cryst::Crystal, ffdata::Vector{Union{FormFactor, Nothing}}, ::Val{NumElem}) where NumElem
    elems = zero(SVector{NumElem, ComplexF64})
    for b2 in 1:nbasis(cryst)
        r2 = position(cryst, b2)
        ffp2 = ffdata[b2]
        ff2 = isnothing(ffp2) ? 1.0 : compute_form(norm(k), ffp2)
        for b1 in 1:nbasis(cryst)
            r1 = position(cryst, b1)
            ffp1 = ffdata[b1] 
            ff1 = isnothing(ffp1) ? 1.0 : compute_form(norm(k), ffp1)
            phase = exp(im*(k ⋅ (r2 - r1)))
            elems += phase * ff1 * ff2 * view(data, :, b1, b2)
        end
    end
    return elems
end
