function phase_averaged_elements(data, k::Vec3, sf::StructureFactor{N, NumCorr, NumBasis}, ffdata::Vector{FormFactor{FFType}}) where {FFType, N, NumCorr, NumBasis}
    cryst = sf.crystal
    elems = zero(SVector{NumCorr, ComplexF64})
    knorm = norm(k)

    ffs = ntuple(b -> compute_form(knorm, ffdata[b]), NumBasis)
    rs = ntuple(b -> position(cryst, b), NumBasis)
    for b2 in 1:NumBasis, b1 in 1:NumBasis
        phase = exp(im*(k ⋅ (rs[b2] - rs[b1])))
        elems += phase * ffs[b1] * ffs[b2] * view(data, :, b1, b2)  # This view allocates
    end

    # for b2 in 1:NumBasis
    #     ff2 = compute_form(knorm, ffdata[b2])
    #     r2 = position(cryst, b2)
    #     for b1 in 1:NumBasis
    #         ff1 = compute_form(knorm, ffdata[b1])
    #         r1 = position(cryst, b1)
    #         phase = exp(im*(k ⋅ (r2 - r1)))
    #         elems += phase * ff1 * ff2 * view(data, :, b1, b2)  # This view allocates
    #     end
    # end

    return elems
end
