function phase_averaged_elements(data, k::Vec3, sf::StructureFactor{N, NumCorr, NumBasis}, ffdata::Vector{FormFactor{FFType}}) where {FFType, N, NumCorr, NumBasis}
    elems = zero(SVector{NumCorr, ComplexF64})
    knorm = norm(k)
    ffs = ntuple(b -> compute_form(knorm, ffdata[b]), NumBasis)
    rs = ntuple(b -> sf.crystal.lat_vecs * sf.crystal.positions[b], NumBasis)

    for b2 in 1:NumBasis, b1 in 1:NumBasis
        phase = exp(im*(k ⋅ (rs[b2] - rs[b1])))
        elems += phase * ffs[b1] * ffs[b2] * view(data, :, b1, b2)  # This view allocates
    end

    return elems
end
