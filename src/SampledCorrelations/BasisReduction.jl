function phase_averaged_elements(data, k::Vec3, sc::SampledCorrelations{N}, ffdata::Vector{FormFactor{FFType}}, ::Val{NCorr}, ::Val{NAtoms}) where {FFType, N, NCorr, NAtoms} 
    elems = zero(MVector{NCorr,ComplexF64})
    knorm = norm(k)
    ffs = ntuple(i -> compute_form(knorm, ffdata[i]), NAtoms)
    rs = ntuple(i -> sc.crystal.latvecs * sc.crystal.positions[i], NAtoms)

    for j in 1:NAtoms, i in 1:NAtoms
        phase = exp(im*(k ⋅ (rs[j] - rs[i])))
        elems .+= phase .* ffs[i] .* ffs[j] .* view(data,:, i, j)
    end

    return SVector{NCorr,ComplexF64}(elems)
end
