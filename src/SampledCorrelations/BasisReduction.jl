function phase_averaged_elements(data, q_absolute::Vec3, crystal::Crystal, ff_atoms, ::Val{NCorr}, ::Val{NAtoms}) where {NCorr, NAtoms}
    elems = zero(MVector{NCorr,ComplexF64})
    ffs = ntuple(i -> compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute), NAtoms)

    # Real space position of each atom within the unit cell
    rs = ntuple(i -> crystal.latvecs * crystal.positions[i], NAtoms)

    for j in 1:NAtoms, i in 1:NAtoms
        phase = exp(im*(q_absolute ⋅ (rs[j] - rs[i])))
        elems .+= phase .* ffs[i] .* ffs[j] .* view(data, :, i, j)
    end

    return SVector{NCorr,ComplexF64}(elems)
end

function error_basis_reduction(data, q_absolute::Vec3, crystal::Crystal, ff_atoms, ::Val{NCorr}, ::Val{NAtoms}) where {NCorr, NAtoms}
    elems = zero(MVector{NCorr,ComplexF64})
    ffs = ntuple(i -> compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute), NAtoms)

    # Real space position of each atom within the unit cell
    rs = ntuple(i -> crystal.latvecs * crystal.positions[i], NAtoms)

    for j in 1:NAtoms, i in 1:NAtoms
        elems .+= ffs[i] .* ffs[j] .* view(data, :, i, j)
    end

    return SVector{NCorr,Float64}(elems)
end
