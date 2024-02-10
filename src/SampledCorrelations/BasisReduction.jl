function phase_averaged_elements(data, q_absolute::Vec3, crystal::Crystal, ff_atoms, ::Val{NCorr}, ::Val{NAtoms}) where {NCorr, NAtoms}
    elems = zero(MVector{NCorr,ComplexF64})

    # Form factor
    ffs = ntuple(i -> compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute), Val{NAtoms}())

    # Overall phase factor for each site
    q = crystal.recipvecs \ q_absolute
    r = crystal.positions
    prefactor = ntuple(i -> ffs[i] * exp(- 2π*im * (q ⋅ r[i])) / sqrt(NAtoms), Val{NAtoms}())

    for j in 1:NAtoms, i in 1:NAtoms
        elems .+= (prefactor[i] * conj(prefactor[j])) .* view(data, :, i, j)
    end

    return SVector{NCorr,ComplexF64}(elems)
end

# Weighted average over variances to propagate error. This is essentially a
# "random phase" approximation, i.e., assumes the interference effects from
# phase averaging will average out. The weighting in the average comes from the
# form factors only. This approximation seems to give results very similar to
# those achieved by directly pulling out statistics during the sampling process.
# But it seems to me that the correct thing to do is simply to add all the
# variances from all the sites that contribute. Revisit this. This is currently
# not used, but would be employed in a parallel intensities-like pipeline for
# error propagation.
function error_basis_reduction(data, q_absolute::Vec3, _::Crystal, ff_atoms, ::Val{NCorr}, ::Val{NAtoms}) where {NCorr, NAtoms}
    elems = zero(MVector{NCorr,ComplexF64})
    ffs = ntuple(i -> compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute), NAtoms)

    for j in 1:NAtoms, i in 1:NAtoms
        elems .+= ffs[i] .* ffs[j] .* view(data, :, i, j)
    end

    return SVector{NCorr,Float64}(elems / (NAtoms*NAtoms))
end
