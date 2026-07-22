function prefactors_for_phase_averaging(q_absolute::Vec3, recipvecs, positions, ffs, ::Val{NCorr}, ::Val{NPos}) where {NCorr, NPos}
    # Form factor
    ffvals = ntuple(i -> compute_form_factor(ffs[i], q_absolute⋅q_absolute), Val{NPos}())

    # Overall phase factor for each position
    q = recipvecs \ q_absolute
    r = positions
    prefactor = ntuple(i -> ffvals[i] * exp(- 2π*im * (q ⋅ r[i])), Val{NPos}())

    return prefactor
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
function error_basis_reduction(data, q_absolute::Vec3, _::Crystal, ffs, ::Val{NCorr}, ::Val{NPos}) where {NCorr, NPos}
    elems = zero(SVector{NCorr, ComplexF64})
    ffvals = ntuple(i -> compute_form_factor(ffs[i], q_absolute⋅q_absolute), NPos)

    for j in 1:NPos, i in 1:NPos
        elems += (ffvals[i] * ffvals[j]) * SVector{NCorr}(view(data, :, i, j))
    end

    return elems / (NPos*NPos)
end
