abstract type AbstractInteractionCUDA end   # Subtype this for actual internal GPU implementations


struct ExternalFieldCUDA <: AbstractInteractionCUDA
    effBs :: CUDA.CuVector{Vec3}  # |S_b|gᵀB for each basis index b
end

function ExternalFieldCUDA(ext_field::ExternalField, sites_info::Vector{SiteInfo})
    # As E = -∑_i 𝐁^T g 𝐒_i, we can precompute effB = g^T S B, so that
    #  we can compute E = -∑_i effB ⋅ 𝐬_i during simulation.
    # However, S_i may be basis-dependent, so we need to store an effB
    #  per sublattice.
    effBs = [site.g' * site.S * ext_field.B for site in sites_info]
    ExternalFieldCUDA(effBs)
end


function energy(spins::CUDA.CuArray{Vec3}, field::ExternalFieldCUDA)
    effBs = field.effBs
    nb = length(effBs)
    effBs = reshape(effBs, nb, ones(Int, ndims(spins)-1)...)
    return -sum(effBs .⋅ spins)
end

"Accumulates the negative local Hamiltonian gradient coming from the external field"
@inline function _accum_neggrad!(B::CUDA.CuArray{Vec3}, field::ExternalFieldCUDA)
    effBs = field.effBs
    nb = length(effBs)
    effBs = reshape(effBs, nb, ones(Int, ndims(B)-1)...)
    B .+= effBs
end
