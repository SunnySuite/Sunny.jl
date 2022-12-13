function get_intensity(sfd::SFData, q, ω::Float64;
    contraction = depolarize,
    c2q_temp = nothing,
)
    (; crystal, Δω, idx_info, site_infos) = sfd
    q = Vec3(q)

    # Just take nearest (put interpolation, etc. here)
    q, qi = nearest_q(sfd, q)
    ωi = round(Int, ω/Δω) + 1  # Could add validity check that this is within bounds

    data_point = raw_data_point(sfd, qi, ωi)
    elems = phase_averaged_elements(data_point, q, crystal, site_infos)
    intensity = contraction(elems, q, idx_info)
    if !isnothing(c2q_temp)
        intensity *= c2q(Δω*round(Int, ω/Δω), c2q_temp)
    end

    return abs(intensity)
end

get_intensity(sf::StructureFactor, q, ω::Float64; kwargs...) = get_intensity(sf.sfdata, q, ω; kwargs...)


## get_intensities with checks

# Approach: 
#    - Construct a list of qis from qs (Remove repeated or just do repeated calculations?)
#    - Extract all of these at once in a big slice
#    - Perform reductions and contractions efficiently
#    - Then interpolate resulting intensities

