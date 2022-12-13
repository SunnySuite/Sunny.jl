export trace, depolarize

# TODO: Add warning
function trace(elems::SVector{N, ComplexF64}, idx_info::SortedDict) where N
    intensity = 0.0im
    for ((i,j), idx) in idx_info 
        if i == j
            intensity += elems[idx]
        end
    end
    return intensity
end
trace(elems, _, idx_info) = trace(elems, idx_info)

# Make sure in dipole mode
function depolarize(elems::SVector{N, ComplexF64}, q::Vec3, idx_info::SortedDict) where N
    q /= norm(q) + 1e-12
    dip_factor = SMatrix{3, 3, Float64, 9}(I(3) - q * q')
    intensity = 0.0
    for ((α, β), idx) in idx_info # Loop from 1 to 6 
        factor = α == β ? 1.0 : 2.0 # Double off-diagonal contribution (if ij is in iteration, ji will not be)
        intensity += factor * dip_factor[α, β] * real(elems[idx])  
    end
    return intensity
end
