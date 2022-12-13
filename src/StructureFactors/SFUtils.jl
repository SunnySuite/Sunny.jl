function c2q(ω, kT)
    (ω == 0.0) && (return 1.0)
    kT == 0.0 ? 1e-12 : kT
    return ω/(kT*(1 - exp(-ω/kT)))
end

function nearest_q(sfd::SFData, q)
    q = convert(Vec3, q)
    data = sfd.data
    Ls = (size(data, 2), size(data, 3), size(data, 4)) # Avoid slice allocation 
    ls = round.(Int, Ls .* q / 2π)
    q = @. 2π*ls/Ls
    # Below is a bit ugly, but couldn't figure out how to do it in another manner
    # without allocation.
    f(l, L) = mod(l, L) + 1
    qi = CartesianIndex{3}(  
        mod(ls[1], Ls[1]) + 1,
        mod(ls[2], Ls[2]) + 1,
        mod(ls[3], Ls[3]) + 1,
    ) 
    return q, qi 
end

function qvals(sf::StructureFactor)
    data = sf.sfdata.data
    Ls = (size(data, 2), size(data, 3), size(data, 4))
    hLs = map(L -> div(L, 2) + 1, Ls)
    qs = zeros(Vec3, Ls)
    for i in CartesianIndices(Ls) 
        vals = @. 2π * (i.I - 1) / Ls
        qs[i] = map((val, l, hL) -> l <= hL ? val : val - 2π, vals, i.I, hLs)
    end
    return qs
end

function ωvals(sf::StructureFactor)
    sfd = sf.sfdata
    Δω = sfd.Δω
    nω = size(sfd.data, 7)
    hω = div(nω, 2) + 1
    ωs = collect(0:(nω-1)) .* Δω
    for i ∈ hω+1:nω
        ωs[i] -= 2ωs[hω]
    end
    return ωs
end


# Note that this is an inefficient slice. It would be nice to have the
# two atom indices immediately after the component index, but this
# conflicts with the efficiency of field calculations and integration.
function raw_data_point(sfd::SFData, qi, ωi)
    data = sfd.data
    nelems, ns = size(data, 1), size(data, 5)
    return SArray{Tuple{nelems, ns, ns}, ComplexF64, 3, ns*ns*nelems}(
        data[:, qi, :, :, ωi]
    )
end