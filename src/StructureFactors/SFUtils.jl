function c2q(ω, kT)
    (ω == 0.0) && (return 1.0)
    kT == 0.0 ? 1e-12 : kT
    return ω/(kT*(1 - exp(-ω/kT)))
end

function nearest_q(sfd::SFData, q)
    data = sfd.data
    Ls = (size(data, 2), size(data, 3), size(data, 4)) # Avoid slice allocation 
    ls = round.(Int, Ls .* q / 2π)
    q = @. 2π*ls/Ls
    qi = (mod(l, L)+1 for (l, L) in zip(ls, Ls))
    return q, qi
end

# Note that this is an inefficient slice. It would be nice to have the
# two atom indices immediately after the component index, but this
# conflicts with the efficiency of field calculations and integration.
function raw_data_point(sfd::SFData, qi, ωi)
    data = sfd.data
    nelems, ns = size(data, 1), size(data, 5)
    return SArray{Tuple{nelems, ns, ns}, ComplexF64, 3, ns*ns*nelems}(
        data[:, qi..., :, :, ωi]
    )
end