abstract type InterpolationScheme{N} end
struct NoInterp <: InterpolationScheme{1} end
struct LinearInterp <: InterpolationScheme{8} end

function stencil_intensities(sf::StructureFactor, qs, iqs, ω, iω, ::InterpolationScheme{N}, contraction::Contraction{T}, temp, ffdata) where {N, T}
    return SVector{N, T}(calc_intensity(sf, qs[n], iqs[n], ω, iω, contraction, temp, ffdata) for n in 1:N)
end

function interpolated_intensity(::StructureFactor, _, _, stencil_intensities, ::NoInterp) 
    return only(stencil_intensities)
end

function interpolated_intensity(::StructureFactor, q_target, qs, stencil_intensities, ::LinearInterp) 
    q000,    _,    _,    _,    _,    _,    _, q111 = qs 
    c000, c100, c010, c110, c001, c101, c011, c111 = stencil_intensities
    x, y, z = q_target
    x0, y0, z0 = q000
    x1, y1, z1 = q111

    xd = (x-x0)/(x1-x0)
    yd = (y-y0)/(y1-y0)
    zd = (z-z0)/(z1-z0)

    c00 = c000*(1-xd) + c100*xd
    c01 = c001*(1-xd) + c101*xd
    c10 = c010*(1-xd) + c110*xd
    c11 = c011*(1-xd) + c111*xd

    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd

    return c0*(1-zd) + c1*zd
end


function stencil_qs(sfd::SFData, q, ::NoInterp)
    Ls = (size(sfd.data, 2), size(sfd.data, 3), size(sfd.data, 4)) 
    l = round.(Int, Ls .* q)
    q = l ./ Ls
    qi = map(i -> mod(l[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex
    return (q,), (qi,)
end


function stencil_qs(sfd::SFData, q, ::LinearInterp)
    Ls = (size(sfd.data, 2), size(sfd.data, 3), size(sfd.data, 4))
    base = map(x -> floor(Int64, x), Ls .* q) 
    offsets = (
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 0), 
        (1, 1, 0), 
        (0, 0, 1),
        (1, 0, 1), 
        (0, 1, 1),
        (1, 1, 1)
    )
    ls = map(x -> x .+ base, offsets) 
    qs = map(x -> x ./ Ls, ls)
    qis = map(ls) do l
        map(i -> mod(l[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex
    end
    return qs, qis
end