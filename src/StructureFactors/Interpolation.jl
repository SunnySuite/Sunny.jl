abstract type InterpolationScheme{NumInterp} end
struct NoInterp <: InterpolationScheme{1} end
struct LinearInterp <: InterpolationScheme{8} end

#=
ddtodo: Explanation of interpolation "API"
=# 

function stencil_intensities(sf::StructureFactor, ks, idcs, ω, iω, ::InterpolationScheme{NumInterp}, contraction::Contraction{T}, temp, ffdata) where {NumInterp, T}
    return SVector{NumInterp, T}(calc_intensity(sf, ks[n], idcs[n], ω, iω, contraction, temp, ffdata) for n in 1:NumInterp)
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


function stencil_points(sf::StructureFactor, q, ::NoInterp)
    Ls = sf.sftraj.sys.latsize 
    m = round.(Int, Ls .* q)
    im = map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex
    return (m,), (im,)
end


function stencil_points(sf::StructureFactor, q, ::LinearInterp)
    Ls = sf.sftraj.sys.latsize 
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
    ms = map(x -> x .+ base, offsets) 
    ims = map(ms) do m
        map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex
    end
    return ms, ims
end