abstract type InterpolationScheme{N} end
struct NoInterp <: InterpolationScheme{1} end
struct LinearInterp <: InterpolationScheme{8} end

function stencil_intensities(sf::StructureFactor, q, iq, ω, iω, ::InterpolationScheme{1}, contraction::Contraction{T}, temp, ffdata) where T
    return SVector{1, T}(calc_intensity(sf, q, iq, ω, iω, contraction, temp, ffdata))
end

function stencil_intensities(sf::StructureFactor, qs, iqs, ω, iω, ::InterpolationScheme{N}, contraction::Contraction{T}, temp, ffdata) where {N, T}
    return SVector{N, T}(calc_intensity(sf, qs[n], iqs[n], ω, iω, contraction, temp, ffdata) for n in 1:N)
end

function interpolated_intensity(::StructureFactor, _, _, stencil_intensities, ::NoInterp) 
    return only(stencil_intensities)
end

function interpolated_intensity(::StructureFactor, q_target, qs, stencil_intensities, ::LinearInterp) 
    (; q000, q111) = qs
    c000, c100, c010, c110, c001, c101, c011, c111 = stencil_intensities
    x, y, z = q_target
    x0, y0, z0 = q000
    x1, y1, z1 = q111

    xd = (x - x0)/(x1 - x0)
    yd = (y - y0)/(y1 - y0)
    zd = (z - z0)/(z1 - z0)

    c00 = c000*(1-xd) + c100*xd
    c01 = c001*(1-xd) + c101*xd
    c10 = c010*(1-xd) + c110*xd
    c11 = c011*(1-xd) + c111*xd

    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd

    return c0*(1-zd) + c1*zd
end


function stencil_qs(sfd::SFData, q, ::NoInterp)
    q = convert(Vec3, q)
    data = sfd.data
    Ls = (size(data, 2), size(data, 3), size(data, 4)) # Avoid slice allocation 
    ls = round.(Int, Ls .* q / 2π)
    q = @. 2π*ls/Ls
    # Below is a bit ugly, but couldn't figure out how to do it in another manner
    # without allocation (e.g., with a comprehension type syntax or map).
    qi = CartesianIndex{3}(  
        mod(ls[1], Ls[1]) + 1,
        mod(ls[2], Ls[2]) + 1,
        mod(ls[3], Ls[3]) + 1,
    ) 
    return q, qi 
end

# Simplify this
function stencil_qs(sfd::SFData, q, ::LinearInterp)
    q = convert(Vec3, q)
    data = sfd.data
    Ls = (size(data, 2), size(data, 3), size(data, 4)) 
    base = floor.(Int, Ls .* q / 2π)
    ls = (;
        l000 = base,
        l100 = base + SVector{3, Int64}(1, 0, 0),
        l010 = base + SVector{3, Int64}(0, 1, 0), 
        l110 = base + SVector{3, Int64}(1, 1, 0), 
        l001 = base + SVector{3, Int64}(0, 0, 1),
        l101 = base + SVector{3, Int64}(1, 0, 1), 
        l011 = base + SVector{3, Int64}(0, 1, 1),
        l111 = base + SVector{3, Int64}(1, 1, 1), 
    )
    qs = (;
        q000 = 2π .* ls.l000 ./ Ls,
        q100 = 2π .* ls.l100 ./ Ls,
        q010 = 2π .* ls.l010 ./ Ls,
        q110 = 2π .* ls.l110 ./ Ls,
        q001 = 2π .* ls.l001 ./ Ls,
        q101 = 2π .* ls.l101 ./ Ls,
        q011 = 2π .* ls.l011 ./ Ls,
        q111 = 2π .* ls.l111 ./ Ls,
    )
    qis = (;
        qi000 = CartesianIndex{3}(mod(ls.l000[1], Ls[1]) + 1,  mod(ls.l000[2], Ls[2]) + 1,  mod(ls.l000[3], Ls[3]) + 1),
        qi100 = CartesianIndex{3}(mod(ls.l100[1], Ls[1]) + 1,  mod(ls.l100[2], Ls[2]) + 1,  mod(ls.l100[3], Ls[3]) + 1),
        qi010 = CartesianIndex{3}(mod(ls.l010[1], Ls[1]) + 1,  mod(ls.l010[2], Ls[2]) + 1,  mod(ls.l010[3], Ls[3]) + 1),
        qi110 = CartesianIndex{3}(mod(ls.l110[1], Ls[1]) + 1,  mod(ls.l110[2], Ls[2]) + 1,  mod(ls.l110[3], Ls[3]) + 1),
        qi001 = CartesianIndex{3}(mod(ls.l001[1], Ls[1]) + 1,  mod(ls.l001[2], Ls[2]) + 1,  mod(ls.l001[3], Ls[3]) + 1),
        qi101 = CartesianIndex{3}(mod(ls.l101[1], Ls[1]) + 1,  mod(ls.l101[2], Ls[2]) + 1,  mod(ls.l101[3], Ls[3]) + 1),
        qi011 = CartesianIndex{3}(mod(ls.l011[1], Ls[1]) + 1,  mod(ls.l011[2], Ls[2]) + 1,  mod(ls.l011[3], Ls[3]) + 1),
        qi111 = CartesianIndex{3}(mod(ls.l111[1], Ls[1]) + 1,  mod(ls.l111[2], Ls[2]) + 1,  mod(ls.l111[3], Ls[3]) + 1),
    ) 
    return qs, qis
end