abstract type InterpolationScheme{NumInterp} end
struct NoInterp <: InterpolationScheme{1} end
struct LinearInterp <: InterpolationScheme{8} end

#=
ddtodo: Explanation of interpolation "API"
=# 

function stencil_intensities(sf::StructureFactor, ks, idcs, ω, iω, ::InterpolationScheme{NumInterp}, contraction::Contraction{T}, temp, ffdata, ::Val{NCorr}, ::Val{NAtoms}) where {NumInterp, T, NCorr, NAtoms}
    return SVector{NumInterp, T}(calc_intensity(sf, ks[n], idcs[n], ω, iω, contraction, temp, ffdata, Val(NCorr), Val(NAtoms)) for n in 1:NumInterp)
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

    # Each of the following lines causes a 32 byte allocation
    Ls = size(sf.samplebuf)[2:4] 
    m = round.(Int, Ls .* q)
    im = map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex{3}

    ## The following lines cause no allocations, but don't seem to be any faster.
    #     _, L1, L2, L3, _, _ = size(sf.samplebuf)
    #     m = (round(Int, L1*q[1]), round(Int, L2*q[2]), round(Int, L3*q[3]))
    #     im = CartesianIndex{3}(mod(m[1], L1)+1, mod(m[2], L2)+1, mod(m[3], L3)+1)

    return (m,), (im,)
end


function stencil_points(sf::StructureFactor, q, ::LinearInterp)
    Ls = size(sf.samplebuf)[2:4] 
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
        map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)) |> CartesianIndex{3}
    end
    return ms, ims
end