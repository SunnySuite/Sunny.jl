function qgrid(sf::StructureFactor; bzsize=(1,1,1))
    Ls = sf.sftraj.sys.size 
    offsets = map(L -> isodd(L) ? 1 : 0, Ls)
    up = Ls .* bzsize
    hi = map(L -> L - div(L, 2), up) .- offsets
    lo = map(L -> 1 - div(L, 2), up) .- offsets
    qs = zeros(Vec3, up...)
    for (k, lz) in enumerate(lo[3]:hi[3]), (j, ly) in enumerate(lo[2]:hi[2]), (i, lx) in enumerate(lo[1]:hi[1])
        qs[i,j,k] = Vec3(lx/Ls[1], ly/Ls[2], lz/Ls[3]) 
    end
    return qs
end

function ωvals(sf::StructureFactor)
    sfd = sf.sfdata
    Δω = sfd.Δω
    nω = size(sfd.data, 7)
    hω = div(nω, 2) + 1
    return collect(0:(hω-1)) .* Δω
end

function ωvals_all(sf::StructureFactor)
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

function classical_to_quantum(ω, kT)
    (ω == 0.0) && (return 1.0)
    kT == 0.0 ? 1e-12 : kT
    return ω/(kT*(1 - exp(-ω/kT)))
end
