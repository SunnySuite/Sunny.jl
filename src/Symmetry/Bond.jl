"""Defines the Bond type"""

"""
    Bond{D}

Represents a class of bond between pairs of sites in a `D`-dimensional crystal,
 separated from each other by a certain vector.
"""
struct Bond{D}
    "The sublattice index of the first site."
    i :: Int
    "The sublattice index of the second site."
    j :: Int
    "The displacement vector between sites, in units of lattice vectors."
    n :: SVector{D, Int}
end

function Bond(i, j, n)
    D = length(n)
    n = convert(SVector{D, Int}, n)
    Bond{D}(i, j, n)
end

"Represents a bond expressed as two fractional coordinates"
struct BondRaw
    r1::SVector{3, Float64}
    r2::SVector{3, Float64}
end

function Bond(cryst::Crystal, b::BondRaw)
    i1 = position_to_index(cryst, b.r1)
    i2 = position_to_index(cryst, b.r2)
    r1 = cryst.positions[i1]
    r2 = cryst.positions[i2]
    n = round.(Int, (b.r2-b.r1) - (r2-r1))
    return Bond{3}(i1, i2, n)
end

function BondRaw(cryst::Crystal, b::Bond{3})
    return BondRaw(cryst.positions[b.i], cryst.positions[b.j]+b.n)
end

function Base.show(io::IO, mime::MIME"text/plain", bond::Bond{3})
    print(io, "Bond($(bond.i), $(bond.j), $(bond.n))")
end

function position_to_index(cryst::Crystal, r::Vec3)
    return findfirst(r′ -> is_same_position(r, r′; symprec=cryst.symprec), cryst.positions)
end

function distance(cryst::Crystal, b::BondRaw)
    return norm(cryst.lat_vecs * (b.r1 - b.r2))
end

function distance(cryst::Crystal, b::Bond{3})
    return distance(cryst, BondRaw(cryst, b))
end

function transform(s::SymOp, b::BondRaw)
    return BondRaw(transform(s, b.r1), transform(s, b.r2))
end

function transform(cryst::Crystal, s::SymOp, b::Bond{3})
    return Bond(cryst, transform(s, BondRaw(cryst, b)))
end

function Base.reverse(b::BondRaw)
    return BondRaw(b.r2, b.r1)
end
