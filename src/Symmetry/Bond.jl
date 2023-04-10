"""
    Bond(i, j, n)

Represents a bond between atom indices `i` and `j`. `n` is a vector of three integers
specifying unit cell displacement in terms of lattice vectors.
"""
struct Bond
    i :: Int
    j :: Int
    n :: SVector{3, Int}
end

# A bond expressed as two positions in fractional coordinates
struct BondPos
    ri::Vec3
    rj::Vec3
end

function Bond(cryst::Crystal, b::BondPos)
    i = position_to_atom(cryst, b.ri)::Int
    j = position_to_atom(cryst, b.rj)::Int
    ri = cryst.positions[i]
    rj = cryst.positions[j]
    n = round.(Int, (b.rj-b.ri) - (rj-ri))
    return Bond(i, j, n)
end

function BondPos(cryst::Crystal, b::Bond)
    return BondPos(cryst.positions[b.i], cryst.positions[b.j]+b.n)
end

function Base.show(io::IO, ::MIME"text/plain", bond::Bond)
    print(io, "Bond($(bond.i), $(bond.j), $(bond.n))")
end

# The displacement vector ``𝐫_j - 𝐫_i`` in global coordinates between atoms
# `b.i` and `b.j`, accounting for the integer offsets `b.n` between unit cells.
function global_displacement(cryst::Crystal, b::BondPos)
    return cryst.latvecs * (b.rj - b.ri)
end

function global_displacement(cryst::Crystal, b::Bond)
    return global_displacement(cryst, BondPos(cryst, b))
end

# The global distance between atoms in bond `b`. Equivalent to
# `norm(global_displacement(cryst, b))`.
function global_distance(cryst::Crystal, b::BondPos)
    return norm(global_displacement(cryst, b))
end

function global_distance(cryst::Crystal, b::Bond)
    return norm(global_displacement(cryst, b))
end

function transform(s::SymOp, b::BondPos)
    return BondPos(transform(s, b.ri), transform(s, b.rj))
end

function transform(cryst::Crystal, s::SymOp, b::Bond)
    return Bond(cryst, transform(s, BondPos(cryst, b)))
end

function Base.reverse(b::Bond)
    return Bond(b.j, b.i, -1 .* b.n)
end

function Base.reverse(b::BondPos)
    return BondPos(b.rj, b.ri)
end
