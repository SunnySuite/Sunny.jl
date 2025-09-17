struct MSymOp
    R :: Mat3
    T :: Vec3
    p :: Int # Parity p == -1 indicates time reversal
end

function MSymOp(str::AbstractString)
    parts = split(str, ",")
    length(parts) == 4 || error("Magnetic symop must have 4 parts")
    (; R, T) = SymOp(join(parts[1:3], ","))
    p = parse(Int, parts[4])
    abs(p) == 1 || error("Time-reversal must be ±1")
    return MSymOp(R, T, p)
end

function MSymOp(s::SymOp)
    return MSymOp(s.R, s.T, +1)
end

function Base.show(io::IO, s::MSymOp)
    print(io, repr(MSymOp), "(\"")
    show(io, "text/plain", s)
    print(io, "\")")
end

function Base.show(io::IO, ::MIME"text/plain", s::MSymOp)
    (; R, T, p) = s
    show(io, MIME("text/plain"), SymOp(R, T))
    @assert abs(p) == 1
    print(io, p == -1 ? ",-1" : ",+1")
end

function Base.:*(s1::MSymOp, s2::MSymOp)
    return MSymOp(s1.R * s2.R, s1.T + s1.R * s2.T, s1.p*s2.p)
end

function Base.inv(s::MSymOp)
    Rinv = inv(s.R)
    MSymOp(Rinv, -Rinv*s.T, s.p)
end

function transform(s::MSymOp, r)
    return s.R * r + s.T
end

# The magnetic dipole is a pseudo-vector. It is invariant to space-inversion, so
# remove the determinant of s.R. The dipole does flip under time-reversal, which
# is marked by s.p = -1.
function transform_dipole(s::MSymOp, μ)
    return s.R * det(s.R) * s.p * μ
end
