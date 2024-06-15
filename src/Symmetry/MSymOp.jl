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
    MSymOp(s1.R * s2.R, s1.T + s1.R * s2.T, s1.p*s2.p)
end
