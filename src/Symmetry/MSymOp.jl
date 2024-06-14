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



# This is a low-level function that will eventually be used to support the
# reading of MCIF files. All keyword parameters follow the MCIF format: their
# components should be interpreted in the coordinate system defined by the
# lattice vectors of the supercell.
function propagate_mcif_dipoles(sys; positions, dipoles, magn_operations, magn_centerings)
    orig_crystal(sys) == sys.crystal || error("System employs a reshaped chemical cell")

    # Defined such that r = D r′ where r is a "usual" position, while r′ is a
    # position given in multiples of the system's (super) lattice vectors
    D = Diagonal(Vec3(sys.latsize))

    # Convert positions, dipoles, and symops to multiples of lattice vectors for
    # the chemical cell
    positions = [D * r for r in positions]
    dipoles = [D * d for d in dipoles]
    magn_operations = [MSymOp(D * s.R * inv(D), D * s.T, s.p) for s in magn_operations]
    magn_centerings = [MSymOp(D * s.R * inv(D), D * s.T, s.p) for s in magn_centerings]

    # Use the zero vector as a marker for unvisited sites
    fill!(sys.dipoles, zero(Vec3))

    for (r, d) in zip(positions, dipoles)
        for s1 in magn_operations, s2 in magn_centerings
            s = s1 * s2
            # Transformed dipole. Because spin is a pseudo-vector, it is
            # invariant to an inversion in space. For this reason, remove the
            # determinant of R. If the parity p = ±1 is negative, this imposes
            # time reversal, which effectively flips the spin.
            d_new = (s.R * det(s.R) * s.p) * d
            r_new = s.R * r + s.T
            site = try 
                position_to_site(sys, r_new)
            catch _
                rethrow(ErrorException("MCIF position $r_new is missing in the chemical cell"))
            end

            state = dipolar_state(sys, site, d_new)
            if !iszero(sys.dipoles[site])
                sys.dipoles[site] ≈ state.s || error("Conflicting dipoles at site $site")
            end
            setspin!(sys, state, site)
        end
    end

    unvisited = [site for site in eachsite(sys) if iszero(sys.dipoles[site])]
    if !isempty(unvisited)
        error("Missing dipoles for sites $(collect(Tuple.(unvisited)))")
    end
end
