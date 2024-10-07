struct SiteSymmetry
    symbol       :: String
    multiplicity :: Int
    wyckoff      :: Char
end

function SiteSymmetry(sg_number::Int, r; symprec)
    g = find_sitegroup(sg_number, r; symprec)
    wp = Crystalline.position(g)
    mult = Crystalline.multiplicity(wp)
    letter = Crystalline.label(wp)[end]
    symbol = sitegroup_symbol(g)
    return SiteSymmetry(symbol, mult, letter)
end

# Does the position r belong to the Wyckoff position wp = (c, F). In other
# words, can F x + c = r be satisfied modulo 1 for some x = [α, β, γ]? Let F⁺
# denote the pseudo-inverse of F and b = r - c. Multiply both sides by F F⁺, and
# use the property F F⁺ F = F, to find the constraint F F⁺ b = b mod 1.
function belongs_to_wyckoff(r, wp::Crystalline.WyckoffPosition{3}; symprec=1e-8)
    (c, F) = Crystalline.parts(wp)
    b = c - r
    return all_integer(F * pinv(F) * b - b; symprec)
end

# Assumes standard spacegroup setting
function find_sitegroup(sg, r; symprec)
    D = 3
    for g in reverse(Crystalline.sitegroups(sg, D))
        if any(belongs_to_wyckoff.(Ref(r), Crystalline.orbit(g); symprec))
            return g
        end
    end
    @assert false
end

# The returned label is missing "dots". That is, it is not an "Oriented
# site-symmetry symbol", per Section 1.2.12 of International Tables for
# Crystallography (2006). Vol. E. [1]. TODO: To load site symmetry labels for
# all Hall numbers, we could parse the file spglib/database/Wyckoff.csv. An
# alternative is to iterate through all spacegroups 1..270 and call the PyXtal
# function `pyxtal.symmetry.Group(sgnum).Wyckoff_positions`.
#
# [1] https://onlinelibrary.wiley.com/iucr/itc/Ea/ch1o2v0001/sec1o2o12/.
function sitegroup_symbol(g)
    pg, _, _ = Crystalline.find_isomorphic_parent_pointgroup(g)
    return Crystalline.label(pg)
end
