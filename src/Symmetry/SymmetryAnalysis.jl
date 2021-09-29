function is_equivalent_by_translation(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    # Displacements between the two bonds
    D1 = b2.r1 - b1.r1
    D2 = b2.r2 - b1.r2
    # Round components of D1 to nearest integers
    n = round.(D1, RoundNearest)
    # If both n ≈ D1 and n ≈ D2, then the bonds are equivalent by translation
    return norm(n - D1) < cryst.symprec && norm(n - D2) < cryst.symprec
end


function is_equivalent_by_symmetry(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    return !isempty(symmetries_between_bonds(cryst::Crystal, b1::BondRaw, b2::BondRaw))
end

function is_equivalent_by_symmetry(cryst::Crystal, b1::Bond{3}, b2::Bond{3})
    return is_equivalent_by_symmetry(cryst, BondRaw(cryst, b1), BondRaw(cryst, b2))
end

# Generate list of all symmetries that transform b2 into b1, along with parity
function symmetries_between_bonds(cryst::Crystal, b1::BondRaw, b2::BondRaw)
    # Fail early if two bonds describe different real-space distances
    # (dimensionless error tolerance is measured relative to the minimum lattice
    # constant ℓ)
    if b1 != b2
        ℓ = minimum(norm, eachcol(cryst.lat_vecs))
        d1 = distance(cryst, b1) / ℓ
        d2 = distance(cryst, b2) / ℓ
        if abs(d1-d2) > cryst.symprec
            return Tuple{SymOp, Bool}[]
        end
    end

    function bonds_equiv(s)
        b2′ = transform(s, b2)
        if is_equivalent_by_translation(cryst, b1, b2′)
            (s, true)
        elseif is_equivalent_by_translation(cryst, b1, reverse(b2′))
            (s, false)
        else
            nothing
        end
    end

    ret = (bonds_equiv(s) for s in cryst.symops)
    return Iterators.filter(!isnothing, ret)
end

"Returns all bonds in `cryst` for which `bond.i == i`"
function all_bonds_for_atom(cryst::Crystal, i::Int, max_dist; min_dist=0.0)
    # be a little generous with the minimum and maximum distances
    ℓ = minimum(norm, eachcol(cryst.lat_vecs))
    max_dist += 4 * cryst.symprec * ℓ
    min_dist -= 4 * cryst.symprec * ℓ

    # columns are the reciprocal vectors
    recip_vecs = 2π * inv(cryst.lat_vecs)'

    # box_lengths[i] represents the perpendicular distance between two parallel
    # boundary planes spanned by lattice vectors a_j and a_k (where indices j
    # and k differ from i)
    box_lengths = [a⋅b/norm(b) for (a,b) = zip(eachcol(cryst.lat_vecs), eachcol(recip_vecs))]
    n_max = round.(Int, max_dist ./ box_lengths, RoundUp)

    bonds = Bond{3}[]

    # loop over neighboring cells
    for n1 in -n_max[1]:n_max[1]
        for n2 in -n_max[2]:n_max[2]
            for n3 in -n_max[3]:n_max[3]
                n = SVector(n1, n2, n3)
                
                # loop over all atoms within neighboring cell
                for j in eachindex(cryst.positions)
                    b = Bond{3}(i, j, n)
                    if min_dist <= distance(cryst, b) <= max_dist
                        push!(bonds, b)
                    end
                end
            end
        end
    end

    return bonds
end


# Calculate score for a bond. Lower would be preferred.
function _score_bond(cryst::Crystal, b)
    # Favor bonds with fewer nonzero elements in basis matrices J
    Js = basis_for_symmetry_allowed_couplings(cryst, b)
    nnz = [count(abs.(J) .> 1e-12) for J in Js]
    score = sum(nnz)

    # Favor bonds with smaller unit cell displacements. Positive
    # displacements are slightly favored over negative displacements.
    # Displacements in x are slightly favored over y, etc.
    score += norm((b.n .- 0.1) .* [0.07, 0.08, 0.09])

    # Favor smaller indices and indices where i < j
    score += 1e-2 * (b.i + b.j) + 1e-2 * (b.i < b.j ? -1 : +1)

    return score
end


function bond_multiplicity(cryst::Crystal, b::Bond{3})
    return length(all_symmetry_related_bonds_for_atom(cryst, b.i, b))
end


"Produces a list of 'canonical' bonds that belong to different symmetry equivalence classes."
function canonical_bonds(cryst::Crystal, max_dist)
    # Atom indices, one for each equivalence class
    canon_atoms = [findfirst(isequal(c), cryst.classes) for c in unique(cryst.classes)]
    
    # Bonds, one for each equivalent class
    cbonds = Bond[]
    for i in canon_atoms
        for b in all_bonds_for_atom(cryst, i, max_dist)
            if !any(is_equivalent_by_symmetry(cryst, b, b′) for b′ in cbonds)
                push!(cbonds, b)
            end
        end
    end

    # Sort by distance
    sort!(cbonds, by=b->distance(cryst, b))

    # Replace each canonical bond by the "best" equivalent bond
    return map(cbonds) do cb
        # Find full set of symmetry equivalent bonds
        equiv_bonds = unique([transform(cryst, s, cb) for s in cryst.symops])
        # Take the bond with lowest score
        scores = [_score_bond(cryst, b) for b in equiv_bonds]
        return equiv_bonds[findmin(scores)[2]]
    end
end

# For each element of `bonds`, return an index into `canonical_bonds` that gives
# the equivalent bond
function equivalent_bond_indices(cryst::Crystal, canonical_bonds, bonds)
    map(bonds) do b
        findfirst(canonical_bonds) do cb
             is_equivalent_by_symmetry(cryst, b, cb)
        end
    end
end

function all_symmetry_related_bonds_for_atom(cryst::Crystal, i::Int, b_ref::Bond{3})
    bs = Bond{3}[]
    dist = distance(cryst, b_ref)
    for b in all_bonds_for_atom(cryst, i, dist; min_dist=dist)
        if is_equivalent_by_symmetry(cryst, b_ref, b)
            push!(bs, b)
        end
    end
    return bs
end

"""
    all_symmetry_related_bonds(cryst::Crystal, b_ref::Bond{3}) :: Vector{Bond{3}}

Construct a list of all bonds which are symmetry-equivalent to the reference bond `b_ref`
 within `cryst.
"""
function all_symmetry_related_bonds(cryst::Crystal, b_ref::Bond{3})
    bs = Bond{3}[]
    for i in eachindex(cryst.positions)
        append!(bs, all_symmetry_related_bonds_for_atom(cryst, i, b_ref))
    end
    bs
end
