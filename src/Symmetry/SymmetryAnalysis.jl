# Wrap x into the range [0,1). To account for finite precision, wrap 1-ϵ to -ϵ.
function wrap_to_unit_cell(x::Float64; tol=1e-12)
    return mod(x+tol, 1) - tol
end

function wrap_to_unit_cell(r::Vec3; tol=1e-12)
    return wrap_to_unit_cell.(r; tol)
end

function is_periodic_copy(r1::Vec3, r2::Vec3; tol=1e-12)
    return all_integer(r1-r2; tol)
end

function is_periodic_copy(b1::BondPos, b2::BondPos; tol=1e-12)
    # Displacements between the two bonds
    D1 = b2.ri - b1.ri
    D2 = b2.rj - b1.rj
    # Round components of D1 to nearest integers
    n = round.(D1, RoundNearest)
    # If both n ≈ D1 and n ≈ D2, then the bonds are equivalent by translation
    return norm(n - D1) < tol && norm(n - D2) < tol
end

function position_to_atom(cryst::Crystal, r::Vec3; tol=1e-12)
    return findfirst(r′ -> is_periodic_copy(r, r′; tol), cryst.positions)
end

function position_to_atom_and_offset(cryst::Crystal, r::Vec3; tol=1e-12)
    i = position_to_atom(cryst, r; tol)
    isnothing(i) && error("Position $r not found in crystal")

    offset = round.(Int, r - wrap_to_unit_cell(r))
    @assert isapprox(cryst.positions[i]+offset, r; atol=tol)
    return (i, offset)
end


# Generate list of SymOps for the pointgroup of atom i
function symmetries_for_pointgroup_of_atom(cryst::Crystal, i::Int)
    ret = SymOp[]
    r = cryst.positions[i]
    for s in cryst.sg.symops
        r′ = transform(s, r)
        if is_periodic_copy(r, r′)
            push!(ret, s)
        end
    end
    return ret
end


# General a list of all symmetries that transform i2 into i1. (Convention for
# definition of `s` is consistent with symmetries_between_bonds())
function symmetries_between_atoms(cryst::Crystal, i1::Int, i2::Int)
    validate_symops(cryst)

    ret = SymOp[]
    r1 = cryst.positions[i1]
    r2 = cryst.positions[i2]
    for s in cryst.sg.symops
        if is_periodic_copy(r1, transform(s, r2))
            push!(ret, s)
        end
    end
    return ret
end


# The list of atoms symmetry-equivalent to i_ref
function all_symmetry_related_atoms(cryst::Crystal, i_ref::Int)
    # The result is the set of atoms sharing the symmetry "class"
    c = cryst.classes[i_ref]
    ret = findall(==(c), cryst.classes)

    # Calculate the result another way, as a consistency check.
    equiv_atoms = Int[]
    r_ref = cryst.positions[i_ref]
    for s in cryst.sg.symops
        push!(equiv_atoms, position_to_atom(cryst, transform(s, r_ref)))
    end
    @assert sort(unique(equiv_atoms)) == ret

    return ret
end


# Generate list of all symmetries that transform b2 into b1, along with parity
function symmetries_between_bonds(cryst::Crystal, b1::BondPos, b2::BondPos)
    validate_symops(cryst)

    # Fail early if two bonds describe different real-space distances
    if b1 != b2
        d1 = global_distance(cryst, b1)
        d2 = global_distance(cryst, b2)
        atol = 1e-12 * opnorm(cryst.latvecs)
        if !isapprox(d1, d2; atol)
            return Tuple{SymOp, Bool}[]
        end
    end

    ret = Tuple{SymOp, Bool}[]
    for s in cryst.sg.symops
        b2′ = transform(s, b2)
        if is_periodic_copy(b1, b2′)
            push!(ret, (s, true))
        elseif is_periodic_copy(b1, reverse(b2′))
            push!(ret, (s, false))
        end
    end
    return ret
end

# Are i1 and i2 symmetry-equivalent sites?
function is_related_by_symmetry(cryst::Crystal, i1::Int, i2::Int)
    return cryst.classes[i1] == cryst.classes[i2]
end

# Is there a symmetry operation that transforms `b1` into either `b2` or its
# reverse?
function is_related_by_symmetry(cryst::Crystal, b1::Bond, b2::Bond)
    return !isempty(symmetries_between_bonds(cryst, BondPos(cryst, b1), BondPos(cryst, b2)))
end

# Collect a list of all periodic images of the positions `rs` that are within
# some distance of a given `pt`. The return values `(idxs, offsets)` will be the
# complete lists of indices and offsets that satisfy,
# 
# r = rs[idxs[α]]
# n = offsets[α]
# dist = norm(latvecs * (r + n - pt))
# @assert min_dist ≤ dist ≤ max_dist

function all_offsets_within_distance(latvecs, rs, pt; min_dist=0, max_dist, nonzeropart=false)
    # box_lengths[i] represents the perpendicular distance between two parallel
    # boundary planes spanned by lattice vectors a_j and a_k (where indices j
    # and k differ from i)
    box_lengths = [a⋅b/norm(b) for (a,b) in zip(eachcol(latvecs), eachrow(inv(latvecs)))]
    n_max = round.(Int, max_dist ./ box_lengths, RoundUp)

    idxs = Int[]
    offsets = Vec3[]

    for (i, r) in enumerate(rs)
        for n1 in -n_max[1]:n_max[1], n2 in -n_max[2]:n_max[2], n3 in -n_max[3]:n_max[3]
            n = Vec3(n1, n2, n3)
            nonzeropart && iszero(n) && continue

            dist = norm(latvecs * (r + n - pt))
            if min_dist <= dist <= max_dist
                push!(idxs, i)
                push!(offsets, n)
            end
        end
    end

    return (idxs, offsets)
end

# Returns all bonds in `cryst` for which `bond.i == i`
function all_bonds_for_atom(cryst::Crystal, i::Int, max_dist; min_dist=0.0)
    atol = 1e-12 * opnorm(cryst.latvecs)
    max_dist += atol
    min_dist -= atol

    idxs, offsets = all_offsets_within_distance(cryst.latvecs, cryst.positions, cryst.positions[i]; min_dist, max_dist)

    return map(zip(idxs, offsets)) do (j, n)
        Bond(i, j, n)
    end
end


# Calculate score for a bond. Lower would be preferred.
function score_bond(cryst::Crystal, b::Bond)
    # Favor bonds with fewer nonzero elements in basis matrices J
    Js = basis_for_symmetry_allowed_couplings(cryst, b)
    nnz = [count(abs.(J) .> 1e-12) for J in Js]
    score = Float64(sum(nnz))

    # Favor bonds with smaller unit cell displacements. Positive
    # displacements are slightly favored over negative displacements.
    # Displacements in x are slightly favored over y, etc.
    score += norm((b.n .- 0.1) .* [0.07, 0.08, 0.09])

    # Favor smaller indices and indices where i < j
    score += 1e-2 * (b.i + b.j) + 1e-2 * (b.i < b.j ? -1 : +1)

    return score
end

# Indices of the unique elements in `a`, ordered by their first appearance.
function unique_indices(a)
    map(x->x[1], unique(x->x[2], enumerate(a)))
end

"""    reference_bonds(cryst::Crystal, max_dist)

Returns a list of [`Bond`](@ref)s, one for each symmetry equivalence class, up
to the `max_dist` cutoff in length units. These reference bonds are
heuristically selected to simplify the expression of symmetry-allowed
interactions."""
function reference_bonds(cryst::Crystal, max_dist::Float64; min_dist=0.0)
    # Bonds, one for each equivalence class
    ref_bonds = Bond[]
    for i in unique_indices(cryst.classes)
        for b in all_bonds_for_atom(cryst, i, max_dist; min_dist)
            if !any(is_related_by_symmetry(cryst, b, b′) for b′ in ref_bonds)
                push!(ref_bonds, b)
            end
        end
    end

    # Sort by distance
    sort!(ref_bonds, by=b->global_distance(cryst, b))

    # Replace each canonical bond by the "best" equivalent bond
    return map(ref_bonds) do rb
        # Find all symmetry equivalent bonds
        equiv_bonds = [transform(cryst, s, rb) for s in cryst.sg.symops]
        # Include also reverse bonds
        equiv_bonds = vcat(equiv_bonds, reverse.(equiv_bonds))
        # Take the bond with lowest score
        return argmin(b -> score_bond(cryst, b), unique(equiv_bonds))
    end
end
reference_bonds(cryst::Crystal, max_dist) = reference_bonds(cryst, convert(Float64, max_dist))

"""
    all_symmetry_related_bonds_for_atom(cryst::Crystal, i::Int, b::Bond)

Returns a list of all bonds that start at atom `i`, and that are symmetry
equivalent to bond `b` or its reverse.
"""
function all_symmetry_related_bonds_for_atom(cryst::Crystal, i::Int, b_ref::Bond)
    bs = Bond[]
    dist = global_distance(cryst, b_ref)
    for b in all_bonds_for_atom(cryst, i, dist; min_dist=dist)
        if is_related_by_symmetry(cryst, b_ref, b)
            push!(bs, b)
        end
    end
    return bs
end

"""
    all_symmetry_related_bonds(cryst::Crystal, b::Bond)

Returns a list of all bonds that are symmetry-equivalent to bond `b` or its
reverse.
"""
function all_symmetry_related_bonds(cryst::Crystal, b_ref::Bond)
    bs = Bond[]
    for i in eachindex(cryst.positions)
        append!(bs, all_symmetry_related_bonds_for_atom(cryst, i, b_ref))
    end
    bs
end

"""    coordination_number(cryst::Crystal, i::Int, b::Bond)

Returns the number times that atom `i` participates in a bond equivalent to `b`.
In other words, the count of bonds that begin at atom `i` and that are
symmetry-equivalent to `b` or its reverse.

Defined as `length(all_symmetry_related_bonds_for_atom(cryst, i, b))`.
"""
function coordination_number(cryst::Crystal, i::Int, b::Bond)
    return length(all_symmetry_related_bonds_for_atom(cryst, i, b))
end
