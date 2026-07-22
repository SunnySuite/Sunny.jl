################################################################################
# Types
################################################################################

# One part of an entangled unit.
struct UnitPart
    atom  :: Int            # Atom index in the uncontracted crystal
    Δcell :: SVector{3,Int} # Cell of this atom relative to the unit centroid's cell
end

# Metadata for a System with entangled units. Hilbert dimension N of the
# uncontracted System{N} is intentionally omitted to avoid dynamic lookup.
struct Entanglement <: AbstractEntanglement
    units           :: Vector{Vector{UnitPart}}         # Grouping of units into parts
    uncontracted    :: System                           # System prior to entanglement
    lifted_spin_ops :: Vector{SVector{3, HermitianC64}} # Product-space spin ops for atoms
end


################################################################################
# Helper functions
################################################################################

# True if the system carries entanglement metadata.
is_entangled(sys::System) = !isnothing(sys.entanglement)

# Concrete-typed view of the entanglement metadata. The explicit type ascription
# avoids JET warnings.
@inline get_entanglement(sys::System) = sys.entanglement :: Entanglement

function uncontracted_system(sys::System)
    is_entangled(sys) ? get_entanglement(sys).uncontracted : sys
end

function nbaresites(sys::System)
    length(eachsite(uncontracted_system(sys)))
end

# Returns the (u, p) indices inverse to units[u][p].atom
function unit_and_part(units, atom)
    for (u, parts) in enumerate(units)
        p = findfirst(part -> part.atom == atom, parts)
        isnothing(p) || return (u, p)
    end
    error("Atom $atom is not assigned to any unit")
end

# Bare atom indices comprising a unit, ordered to match the units entries.
function atoms_in_unit(units, u)
    [atom for (; atom) in units[u]]
end
function atoms_in_unit(sys::System, u)
    is_entangled(sys) ? atoms_in_unit(get_entanglement(sys).units, u) : [u]
end

# Given a site representing an entangled unit, return an iterator over the
# corresponding "bare sites" of the uncontracted system.
function bare_sites_for_unit(ent::Entanglement, unit_site)
    (; uncontracted, units) = ent
    (; dims) = uncontracted
    parts = units[to_atom(unit_site)]
    return Iterators.map(parts) do (; atom, Δcell)
        # Mirrors `bonded_site`
        CartesianIndex(altmod1.(to_cell(unit_site) .+ Δcell, dims)..., atom)
    end
end

function lifted_spin_op(sys::System, i)
    if is_entangled(sys)
        return get_entanglement(sys).lifted_spin_ops[i]
    else
        return spin_matrices_of_dim(; N=sys.Ns[i])
    end
end

# Needed because uncontracted.dipoles and related state are mutably updated
function clone_entanglement(ent::Entanglement)
    return Entanglement(ent.units, clone_system(ent.uncontracted),
                        ent.lifted_spin_ops)
end
clone_entanglement(::Nothing) = nothing

# Conversion between user-supplied `groupings` and regularized `units`
function groupings_to_units(groupings)
    map(groupings) do grouping
        map(grouping) do elem
            elem isa Integer ? UnitPart(elem, zero(Vec3)) : UnitPart(elem...)
        end
    end
end
function units_to_groupings(units)
    map(units) do parts
        [(atom, Vector(Δcell)) for (; atom, Δcell) in parts]
    end
end


################################################################################
# Pair-coupling contraction
################################################################################

# Convert a bond between physical atoms into the corresponding bond between
# contracted units. If `bond.n` points from atom i in one cell to atom j in a
# neighboring cell, then the unit-to-unit cell displacement must be corrected by
# the two atoms' cell offsets inside their units.
function contracted_bond(bond::Bond, units)
    ui, pi = unit_and_part(units, bond.i)
    uj, pj = unit_and_part(units, bond.j)
    n = SVector{3,Int}(bond.n)
    Δcell_i = units[ui][pi].Δcell
    Δcell_j = units[uj][pj].Δcell
    n_new = n + Δcell_i - Δcell_j
    return Bond(ui, uj, n_new)
end

# Checks whether a bond given in terms of the original crystal lies inside a
# single unit of the contracted system.
function bond_is_in_unit(bond::Bond, units)
    newbond = contracted_bond(bond, units)
    return newbond.i == newbond.j && all(iszero, newbond.n)
end

# Converts what was a pair coupling between two different sites of a single unit
# in the original system into an "onsite coupling" for the entangled unit.
# Accumulates to `op` in place.
function accum_intra_unit_pair_coupling!(op, pc, uncontracted_Ns, units)
    pc.isculled && return

    ui, pi = unit_and_part(units, pc.bond.i)
    uj, pj = unit_and_part(units, pc.bond.j)
    @assert ui == uj

    Ns_unit = uncontracted_Ns[atoms_in_unit(units, ui)]
    Ni, Nj = Ns_unit[[pi, pj]]

    # Add scalar part
    op .+= pc.scalar * I(size(op, 1))

    # Add the bilinear, biquadratic, and general parts, each decomposed into a
    # sum of tensor product pairs a ⊗ b. Both sites lie in the same unit, so both
    # factors embed into the same product space.
    for (a, b) in pair_coupling_tensor_data(pc, Ni, Nj)
        A = local_op_to_product_space(a, pi, Ns_unit)
        B = local_op_to_product_space(b, pj, Ns_unit)
        op .+= A * B
    end

    return
end

# Promote an atomistic-level pair coupling to an inter-unit pair coupling on the
# contracted system.
function promote_inter_unit_pair_coupling(pc, uncontracted_Ns, units)
    (; i, j) = pc.bond
    ui, pi = unit_and_part(units, i)
    uj, pj = unit_and_part(units, j)

    Ns1 = uncontracted_Ns[atoms_in_unit(units, ui)]
    Ns2 = uncontracted_Ns[atoms_in_unit(units, uj)]
    Ni = uncontracted_Ns[1, 1, 1, i]
    Nj = uncontracted_Ns[1, 1, 1, j]

    # Decompose the pair operator as ∑ₖ aₖ ⊗ bₖ. Then promote each factor to the
    # tensor product space, aₖ → Aₖ.
    data = map(pair_coupling_tensor_data(pc, Ni, Nj)) do (a, b)
        A = Hermitian(local_op_to_product_space(Matrix(a), pi, Ns1))
        B = Hermitian(local_op_to_product_space(Matrix(b), pj, Ns2))
        (A, B)
    end

    newbond = contracted_bond(pc.bond, units)
    gen1 = spin_matrices_of_dim(; N=prod(Ns1))
    gen2 = spin_matrices_of_dim(; N=prod(Ns2))
    return (; newbond, pc.scalar, tensordec=TensorDecomposition(gen1, gen2, data))
end

# Map ModelParam of an uncontracted system to that of an entangled system.
function contract_param(param::ModelParam, uncontracted_Ns, units)
    nunits = length(units)

    # Intra-unit terms become onsite (on-bond) operators, one per touched unit.
    onsites = Tuple{Int, OnsiteCoupling}[]
    for u in 1:nunits
        Ns = uncontracted_Ns[atoms_in_unit(units, u)]
        unit_operator = zeros(ComplexF64, prod(Ns), prod(Ns))
        relevant = atoms_in_unit(units, u)

        # Bare onsite couplings embedded into the unit's product space.
        for (atom, oc) in param.onsites
            atom in relevant || continue
            _, p = unit_and_part(units, atom)
            unit_operator += local_op_to_product_space(oc, p, Ns)
        end

        # Intra-unit pair couplings fold into the same onsite operator.
        for pc in param.pairs
            u′, _ = unit_and_part(units, pc.bond.i)
            if bond_is_in_unit(pc.bond, units) && u′ == u
                accum_intra_unit_pair_coupling!(unit_operator, pc, uncontracted_Ns, units)
            end
        end

        iszero(unit_operator) || push!(onsites, (u, Hermitian(unit_operator)))
    end

    # Inter-unit pair couplings become pair couplings between contracted units.
    pairs = PairCoupling[]
    for pc in param.pairs
        bond_is_in_unit(pc.bond, units) && continue
        (; newbond, scalar, tensordec) = promote_inter_unit_pair_coupling(pc, uncontracted_Ns, units)
        bilin = biquad = 0.0 # Dipoles are NaN for entangled systems
        push!(pairs, PairCoupling(newbond, scalar, bilin, biquad, tensordec))
    end

    return ModelParam(param.label, param.val, onsites, pairs)
end


################################################################################
# Symmetry analysis
################################################################################

# A canonical key for a unit's parts. Two units share a key iff they are related
# by lattice translation, independent of how their parts are permuted.
function unit_position_key(parts)
    sorted = sort([(atom, Δcell...) for (; atom, Δcell) in parts])
    (_, o1, o2, o3) = sorted[1]
    return [(a, x-o1, y-o2, z-o3) for (a, x, y, z) in sorted]
end

# The k-atom analog of `transform(cryst, s, bond)`. Applies symop `s` to each
# part and reindexes to (atom, cell). The order of the returned parts records
# the permutation induced by `s`. This ordering would be needed for future
# functionality that transforms product-space operators.
function transform_unit(cryst::Crystal, s::SymOp, parts)
    map(parts) do (; atom, Δcell)
        a, cell = position_to_atom_and_offset(cryst, transform(s, cryst.positions[atom] + Δcell))
        UnitPart(a, SVector{3,Int}(cell))
    end
end

# True if every symop maps the set of units onto itself (each unit's image is
# another unit in the set, up to lattice translation and part reordering).
function are_units_symmetry_consistent(crystal, units)
    keys = Set(unit_position_key(parts) for parts in units)
    for parts in units, s in crystal.sg.symops
        unit_position_key(transform_unit(crystal, s, parts)) in keys || return false
    end
    return true
end

# Shift each unit's cell offsets by a common cell so its centroid lands in
# [0, 1) (tolerantly: -ϵ maps to 0).
function center_units(cryst::Crystal, units)
    map(units) do parts
        centroid = sum(cryst.positions[atom] + Δcell for (; atom, Δcell) in parts) / length(parts)
        cell = round.(Int, centroid - wrap_to_unit_cell(centroid))
        [UnitPart(atom, Δcell - cell) for (; atom, Δcell) in parts]
    end
end

# Reassign a unit's cell offsets to greedily compactify its atoms (existing
# offsets are ignored and recomputed from scratch).
function compactify_cell_offsets(cryst::Crystal, parts)
    atoms = [atom for (; atom) in parts]
    # Place atoms one-by-one. The periodic image of each new atom is selected by
    # minimizing distance to an already-placed atom.
    n = length(atoms)
    Δcells = Vector{SVector{3,Int}}(undef, n)
    Δcells[1] = zero(SVector{3,Int})
    placed = [1]
    for _ in 2:n
        best = (dist=Inf, k=0, Δ=zero(SVector{3,Int}))
        for i in placed, k in 1:n
            k in placed && continue
            ri = cryst.positions[atoms[i]] + Δcells[i]
            Δ = round.(Int, ri - cryst.positions[atoms[k]])
            d = norm(cryst.latvecs * (cryst.positions[atoms[k]] + Δ - ri))
            d < best.dist - 1e-10 && (best = (dist=d, k=k, Δ=SVector{3,Int}(Δ)))
        end
        Δcells[best.k] = best.Δ
        push!(placed, best.k)
    end
    return [UnitPart(a, Δ) for (a, Δ) in zip(atoms, Δcells)]
end

# Heuristic attempt to rescue a failed grouping by finding a new grouping with
# compactified cell offsets. Returns `nothing` if the proposed grouping would be
# symmetry-inconsistent.
function compact_grouping_suggestion(crystal, units)
    compacted = [compactify_cell_offsets(crystal, parts) for parts in units]
    # Reject if the recompactified offsets still break symmetry
    if are_units_symmetry_consistent(crystal, compacted)
        return units_to_groupings(center_units(crystal, compacted))
    else
        return nothing
    end
end

# Validate the entangled units in stages for better error messages.
function validate_units(crystal, units, enforce_symmetry)
    unit_atoms = [atoms_in_unit(units, u) for u in eachindex(units)]

    # Every atom must be assigned to a unit
    atoms_flat = collect(Iterators.flatten(unit_atoms))
    atoms_missing = setdiff(1:natoms(crystal), atoms_flat)
    isempty(atoms_missing) || error("Atoms $atoms_missing have not been assigned to a unit")

    # Atoms cannot be repeated
    is_repeated(a) = count(==(a), atoms_flat) > 1
    atoms_repeated = filter(is_repeated, unique(atoms_flat))
    isempty(atoms_repeated) || error("Atoms $atoms_repeated are repeatedly assigned")

    # The partition must be consistent with the crystal symmetry
    if enforce_symmetry
        # Atom indices within each unit must be invariant under symops
        unit_sets = Set(sort(atoms) for atoms in unit_atoms)
        for s in crystal.sg.symops, atoms in unit_atoms
            image = sort([position_to_atom_and_offset(crystal, transform(s, crystal.positions[a]))[1] for a in atoms])
            image in unit_sets || error("Groupings would be split by spacegroup symmetries; pass `enforce_symmetry=false` to override")
        end

        # Every unit must transform into another unit under symops. This is a
        # more stringent check, because it also depends on cell offsets.
        if !are_units_symmetry_consistent(crystal, units)
            suggestion = compact_grouping_suggestion(crystal, units)
            if isnothing(suggestion)
                error("Symmetry-inconsistency detected; consider different cell offsets or pass `enforce_symmetry=false` to override")
            else
                error("Offsets are symmetry-inconsistent; consider groupings $suggestion")
            end
        end
    end
end

# Maps the provided crystal to a new "contracted" one with a site per entangled
# unit, positioned at each unit's centroid.
function contract_crystal(crystal, units)
    # Unit centroids (coordinates in [0, 1) by construction)
    new_positions = map(units) do parts
        sum(crystal.positions[atom] + Δcell for (; atom, Δcell) in parts) / length(parts)
    end

    # Build the contracted spacegroup, dropping the symops to mark unavailability
    # of symmetry analysis (in analogy to `reshape_crystal`).
    if are_units_symmetry_consistent(crystal, units)
        # Determine equivalence classes from unit centroids
        new_classes = zeros(Int, length(new_positions))
        for u in eachindex(new_positions)
            if iszero(new_classes[u])
                c = maximum(new_classes) + 1
                for s in crystal.sg.symops
                    r = transform(s, new_positions[u])
                    v = findfirst(p -> is_periodic_copy(r, p), new_positions)
                    new_classes[v] = c
                end
            end
        end
        # Although symmetry is unbroken, we nonetheless empty the symops because
        # symmetry analysis of entangled interactions is not yet implemented.
        new_sg = Spacegroup(SymOp[], crystal.sg.label, crystal.sg.number, crystal.sg.setting)
    else
        # For broken symmetries, revert to P1 with no symops.
        new_classes = collect(eachindex(new_positions))
        new_sg = Spacegroup(SymOp[], spacegroup_label(1), 1, one(SymOp))
    end

    new_types = fill("", length(new_positions))
    return Crystal(nothing, crystal.latvecs, crystal.recipvecs, new_positions,
                   new_types, new_classes, new_sg)
end


################################################################################
# Entanglement construction
################################################################################

# Populate the entanglement metadata of a (possibly reshaped) contracted `sys`,
# given the corresponding `uncontracted` system. The `orig_units` defines the
# entanged units in the context of the "original" chemical crystal.
function rebuild_entanglement!(sys::System, uncontracted, orig_units)
    nunits = natoms(sys.crystal)
    nparts = length(first(orig_units))
    natoms_bare = natoms(uncontracted.crystal)
    @assert natoms_bare == nunits * nparts

    # Build the units (unit → parts) for the possibly reshaped sys.crystal
    units = [Vector{UnitPart}(undef, nparts) for _ in 1:nunits]
    for a in 1:natoms_bare
        orig_a = map_atom_to_other_crystal(uncontracted.crystal, a, orig_crystal(uncontracted))
        orig_u, p = unit_and_part(orig_units, orig_a)

        # Global displacement from the unit centroid to the part's atom
        orig_cryst = orig_crystal(sys)
        orig_Δpos = orig_crystal(uncontracted).positions[orig_a] +
                    orig_units[orig_u][p].Δcell - orig_cryst.positions[orig_u]
        global_Δpos = orig_cryst.latvecs * orig_Δpos

        # Global position of the unit centroid
        global_pos = global_position_at(uncontracted, (1, 1, 1, a)) - global_Δpos
        # New unit index and the cell hosting it
        u, unit_cell = position_to_atom_and_offset(sys.crystal, sys.crystal.latvecs \ global_pos)
        # Cell displacement from this unit to an atom in cell (0, 0, 0)
        Δcell = -unit_cell

        units[u][p] = UnitPart(a, Δcell)
    end

    lifted_spin_ops = map(1:natoms_bare) do atom
        S_local = spin_matrices_of_dim(; N=uncontracted.Ns[1, 1, 1, atom])
        u, p = unit_and_part(units, atom)
        Ns_unit = uncontracted.Ns[atoms_in_unit(units, u)]
        ntuple(α -> local_op_to_product_space(S_local[α], p, Ns_unit), 3)
    end

    sys.entanglement = Entanglement(units, uncontracted, lifted_spin_ops)
    return nothing
end

# Copy per-site state from an uncontracted system into the entangled system.
function transfer_uncontracted_state!(sys::System, uncontracted::System)
    ent = get_entanglement(sys)
    dst = ent.uncontracted  # entanglement's own clone, distinct from arg
    @. dst.dipoles   = uncontracted.dipoles
    @. dst.coherents = uncontracted.coherents
    @. dst.extfield  = uncontracted.extfield
    @. dst.κs        = uncontracted.κs
    @. dst.gs        = uncontracted.gs
    for unit in eachsite(sys)
        Zs = (uncontracted.coherents[bs] for bs in bare_sites_for_unit(ent, unit))
        sys.coherents[unit] = kron(Zs...)
    end
end

"""
    entangle_system(sys::System, groupings; enforce_symmetry=true)

Create a new [`System`](@ref) with `groupings` of atoms contracted into
"entangled units". The input `sys` should be in `:SUN` mode with interactions
already populated. The returned system will contain these same interactions in
tensor-product form. This formalism is most useful when the couplings within one
unit are relatively strong compared to the couplings between different units.

Each element of `groupings` contains a list of atomic parts. For example,
`groupings = [[1, 2], [3, 4]]` would return a system in which atoms `[1, 2]` and
atoms `[3, 4]` are dimerized. If any entangled unit straddles the cell boundary,
then all atoms must be paired with cell offsets. For example, replacing `[1, 2]`
with `[(1, [0, 0, 0]), (2, [-1, 0, 0])]` would indicate that atom `2`
participates in the dimer via its periodic image in the neighboring cell
``-𝐚_1``.

By default, the grouping of atoms into units must respect all spacegroup
symmetries. Select `enforce_symmetry=false` to disable this symmetry check.
"""
function entangle_system(uncontracted::System{M}, groupings::Vector{<: Vector{<: Union{Integer, Tuple{Integer, Vector{<: Integer}}}}}; enforce_symmetry=true) where M
    is_homogeneous(uncontracted) || error("Entanglement not supported on inhomogeneous systems")
    is_entangled(uncontracted)  && error("Cannot entangle a system twice")

    # If `uncontracted` has already been reshaped, then entangle on the original
    # unreshaped system and repeat the reshaping.
    if !isnothing(uncontracted.origin)
        shape = cell_shape(uncontracted) * diagm(Vec3(uncontracted.dims))
        sys = reshape_supercell(entangle_system(uncontracted.origin, groupings; enforce_symmetry), shape)
        transfer_uncontracted_state!(sys, uncontracted)
        return sys
    end

    # Convert to regularized format
    units = groupings_to_units(groupings)

    # Check symmetry consistency
    validate_units(uncontracted.crystal, units, enforce_symmetry)

    # Build contracted crystal with units centered on [0, 1)
    units = center_units(uncontracted.crystal, units)
    contracted_crystal = contract_crystal(uncontracted.crystal, units)

    # Local Hilbert space dimension per unit. TODO: Relax uniformity constraint.
    nunits = natoms(contracted_crystal)
    Ns_contracted = [prod(uncontracted.Ns[atoms_in_unit(units, u)]) for u in 1:nunits]
    @assert allequal(Ns_contracted) "The dimensions of the local contracted Hilbert spaces must be uniform"
    N = first(Ns_contracted)

    # Build the contracted System directly (cf. `reshape_supercell_aux`). Dipole
    # and external field data are undefined here (NaN), and will instead be
    # derived from the associated uncontracted system.
    dims      = uncontracted.dims
    Ns        = fill(N, dims..., nunits)
    κs        = ones(Float64, dims..., nunits)
    gs        = fill(Mat3(I), dims..., nunits)
    ints      = empty_interactions(:SUN, nunits, N)
    extfield  = fill(Vec3(NaN, NaN, NaN), dims..., nunits)
    dipoles   = fill(Vec3(NaN, NaN, NaN), dims..., nunits)
    coherents = zeros(CVec{N}, dims..., nunits)
    sys = System(nothing, :SUN, contracted_crystal, dims, Ns, κs, gs, ModelParam[],
                 Symbol[], ints, nothing, extfield, dipoles, coherents, Array{Vec3, 4}[],
                 Array{CVec{N}, 4}[], copy(uncontracted.rng), nothing)

    # Contract each labeled parameter independently (contraction is linear in
    # coupling strength), then let the ordinary machinery fill the interactions.
    sys.params = [contract_param(p, uncontracted.Ns, units) for p in uncontracted.params]
    repopulate_couplings_from_params!(sys)

    # Set entanglement metadata in sys. Keep a clone of the uncontracted system,
    # e.g., as a helper for dipole-level physics.
    uncontracted = clone_system(uncontracted)
    rebuild_entanglement!(sys, uncontracted, units)

    # Map uncontracted spin states to the product-space representations in sys
    transfer_uncontracted_state!(sys, uncontracted)

    return sys
end

function entangle_system(::System{0}, _)
    error("Cannot entangle units of a `:dipole`-mode `System`. Use `:SUN` mode.")
end


################################################################################
# Coherent states to bare dipoles
################################################################################

# ⟨Z|S|Z⟩ for an embedded spin-operator triple S (an entry of
# `lifted_spin_ops`).
expected_dipole(Z, S) = Vec3(real(dot(Z, S[1], Z)), real(dot(Z, S[2], Z)), real(dot(Z, S[3], Z)))

# Helper for set_energy_grad_coherents! that accumulates Zeeman and Ewald
# interactions. These need special handling for entangled systems because they
# refer to the uncontracted ("bare") dipoles.
function accum_bare_field_grad_coherents!(HZ, Z::Array{CVec{N}, 4}, ent::Entanglement) where N
    (; uncontracted, lifted_spin_ops) = ent

    # Reuse the bare system's persistent scratch buffers.
    S_bare, dE_dS_bare = get_dipole_buffers(uncontracted, 2)
    fill!(dE_dS_bare, zero(Vec3))

    # Expected dipoles of the uncontracted system for the provided coherents Z.
    for u in CartesianIndices(Z)
        for bs in bare_sites_for_unit(ent, u)
            S_bare[bs] = expected_dipole(Z[u], lifted_spin_ops[to_atom(bs)])
        end
    end

    # Energy gradient dE/dS associated with Zeeman and Ewald couplings, for each
    # bare atom of the uncontracted system.
    for site in eachsite(uncontracted)
        dE_dS_bare[site] += uncontracted.gs[site]' * uncontracted.extfield[site]
    end
    if !isnothing(uncontracted.ewald)
        accum_ewald_grad!(dE_dS_bare, S_bare, uncontracted)
    end

    # Accumulate gradients for entangled units u
    for u in CartesianIndices(Z)
        for bs in bare_sites_for_unit(ent, u)
            S = lifted_spin_ops[to_atom(bs)]
            h = dE_dS_bare[bs]
            for β in 1:3
                HZ[u] += h[β] * mul_svec(S[β], Z[u])
            end
        end
    end
    return
end

# Sync the bare dipoles with the coherent state of an unit at `site`
function sync_bare_dipoles_at!(sys::System, site)
    Z = sys.coherents[site]
    ent = get_entanglement(sys)
    (; uncontracted, lifted_spin_ops) = ent
    for bs in bare_sites_for_unit(ent, site)
        uncontracted.dipoles[bs] = expected_dipole(Z, lifted_spin_ops[to_atom(bs)])
    end
    return
end


################################################################################
# Measurements
################################################################################

# Transform an atom-level MeasureSpec (for `sys.entanglement.uncontracted`) into
# a unit-level MeasureSpec (for the contracted `sys`).
function entangled_measure(measure, sys::System)
    @assert is_entangled(sys)                 # System is entangled
    @assert num_parts_per_unit(measure) == 1  # Measure is for uncontracted system

    (; uncontracted, units) = get_entanglement(sys)

    nobs = num_observables(measure)
    dims = sys.dims
    nunits = length(units)
    nparts = length(units[1])  # uniform by construction

    Op = eltype(measure.observables)
    new_obs     = Array{Op, 6}(undef, nobs, dims..., nunits, nparts)
    new_offsets = zeros(Vec3, nunits, nparts)
    new_ff      = Array{FormFactor, 3}(undef, nobs, nunits, nparts)

    for u in 1:nunits
        Ns_unit = uncontracted.Ns[atoms_in_unit(units, u)]
        for (p, (; atom, Δcell)) in enumerate(units[u])
            # Displacement from the unit centroid to the part's atom
            new_offsets[u, p] = uncontracted.crystal.positions[atom] + Δcell - sys.crystal.positions[u]
            for μ in 1:nobs
                new_ff[μ, u, p] = measure.formfactors[μ, atom, 1]
                for c in CartesianIndices(dims)
                    c′ = altmod1.(Tuple(c) .+ Δcell, dims)
                    A = measure.observables[μ, c′..., atom, 1]
                    A_product = local_op_to_product_space(A, p, Ns_unit)
                    new_obs[μ, c, u, p] = Hermitian(A_product)
                end
            end
        end
    end

    return MeasureSpec(new_obs, measure.corr_pairs, measure.combiner, new_ff; offsets=new_offsets)
end
