################################################################################
# Types
################################################################################

# Index to a part of an entangled unit.
struct UnitAndPartIndex
    unit :: Int                   # Unit index in contracted crystal
    part :: Int                   # Part index within unit
end

# A member atom within an entangled unit. In principle, Δcell is redundant -- it
# could be recalculated from Δpos and the unit centers.
struct UnitMember
    atom  :: Int            # Atom index in the associated crystal
    Δpos  :: Vec3           # Fractional position relative to unit center
    Δcell :: SVector{3,Int} # Cell of this atom relative to the unit's cell
end

# Bidirectional mapping between entangled units and their atomic members
struct UnitMap
    atom_to_unit    :: Vector{UnitAndPartIndex}   # bare atom → (unit, part)
    unit_to_members :: Vector{Vector{UnitMember}} # unit → bare atom offsets, ordered by part
end

# Metadata for a System with entangled units. Hilbert dimension N of the
# uncontracted System{N} is intentionally omitted to avoid dynamic lookup.
struct Entanglement <: AbstractEntanglement
    uncontracted          :: System                          # System prior to entanglement
    unit_map              :: UnitMap                         # Map on the reshaped uncontracted crystal
    bare_dipole_operators :: Vector{NTuple{3, HermitianC64}} # Product-space spin ops per physical atom
end


################################################################################
# Helper functions
################################################################################

# True if the system carries entanglement metadata.
is_entangled(sys::System) = !isnothing(sys.entanglement)

# Concrete-typed view of the entanglement metadata. The explicit type ascription
# avoids JET warnings.
@inline get_entanglement(sys::System) = sys.entanglement :: Entanglement

# Atom indices comprising a unit, ordered to match the unit_to_members field of
# a UnitMap.
atoms_in_unit(unit_map::UnitMap, i) = [m.atom for m in unit_map.unit_to_members[i]]

# UnitMember data (Δpos and Δcell) for a bare atom.
function unit_offsets_for_atom(unit_map::UnitMap, atom)
    (; unit, part) = unit_map.atom_to_unit[atom]
    return unit_map.unit_to_members[unit][part]
end

# Given a site representing an entangled unit, obtain the corresponding "bare
# site" for the uncontracted system. Implementation mirrors `bonded_site`.
function bare_site_for_unit_member(unit_site, member::UnitMember, dims)
    (; atom, Δcell) = member
    CartesianIndex(altmod1.(to_cell(unit_site) .+ Δcell, dims)..., atom)
end

# Given a site representing an entangled unit, return an iterator over the
# corresponding "bare sites" of the uncontracted system.
function bare_sites_for_unit(ent::Entanglement, unit_site)
    (; uncontracted, unit_map) = ent
    members = unit_map.unit_to_members[to_atom(unit_site)]
    return (bare_site_for_unit_member(unit_site, m, uncontracted.dims) for m in members)
end

# Needed because uncontracted.dipoles and related state are mutably updated
function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.uncontracted), ent.unit_map,
                        ent.bare_dipole_operators)
end
clone_entanglement(::Nothing) = nothing


################################################################################
# Pair-coupling contraction
################################################################################

# Convert a bond between physical atoms into the corresponding bond between
# contracted units. If `bond.n` points from atom i in one cell to atom j in a
# neighboring cell, then the unit-to-unit cell displacement must be corrected by
# the two atoms' cell offsets inside their units.
function contracted_bond(bond::Bond, unit_map::UnitMap)
    ui = unit_map.atom_to_unit[bond.i]
    uj = unit_map.atom_to_unit[bond.j]
    n = SVector{3,Int}(bond.n)
    Δcell_i = unit_offsets_for_atom(unit_map, bond.i).Δcell
    Δcell_j = unit_offsets_for_atom(unit_map, bond.j).Δcell
    n_new = n + Δcell_i - Δcell_j
    return Bond(ui.unit, uj.unit, n_new)
end

# Checks whether a bond given in terms of the original crystal lies inside a
# single unit of the contracted system.
function bond_is_in_unit(bond::Bond, unit_map::UnitMap)
    newbond = contracted_bond(bond, unit_map)
    return newbond.i == newbond.j && all(iszero, newbond.n)
end

# Converts what was a pair coupling between two different sites of a single unit
# in the original system into an "onsite coupling" for the entangled unit.
# Accumulates to `op` in place.
function accum_intra_unit_pair_coupling!(op, pc, uncontracted_Ns, unit_map::UnitMap)
    pc.isculled && return

    ui = unit_map.atom_to_unit[pc.bond.i]
    uj = unit_map.atom_to_unit[pc.bond.j]
    Ns_unit = uncontracted_Ns[atoms_in_unit(unit_map, ui.unit)]
    Ni, Nj = Ns_unit[[ui.part, uj.part]]

    # Add scalar part
    op .+= pc.scalar * I(size(op, 1))

    # Add the bilinear, biquadratic, and general parts, each decomposed into a
    # sum of tensor product pairs a ⊗ b. Both sites lie in the same unit, so both
    # factors embed into the same product space.
    for (a, b) in pair_coupling_tensor_data(pc, Ni, Nj)
        A = local_op_to_product_space(a, ui.part, Ns_unit)
        B = local_op_to_product_space(b, uj.part, Ns_unit)
        op .+= A * B
    end

    return
end

# Promote an atomistic-level pair coupling to an inter-unit pair coupling on the
# contracted system.
function promote_inter_unit_pair_coupling(pc, uncontracted_Ns, unit_map::UnitMap)
    (; i, j) = pc.bond
    ui = unit_map.atom_to_unit[i]
    uj = unit_map.atom_to_unit[j]

    Ns1 = uncontracted_Ns[atoms_in_unit(unit_map, ui.unit)]
    Ns2 = uncontracted_Ns[atoms_in_unit(unit_map, uj.unit)]
    Ni = uncontracted_Ns[1, 1, 1, i]
    Nj = uncontracted_Ns[1, 1, 1, j]

    # Decompose the pair operator as ∑ₖ aₖ ⊗ bₖ. Then promote each factor to the
    # tensor product space, aₖ → Aₖ.
    data = map(pair_coupling_tensor_data(pc, Ni, Nj)) do (a, b)
        A = Hermitian(local_op_to_product_space(Matrix(a), ui.part, Ns1))
        B = Hermitian(local_op_to_product_space(Matrix(b), uj.part, Ns2))
        (A, B)
    end

    newbond = contracted_bond(pc.bond, unit_map)
    gen1 = spin_matrices_of_dim(; N=prod(Ns1))
    gen2 = spin_matrices_of_dim(; N=prod(Ns2))
    return (; newbond, pc.scalar, tensordec=TensorDecomposition(gen1, gen2, data))
end

# Contract ModelParam data. Contraction is linear in the coupling strength, so
# each labeled parameter is contracted independently at unit strength; the
# ordinary `repopulate_couplings_from_params!` then rescales by `param.val`.
function contract_param(param::ModelParam, uncontracted_Ns, unit_map::UnitMap, Ns_units)
    nunits = length(unit_map.unit_to_members)

    # Intra-unit terms become onsite (on-bond) operators, one per touched unit.
    onsites = Tuple{Int, OnsiteCoupling}[]
    for u in 1:nunits
        Ns = Ns_units[u]
        unit_operator = zeros(ComplexF64, prod(Ns), prod(Ns))
        relevant = atoms_in_unit(unit_map, u)

        # Bare onsite couplings embedded into the unit's product space.
        for (atom, oc) in param.onsites
            atom in relevant || continue
            (; part) = unit_map.atom_to_unit[atom]
            unit_operator += local_op_to_product_space(oc, part, Ns)
        end

        # Intra-unit pair couplings fold into the same onsite operator.
        for pc in param.pairs
            if bond_is_in_unit(pc.bond, unit_map) && unit_map.atom_to_unit[pc.bond.i].unit == u
                accum_intra_unit_pair_coupling!(unit_operator, pc, uncontracted_Ns, unit_map)
            end
        end

        iszero(unit_operator) || push!(onsites, (u, Hermitian(unit_operator)))
    end

    # Inter-unit pair couplings become pair couplings between contracted units.
    pairs = PairCoupling[]
    for pc in param.pairs
        bond_is_in_unit(pc.bond, unit_map) && continue
        (; newbond, scalar, tensordec) = promote_inter_unit_pair_coupling(pc, uncontracted_Ns, unit_map)
        bilin = biquad = 0.0 # Dipoles are NaN for entangled systems
        push!(pairs, PairCoupling(newbond, scalar, bilin, biquad, tensordec))
    end

    return ModelParam(param.label, param.val, onsites, pairs)
end


################################################################################
# Entanglement construction
################################################################################

# Maps the provided crystal to a new "contracted" one. `unit_atoms[u]` and
# `unit_Δcells[u]` are atom indices and cell offsets for entangled unit `u`,
# indexed for the original chemical cell.
function contract_crystal(crystal, unit_atoms, unit_Δcells)
    @assert sum(length, unit_atoms; init=0) == natoms(crystal) "Invalid entangled unit specification."

    atom_to_unit = Vector{UnitAndPartIndex}(undef, natoms(crystal))
    unit_to_members = Vector{UnitMember}[]
    new_positions = Vec3[]

    for (unit_idx, (atoms, Δcells)) in enumerate(zip(unit_atoms, unit_Δcells))
        # Position of the atomic parts
        parts_pos = [crystal.positions[atom] + Δcell for (atom, Δcell) in zip(atoms, Δcells)]
        # Unit position is the centroid of parts positions
        unit_pos = sum(parts_pos) / length(parts_pos)
        # Cell of the unit, computed tolerantly (-ϵ maps to 0)
        unit_cell = round.(Int, unit_pos - wrap_to_unit_cell(unit_pos))
        # Unit position wrapped to [0, 1]
        wrapped_unit_pos = unit_pos - Vec3(unit_cell...)
        push!(new_positions, wrapped_unit_pos)

        # Record each atom's part, its cell offset, and its position offset
        members = map(enumerate(zip(atoms, Δcells, parts_pos))) do (part, (atom, Δcell, pos))
            atom_to_unit[atom] = UnitAndPartIndex(unit_idx, part)
            UnitMember(atom, pos - unit_pos, Δcell - unit_cell)
        end
        push!(unit_to_members, members)
    end

    # Space group = 1 allows "invalid" anisotropies from pair → onsite conversions
    new_crystal = Crystal(crystal.latvecs, new_positions, 1)
    unit_map = UnitMap(atom_to_unit, unit_to_members)

    return new_crystal, unit_map
end

# Populate the entanglement metadata of a (possibly reshaped) contracted `sys`,
# given the corresponding `uncontracted` system. The `orig_unit_map` defines the
# entanged units in the context of the "original" chemical crystal.
function rebuild_entanglement!(sys::System, uncontracted, orig_unit_map::UnitMap)
    nunits = natoms(sys.crystal)
    atoms_per_unit = length(first(orig_unit_map.unit_to_members))
    natoms_bare = natoms(uncontracted.crystal)
    @assert natoms_bare == nunits * atoms_per_unit

    # Build a unit_map for the possibly reshaped sys.crystal
    atom_to_unit = UnitAndPartIndex[]
    unit_to_members = [Vector{UnitMember}(undef, atoms_per_unit) for _ in 1:nunits]
    for a in 1:natoms_bare
        orig_a = map_atom_to_other_crystal(uncontracted.crystal, a, orig_crystal(uncontracted))
        part = orig_unit_map.atom_to_unit[orig_a].part

        # Fractional displacement from the unit center to member atom
        orig_Δpos = unit_offsets_for_atom(orig_unit_map, orig_a).Δpos
        global_Δpos = orig_crystal(sys).latvecs * orig_Δpos
        Δpos = uncontracted.crystal.latvecs \ global_Δpos

        # Cell offset from the unit center to member atom at cell [0, 0, 0]
        unit_center = global_position_at(uncontracted, (1, 1, 1, a)) - global_Δpos
        u, unit_cell = position_to_atom_and_offset(sys.crystal, sys.crystal.latvecs \ unit_center)
        Δcell =  Vec3(0, 0, 0) - unit_cell

        # Store offset data
        unit_to_members[u][part] = UnitMember(a, Δpos, Δcell)

        # Store unit part data
        push!(atom_to_unit, UnitAndPartIndex(u, part))
    end
    unit_map = UnitMap(atom_to_unit, unit_to_members)

    bare_dipole_operators = map(1:natoms_bare) do atom
        S_local = spin_matrices_of_dim(; N=uncontracted.Ns[1, 1, 1, atom])
        (; unit, part) = unit_map.atom_to_unit[atom]
        Ns_unit = uncontracted.Ns[atoms_in_unit(unit_map, unit)]
        ntuple(α -> local_op_to_product_space(S_local[α], part, Ns_unit), 3)
    end

    sys.entanglement = Entanglement(uncontracted, unit_map, bare_dipole_operators)
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
    entangle_system(sys::System{N}, groupings)

Create a new [`System`](@ref) with `groupings` of atoms contracted into
"entangled units". This formalism becomes especially useful when the exchange
interactions within a unit are relatively strong compared to the exchange
between distinct units.

Often, `groupings` will include just atom indices for the chemical cell. For
example, selecting `groupings = [[1, 2], [3, 4]]` would return a system in which
atoms `[1, 2]` and atoms `[3, 4]` are dimerized. If any entangled unit straddles
the cell boundary, however, then cell offsets must be provided for all atoms.
For example, replacing `[1, 2]` with `[(1, [0, 0, 0]), (2, [-1, 0, 0])]` would
indicate that atom `2` participates in the dimer via its periodic image in the
neighboring cell ``-𝐚_1``.

The input system should be in `:SUN` mode with its interactions already
populated. These interactions will be transferred to the entangled system in
tensor-product form.
"""
function entangle_system(uncontracted::System{M}, groupings::Vector{<: Vector{<: Union{Integer, Tuple{Integer, Vector{<: Integer}}}}}) where M
    is_homogeneous(uncontracted) || error("Entanglement not supported on inhomogeneous systems")
    is_entangled(uncontracted)  && error("Cannot entangle a system twice")

    # If `uncontracted` has already been reshaped, then entangle on the original
    # (unreshaped) system, and then re-perform the reshaping. 
    if !isnothing(uncontracted.origin)
        shape = cell_shape(uncontracted) * diagm(Vec3(uncontracted.dims))
        sys = reshape_supercell(entangle_system(uncontracted.origin, groupings), shape)
        transfer_uncontracted_state!(sys, uncontracted)
        return sys
    end

    # Preprocess `groupings` into parallel per-unit lists of atom indices and
    # cell offsets. If omitted, cell offsets default to zero.
    grouping_type = eltype(eltype(groupings))
    if grouping_type <: Integer
        unit_atoms = [collect(g) for g in groupings]
        unit_Δcells = [[zero(SVector{3,Int}) for _ in g] for g in groupings]
    else
        @assert grouping_type <: Tuple{Integer, Vector{<: Integer}}
        unit_atoms = [[Int(x) for (x, _) in g] for g in groupings]
        unit_Δcells = [[SVector{3,Int}(y) for (_, y) in g] for g in groupings]
    end

    # Construct contracted (P1) crystal and the chemical-cell map
    contracted_crystal, unit_map = contract_crystal(uncontracted.crystal, unit_atoms, unit_Δcells)

    # Local Hilbert space dimensions per unit (all must be equal). (TODO:
    # Determine if alternative behavior preferable in mixed case.)
    Ns_units = [uncontracted.Ns[atoms_in_unit(unit_map, u)] for u in eachindex(unit_map.unit_to_members)]
    Ns_contracted = map(prod, Ns_units)
    @assert allequal(Ns_contracted) "After contraction, the dimensions of the local Hilbert spaces on each generalized site must all be equal."
    N = first(Ns_contracted)

    # Build the contracted System directly (cf. `reshape_supercell_aux`). Dipole
    # and external field data are undefined here (NaN), and will instead be
    # derived from the associated uncontracted system.
    nunits    = natoms(contracted_crystal)
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
    sys.params = [contract_param(p, uncontracted.Ns, unit_map, Ns_units) for p in uncontracted.params]
    repopulate_couplings_from_params!(sys)

    # Set entanglement metadata in sys. Keep a clone of the uncontracted system,
    # e.g., as a helper for dipole-level physics.
    uncontracted = clone_system(uncontracted)
    rebuild_entanglement!(sys, uncontracted, unit_map)

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

# ⟨Z|S|Z⟩ for an embedded spin-operator triple S (an entry of `bare_dipole_operators`).
@inline bare_expected_dipole(Z, S) = Vec3(real(dot(Z, S[1], Z)), real(dot(Z, S[2], Z)), real(dot(Z, S[3], Z)))

# Translate the physical (bare-system) dipole-level fields into the entangled
# unit coherent gradient:
#
#   dE/dZ̄_u += Σ_{a ∈ u} Σ_β (dE/dSₐ)_β Sₐ^{β,embed} Z_u,
#
# where dE/dSₐ = gₐ' B (Zeeman) + (long-range dipole-dipole field), both
# produced per physical atom on the bare system. A unit's atoms feel independent
# fields, so this cannot be represented by a single per-unit dE/dS vector —
# hence the bare-atom granularity here rather than the contracted
# `set_energy_grad_dipoles!` path. `N` is the local dimension of entangled units
# (not of the `uncontracted`).
function accum_bare_field_grad_coherents!(HZ, Z::Array{CVec{N}, 4}, ent::Entanglement) where N
    (; uncontracted, bare_dipole_operators) = ent

    # Reuse the bare system's persistent scratch buffers.
    S_bare, dE_dS_bare = get_dipole_buffers(uncontracted, 2)
    fill!(dE_dS_bare, zero(Vec3))

    # Physical dipoles ⟨Sₐ⟩ evaluated from the *passed* `Z` (which may be a
    # trial state not yet synced into `uncontracted.dipoles`), so the
    # dipole-dependent Ewald field is consistent with the state whose gradient
    # is requested.
    for u in CartesianIndices(Z)
        Zu = Z[u]
        for bs in bare_sites_for_unit(ent, u)
            S_bare[bs] = bare_expected_dipole(Zu, bare_dipole_operators[to_atom(bs)])
        end
    end

    # Per-atom dE/dSₐ on the bare system: Zeeman gₐ' B plus, if enabled, the
    # long-range dipole-dipole field (which itself includes the g-tensor).
    for site in eachsite(uncontracted)
        dE_dS_bare[site] += uncontracted.gs[site]' * uncontracted.extfield[site]
    end
    if !isnothing(uncontracted.ewald)
        accum_ewald_grad!(dE_dS_bare, S_bare, uncontracted)
    end

    # Sum each physical atom's contribution into its containing unit's coherent
    # gradient.
    for u in CartesianIndices(Z)
        Zu = Z[u]
        for bs in bare_sites_for_unit(ent, u)
            S = bare_dipole_operators[to_atom(bs)]
            h = dE_dS_bare[bs]
            for β in 1:3
                HZ[u] += h[β] * mul_svec(S[β], Zu)
            end
        end
    end
    return
end

# Sync the bare dipoles with the coherent state of an unit at `site`
function sync_bare_dipoles_at!(sys::System, site)
    Z = sys.coherents[site]
    ent = get_entanglement(sys)
    (; uncontracted, bare_dipole_operators) = ent
    for bs in bare_sites_for_unit(ent, site)
        uncontracted.dipoles[bs] = bare_expected_dipole(Z, bare_dipole_operators[to_atom(bs)])
    end
    return
end


################################################################################
# Measurements
################################################################################

# Transform an atom-level MeasureSpec (indexed to `uncontracted`) into a
# unit-level MeasureSpec (indexed to the contracted `sys`). For each unit and
# each of its `atoms_per_unit` subsites, the atom operator is embedded into the
# product-space Hilbert space via `local_op_to_product_space`. Position offsets
# and form factors come from `unit_map`.
function entangled_measure(measure, sys::System)
    @assert is_entangled(sys)               # System is entangled
    @assert size(measure.operators, 1) == 1 # Measure is for uncontracted system

    (; uncontracted, unit_map) = get_entanglement(sys)

    nobs = num_observables(measure)
    dims = sys.dims
    nunits = length(unit_map.unit_to_members)
    atoms_per_unit = length(unit_map.unit_to_members[1])  # uniform by construction

    Op = eltype(measure.operators)
    new_ops     = Array{Op, 6}(undef, atoms_per_unit, nobs, dims..., nunits)
    new_offsets = zeros(Vec3, atoms_per_unit, nunits)
    new_ff      = Array{FormFactor, 3}(undef, atoms_per_unit, nobs, nunits)

    for u in 1:nunits
        Ns_unit = uncontracted.Ns[atoms_in_unit(unit_map, u)]
        for (k, member) in enumerate(unit_map.unit_to_members[u])
            (; atom, Δpos, Δcell) = member
            new_offsets[k, u] = Δpos
            for μ in 1:nobs
                new_ff[k, μ, u] = measure.formfactors[1, μ, atom]
                for c in CartesianIndices(dims)
                    c′ = altmod1.(Tuple(c) .+ Δcell, dims)
                    A = measure.operators[1, μ, c′..., atom]
                    A_product = local_op_to_product_space(A, k, Ns_unit)
                    new_ops[k, μ, c, u] = Hermitian(A_product)
                end
            end
        end
    end

    return MeasureSpec(new_ops, measure.corr_pairs, measure.combiner, new_ff; offsets=new_offsets)
end
