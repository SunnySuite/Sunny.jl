################################################################################
# Types
################################################################################

# Base definition of the entangled units. Each unit is a list of atoms (for the
# original chemical cell) and their corresponding cell offsets.
const Groupings = Vector{Tuple{Int, SVector{3,Int}}}

# Index to a part of an entangled
struct UnitPart
    unit :: Int                   # Unit index in contracted crystal
    part :: Int                   # Part index within unit
end

# A member atom within an entangled unit. Both fields are expressed in the basis
# of the crystal of the containing `UnitMap`, and record this atom's position
# relative to the unit center: `relpos` is the fractional position offset (in
# the same fractional coordinates as `crystal.positions`), and `cell_offset` is
# the integer cell (nonzero ⇒ the unit straddles a cell boundary). Because they
# are basis-relative, the same physical member has different values in a
# `UnitMap` keyed to the chemical origin crystal versus one keyed to a reshaped
# crystal.
struct MemberAtom
    atom        :: Int            # Atom index in the associated crystal
    relpos      :: Vec3           # Fractional position relative to unit center
    cell_offset :: SVector{3,Int} # Cell of this atom relative to the unit's cell (for straddling)
end

# Bidirectional mapping between atoms and entangled unit indices, keyed to a
# particular crystal (see the caveat on `MemberAtom`). An `Entanglement` stores
# one `unit_map` keyed to the reshaped uncontracted crystal, rebuilt
# geometrically per reshape by `rebuild_entanglement!`. Reshape-invariant uses
# (e.g. `contract_param`, which folds couplings in the chemical/origin frame)
# take the `unit_map` of the invariant origin system, `something(sys.origin, sys)`.
struct UnitMap
    atom_to_unit    :: Vector{UnitPart}           # atom → (unit, part)
    unit_to_members :: Vector{Vector{MemberAtom}} # unit → its members, ordered by part
end

# Metadata for a System with entangled units. Hilbert dimension N of the
# uncontracted system is intentionally omitted to avoid dynamic lookup.
struct Entanglement <: AbstractEntanglement
    uncontracted          :: System                          # System prior to entanglement
    unit_map              :: UnitMap                         # Map on the reshaped uncontracted crystal
    bare_dipole_operators :: Vector{NTuple{3, HermitianC64}} # Product-space spin ops per physical atom
end

# True if the system carries entanglement metadata.
is_entangled(sys::System) = !isnothing(sys.entanglement)

# Needed because uncontracted.dipoles and related state are mutably updated
function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.uncontracted), ent.unit_map,
                        ent.bare_dipole_operators)
end

clone_entanglement(::Nothing) = nothing



# Concrete-typed view of the entanglement metadata. The explicit type ascription
# avoids JET warnings.
@inline get_entanglement(sys::System) = sys.entanglement :: Entanglement


################################################################################
# Crystal contraction
################################################################################
# Takes a crystal and a list of entangled units. Each member of each unit is
# `(atom, cell_offset)`, where `atom` is a chemical-cell atom index and
# `cell_offset` chooses the periodic image of that atom in the unwrapped unit
# definition. Nonzero offsets allow a unit to straddle a periodic boundary.
function contract_crystal(crystal, groupings::Vector{Groupings})
    @assert sum(length, groupings; init=0) == natoms(crystal) "Invalid entangled unit specification."

    atom_to_unit = Vector{UnitPart}(undef, natoms(crystal))
    unit_to_members = Vector{MemberAtom}[]
    new_positions = Vec3[]

    for unit_spec in groupings
        # Compute unwrapped positions (applying user-specified cell offsets)
        unwrapped_positions = [
            crystal.positions[atom] + Vec3(cell_offset...)
            for (atom, cell_offset) in unit_spec
        ]
        center_unwrapped = sum(unwrapped_positions) / length(unwrapped_positions)
        # Integer cell of the center, computed tolerantly (-ϵ maps to 0)
        center_cell = round.(Int, center_unwrapped - wrap_to_unit_cell(center_unwrapped))

        # Store wrapped unit center in the contracted crystal
        push!(new_positions, center_unwrapped - Vec3(center_cell...))

        # Record each atom's part, its cell relative to the unit center, and its
        # cartesian displacement from the center.
        unit_idx = length(unit_to_members) + 1
        members = map(enumerate(zip(unit_spec, unwrapped_positions))) do (part, ((atom, cell_offset), pos))
            atom_to_unit[atom] = UnitPart(unit_idx, part)
            MemberAtom(atom, pos - center_unwrapped, cell_offset - center_cell)
        end
        push!(unit_to_members, members)
    end

    # Space group = 1 allows "invalid" anisotropies from pair → onsite conversions
    new_crystal = Crystal(crystal.latvecs, new_positions, 1)
    unit_map = UnitMap(atom_to_unit, unit_to_members)

    return new_crystal, unit_map
end


# Hilbert-space dimensions N for a list of physical atoms, in the given order.
Ns_for_atoms(uncontracted, atoms) = [uncontracted.Ns[a] for a in atoms]

# Hilbert-space dimensions N for each part of unit `u`, ordered to match the
# `unit_to_members` field of a `UnitMap`.
Ns_at_unit(uncontracted, unit_map, u) = Ns_for_atoms(uncontracted, atoms_in_unit(unit_map, u))

# Pull out atom indices comprising a unit
atoms_in_unit(unit_map::UnitMap, i) = [member.atom for member in unit_map.unit_to_members[i]]

# The `MemberAtom` record for a given atom (inverse of `atom_to_unit`).
function member_for_atom(unit_map::UnitMap, atom)
    (; unit, part) = unit_map.atom_to_unit[atom]
    return unit_map.unit_to_members[unit][part]
end

# Physical (bare-system) site of a unit member, given the containing unit's site
# and the member's `MemberAtom`. The member atom sits in cell `cell_offset`
# relative to the unit's cell (periodically wrapped), so a unit may straddle
# cell boundaries. Mirrors `bonded_site`.
@inline function member_site(unit_site, member::MemberAtom, dims)
    CartesianIndex(altmod1.(to_cell(unit_site) .+ Tuple(member.cell_offset), dims)..., member.atom)
end

# Physical (bare-system) sites comprising the entangled unit at contracted
# `unit_site` of `sys`. Accepts CartesianIndex, tuple, or integer atom index.
function entangled_unit_members(sys::System, unit_site::Union{CartesianIndex,NTuple{4,Int}})
    (; uncontracted, unit_map) = get_entanglement(sys)
    unit_site = to_cartesian(unit_site)
    return [member_site(unit_site, member, uncontracted.dims) for member in unit_map.unit_to_members[to_atom(unit_site)]]
end

function entangled_unit_members(sys::System, atom::Int)
    return entangled_unit_members(sys, CartesianIndex(1, 1, 1, atom))
end


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
    n_new = n + member_for_atom(unit_map, bond.i).cell_offset - member_for_atom(unit_map, bond.j).cell_offset
    return Bond(ui.unit, uj.unit, collect(n_new))
end

# Checks whether a bond given in terms of the original crystal lies inside a
# single unit of the contracted system.
function bond_is_in_unit(bond::Bond, unit_map::UnitMap)
    newbond = contracted_bond(bond, unit_map)
    return newbond.i == newbond.j && all(iszero, newbond.n)
end

# Accumulate the operator form of a `PairCoupling` into `op`. The two sites'
# local operators are lifted into `op`'s Hilbert space by the caller-supplied
# closures `embed_i`/`embed_j`; `Ni`/`Nj` are the local Hilbert space dimensions
# of the two original sites. This is the shared kernel of the intra-unit and
# inter-unit conversions below, which differ only in that embedding.
function accum_bond_operator!(op, pc::PairCoupling, embed_i, embed_j, Ni, Nj)
    (; scalar, bilin, biquad, general) = pc
    Ntot = size(op, 1)

    # Scalar part
    op .+= scalar * I(Ntot)

    # Bilinear part
    J = bilin isa Float64 ? bilin*I(3) : bilin
    Si = [embed_i(Sa) for Sa in spin_matrices((Ni-1)/2)]
    Sj = [embed_j(Sb) for Sb in spin_matrices((Nj-1)/2)]
    op .+= Si' * J * Sj

    # Biquadratic part
    K = biquad isa Float64 ? diagm(biquad * scalar_biquad_metric) : biquad
    Oi = [embed_i(Oa) for Oa in stevens_matrices_of_dim(2; N=Ni)]
    Oj = [embed_j(Ob) for Ob in stevens_matrices_of_dim(2; N=Nj)]
    op .+= Oi' * K * Oj

    # General part
    for (A, B) in general.data
        op .+= embed_i(A) * embed_j(B)
    end

    return op
end

# Converts what was a pair coupling between two different sites of a single unit
# in the original system into an on-bond operator (an onsite operator in terms
# of the "units"). Accumulates into `op` in place.
function accum_pair_coupling_into_bond_operator_in_unit!(op, pc, uncontracted, unit_map::UnitMap)
    pc.isculled && return

    i, j = pc.bond.i, pc.bond.j
    ui = unit_map.atom_to_unit[i]
    uj = unit_map.atom_to_unit[j]
    Ns_unit = Ns_at_unit(uncontracted, unit_map, ui.unit)
    Ni = uncontracted.Ns[1, 1, 1, i]
    Nj = uncontracted.Ns[1, 1, 1, j]

    # Both sites lie in the same unit, so both embed into the same product space.
    embed_i = A -> local_op_to_product_space(A, ui.part, Ns_unit)
    embed_j = B -> local_op_to_product_space(B, uj.part, Ns_unit)
    accum_bond_operator!(op, pc, embed_i, embed_j, Ni, Nj)
    return
end

# Converts a pair coupling between two distinct units in the original system into
# a pair coupling between the corresponding units of the contracted system.
function pair_coupling_into_bond_operator_between_units(pc, uncontracted, unit_map::UnitMap)
    (; i, j) = pc.bond
    ui = unit_map.atom_to_unit[i]
    uj = unit_map.atom_to_unit[j]

    Ns1 = Ns_at_unit(uncontracted, unit_map, ui.unit)
    Ns2 = Ns_at_unit(uncontracted, unit_map, uj.unit)
    N1 = uncontracted.Ns[1, 1, 1, i]
    N2 = uncontracted.Ns[1, 1, 1, j]
    dim1 = prod(Ns1)
    dim2 = prod(Ns2)
    N = dim1 * dim2

    # Site i acts on unit1 (first tensor factor), site j on unit2 (second).
    embed_i = A -> kron(local_op_to_product_space(A, ui.part, Ns1), I(dim2))
    embed_j = B -> kron(I(dim1), local_op_to_product_space(B, uj.part, Ns2))

    bond_operator = zeros(ComplexF64, N, N)
    accum_bond_operator!(bond_operator, pc, embed_i, embed_j, N1, N2)
    newbond = contracted_bond(pc.bond, unit_map)
    return (; newbond, bond_operator)
end

# Contract a single bare-system `ModelParam` into the corresponding `ModelParam`
# for the contracted system. Contraction is linear in the coupling strength, so
# each labeled parameter is contracted independently at unit strength; the
# ordinary `repopulate_couplings_from_params!` then rescales by `param.val`. The
# `uncontracted` system supplies only Hilbert-space dimensions (value-independent).
function contract_param(param::ModelParam, uncontracted::System, unit_map::UnitMap, Ns_units)
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
                accum_pair_coupling_into_bond_operator_in_unit!(unit_operator, pc, uncontracted, unit_map)
            end
        end

        iszero(unit_operator) || push!(onsites, (u, Hermitian(unit_operator)))
    end

    # Inter-unit pair couplings become pair couplings between contracted units.
    # Both bond directions are present in `param.pairs` (from symmetry
    # propagation on the bare system), so no re-propagation is needed here.
    pairs = PairCoupling[]
    for pc in param.pairs
        bond_is_in_unit(pc.bond, unit_map) && continue
        (; newbond, bond_operator) = pair_coupling_into_bond_operator_between_units(pc, uncontracted, unit_map)
        Ni = prod(Ns_units[newbond.i])
        Nj = prod(Ns_units[newbond.j])
        # `extract_parts=false` because dipoles are NaN for entangled systems.
        scalar, bilin, biquad, tensordec = decompose_general_coupling(bond_operator, Ni, Nj; extract_parts=false)
        push!(pairs, PairCoupling(newbond, scalar, bilin, biquad, tensordec))
    end

    return ModelParam(param.label, param.val, onsites, pairs)
end


################################################################################
# Public API
################################################################################
"""
    entangle_system(sys::System{N}, groupings)

Create a new [`System`](@ref), where `groupings` specifies how atoms are
collected into "entangled units". Each unit member is `(atom, cell_offset)`,
where `atom` is a chemical-cell atom index and `cell_offset::SVector{3,Int}`
chooses the periodic image of that atom in the unwrapped unit definition. Sunny
will model each such unit within a tensor-product Hilbert space to allow for
local quantum entanglement. This feature currently requires a single unit type
(all dimers, all trimers, etc).

For example, `[(1, [0, 0, 0]), (2, [1, 0, 0])]` defines a dimer whose second
atom is taken from the neighboring +a1 cell.

The input system should be in `:SUN` mode and have 1×1×1 dimensions. All its
interactions will be transferred to the entangled system. Subsequent reshapings
of this entangled system are then allowed.
"""
function entangle_system(uncontracted::System{M}, groupings) where M
    isnothing(uncontracted.origin) || error("Entangle a single-cell system first, then reshape")

    if eltype(eltype(groupings)) <: Integer
        groupings = [[(x, SA[0, 0, 0]) for x in g] for g in groupings]
    end

    groupings = [[(x, SVector{3, Int}(y)) for (x, y) in g] for g in groupings]

    # Construct contracted (P1) crystal and the chemical-cell map
    contracted_crystal, unit_map = contract_crystal(uncontracted.crystal, groupings)

    # The (uniform) external field couples to each unit's total moment
    @assert allequal(@view uncontracted.extfield[:,:,:,:]) "Entangled units require a uniform applied field."
    B = uncontracted.extfield[1,1,1,1]

    # Local Hilbert space dimensions per unit (all must be equal). (TODO:
    # Determine if alternative behavior preferable in mixed case.)
    Ns_units = [Ns_at_unit(uncontracted, unit_map, u) for u in eachindex(unit_map.unit_to_members)]
    Ns_contracted = map(prod, Ns_units)
    @assert allequal(Ns_contracted) "After contraction, the dimensions of the local Hilbert spaces on each generalized site must all be equal."
    N = first(Ns_contracted)

    # Build the contracted `System{N}` fields directly (cf. `reshape_supercell_aux`).
    # An entangled system is an ordinary single-cell system that carries
    # entanglement metadata: unit g-factors are the identity, κ = 1, and dipoles
    # are undefined (NaN) because a unit's moment is a derived quantity. The
    # couplings are installed below through the ordinary params backbone.
    nunits    = natoms(contracted_crystal)
    dims      = uncontracted.dims
    Ns        = reshape(collect(Ns_contracted), 1, 1, 1, :)
    κs        = ones(Float64, dims..., nunits)
    gs        = fill(Mat3(I), dims..., nunits)
    ints      = empty_interactions(:SUN, nunits, N)
    extfield  = fill(B, dims..., nunits)
    dipoles   = fill(Vec3(NaN, NaN, NaN), dims..., nunits)
    coherents = zeros(CVec{N}, dims..., nunits)
    sys_entangled = System(nothing, :SUN, contracted_crystal, dims, Ns, κs, gs,
                           ModelParam[], Symbol[], ints, nothing, extfield,
                           dipoles, coherents, Array{Vec3, 4}[], Array{CVec{N}, 4}[],
                           copy(uncontracted.rng), nothing)

    # Contract each labeled parameter independently (contraction is linear in
    # coupling strength), then let the ordinary machinery fill the interactions.
    sys_entangled.params = [contract_param(p, uncontracted, unit_map, Ns_units) for p in uncontracted.params]
    repopulate_couplings_from_params!(sys_entangled)

    # Clone the physical system and derive all entanglement metadata from the
    # immutable `groupings` truth. Both `sys_entangled` and the physical
    # `uncontracted` start as single chemical cells (`origin = nothing`) and
    # reshape in tandem through the ordinary backbone.
    uncontracted = clone_system(uncontracted)
    rebuild_entanglement!(sys_entangled, uncontracted, unit_map)

    # Initialize coherent states of each entangled unit to the tensor product of
    # bare coherent states.
    for unit in eachsite(sys_entangled)
        members = entangled_unit_members(sys_entangled, unit)
        Z = kron((uncontracted.coherents[bs] for bs in members)...)
        set_coherent!(sys_entangled, Z, unit)
    end

    return sys_entangled
end

function entangle_system(::System{0}, _)
    error("Cannot entangle units of a `:dipole`-mode `System`. Use `:SUN` mode.")
end

# Populate the entanglement metadata of a contracted `sys` (which may be
# reshaped), given the physical `uncontracted` system reshaped to the same
# geometry and the reshape-invariant `orig_unit_map` (the `UnitMap` keyed to the
# chemical/origin crystal). This `orig_unit_map` depends only on the origin
# crystal and the `groupings` truth, both fixed for the system's lifetime: on
# the first call it comes from `contract_crystal`, and on every reshape it is
# the (invariant) origin's own `unit_map`, so it is never recomputed. Only the
# geometric `unit_map` and operator tables are rebuilt here. This is the single
# entry point used by both `entangle_system` and every reshaping variant.
#
# The mapping is derived from positions in the shared lattice frame (contracted
# and physical crystals share `latvecs`); each member's `cell_offset` records
# the cell of its atom relative to the unit's cell, so a unit's atoms may lie in
# different cells of the reshaped system (a "straddling" unit).
function rebuild_entanglement!(sys::System, uncontracted, orig_unit_map::UnitMap)
    # Chemical-cell atom indices are resolved against the origin crystal
    # (`orig_unit_map` is keyed by them); the origin is invariant across reshapes.
    bare_origin = something(uncontracted.origin, uncontracted)

    # Build the reshaped, homogeneous `unit_map` sized to the reshaped
    # physical crystal. For each atom `a` of the reshaped chemical cell, locate
    # the contracted site sitting at its unit center (which may fall in a
    # neighboring cell ⇒ straddling). `inverse[u]` is ordered by part to align
    # with the product-space factor order.
    nunits = natoms(sys.crystal)
    atoms_per_unit = length(first(orig_unit_map.unit_to_members))
    natoms_bare = natoms(uncontracted.crystal)
    @assert natoms_bare == nunits * atoms_per_unit
    atom_to_unit = UnitPart[]
    unit_to_members = [Vector{MemberAtom}(undef, atoms_per_unit) for _ in 1:nunits]
    for a in 1:natoms_bare
        a_chem = map_atom_to_other_crystal(uncontracted.crystal, a, bare_origin.crystal)
        chem_member = member_for_atom(orig_unit_map, a_chem)
        part = orig_unit_map.atom_to_unit[a_chem].part
        # The physical displacement from the unit center is reshape-invariant, so
        # convert the chem-frame relative position to Cartesian to locate the unit
        # center, then re-express it in the reshaped (uncontracted) frame.
        cart_disp = bare_origin.crystal.latvecs * chem_member.relpos
        center = global_position_at(uncontracted, CartesianIndex(1, 1, 1, a)) - cart_disp
        # Locate the unit at this center. Take the *unwrapped* cell (rather than
        # `position_to_site`, which wraps modulo `dims`) so that `cell_offset`
        # records the genuine straddle even on the 1×1×1 origin, where every
        # wrapped cell collapses to (1,1,1). Periodicity is reapplied at use in
        # `member_site`.
        u, unit_cell = position_to_atom_and_offset(sys.crystal, sys.crystal.latvecs \ center)
        push!(atom_to_unit, UnitPart(u, part))
        # Fractional position relative to the unit center (in the reshaped frame;
        # used by `entangled_measure`) and the integer cell of atom `a` relative
        # to the unit's cell (nonzero ⇒ straddling; used by `member_site`). Atom
        # `a` sits in cell (1,1,1)==0 (0-based), so its cell relative to the
        # unit's is `-unit_cell`.
        relpos = uncontracted.crystal.latvecs \ cart_disp
        cell_offset = -SVector(unit_cell)
        unit_to_members[u][part] = MemberAtom(a, relpos, cell_offset)
    end
    unit_map = UnitMap(atom_to_unit, unit_to_members)

    bare_dipole_operators =  map(1:natoms_bare) do atom
        S_local = spin_matrices_of_dim(; N=uncontracted.Ns[1, 1, 1, atom])
        (; unit, part) = unit_map.atom_to_unit[atom]
        ntuple(α -> local_op_to_product_space(S_local[α], part, Ns_at_unit(uncontracted, unit_map, unit)), 3)
    end

    sys.entanglement = Entanglement(uncontracted, unit_map, bare_dipole_operators)
    return sys
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
    (; uncontracted, unit_map, bare_dipole_operators) = ent
    dims = uncontracted.dims

    # Reuse the bare system's persistent scratch buffers.
    S_bare, dE_dS_bare = get_dipole_buffers(uncontracted, 2)
    fill!(dE_dS_bare, zero(Vec3))

    # Physical dipoles ⟨Sₐ⟩ evaluated from the *passed* `Z` (which may be a
    # trial state not yet synced into `uncontracted.dipoles`), so the
    # dipole-dependent Ewald field is consistent with the state whose gradient
    # is requested.
    for u in CartesianIndices(Z)
        Zu = Z[u]
        for member in unit_map.unit_to_members[to_atom(u)]
            bs = member_site(u, member, dims)
            S_bare[bs] = bare_expected_dipole(Zu, bare_dipole_operators[member.atom])
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
    # gradient. A unit's members may straddle cells; `member_site` routes them.
    for u in CartesianIndices(Z)
        Zu = Z[u]
        for member in unit_map.unit_to_members[to_atom(u)]
            bs = member_site(u, member, dims)
            S = bare_dipole_operators[member.atom]
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
    (; uncontracted, unit_map, bare_dipole_operators) = get_entanglement(sys)
    for member in unit_map.unit_to_members[to_atom(site)]
        bs = member_site(site, member, uncontracted.dims)
        uncontracted.dipoles[bs] = bare_expected_dipole(Z, bare_dipole_operators[member.atom])
    end
    return
end


################################################################################
# Measurements
################################################################################

# Transform an atom-level MeasureSpec (indexed to `uncontracted`) into a
# unit-level MeasureSpec indexed to the contracted `sys`. For each unit and each
# of its `atoms_per_unit` subsites, the atom operator is embedded into the
# product-space Hilbert space via `local_op_to_product_space`. Position offsets
# and form factors come from `unit_map`. Observables are uniform across
# cells (g-factors are uniform), so the per-cell operator is broadcast.
function entangled_measure(measure, sys::System)
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
        Ns_unit = Ns_at_unit(uncontracted, unit_map, u)
        for (k, member) in enumerate(unit_map.unit_to_members[u])
            atom = member.atom  # atom index within a chemical cell of uncontracted
            # `MeasureSpec.offsets` is added to reshaped-fractional positions `r`
            # in the SWT phase factor; `relpos` is already in that frame.
            new_offsets[k, u] = member.relpos
            for μ in 1:nobs
                new_ff[k, μ, u] = measure.formfactors[1, μ, atom]
                for c in CartesianIndices(dims)
                    A = measure.operators[1, μ, c, atom]
                    A_product = local_op_to_product_space(A, k, Ns_unit)
                    new_ops[k, μ, c, u] = Hermitian(A_product)
                end
            end
        end
    end

    return MeasureSpec(new_ops, measure.corr_pairs, measure.combiner, new_ff; offsets=new_offsets)
end
