################################################################################
# Types
################################################################################

# Base definition of the entangled units. Each unit is a list of atoms (for the
# original chemical cell) and their corresponding cell offsets.
const Groupings = Vector{Tuple{Int, SVector{3,Int}}}

# A member atom within an entangled unit.
struct MemberAtom
    atom        :: Int            # Atom index in the uncontracted crystal
    offset      :: Vec3           # Position offset from unit center (fractional coords)
    cell_offset :: SVector{3,Int} # Cell of this atom relative to the unit's cell (for straddling)
end

# Index to a part of an entangled
struct UnitPart
    unit :: Int                   # Unit index in contracted crystal
    part :: Int                   # Part index within unit 
end

# Bidirectional mapping between atoms (in uncontracted crystal) and entangled
# unit indices (in contracted crystal).
struct UnitMap
    atom_to_unit :: Vector{UnitPart}
    unit_to_members :: Vector{Vector{MemberAtom}}
end

# Metadata for a System with entangled units.  Hilbert dimension N of the
# uncontracted system is intentionally omitted to avoid dynamic lookup.
struct Entanglement <: AbstractEntanglement
    uncontracted          :: System                          # System prior to entanglement
    groupings             :: Vector{Groupings}               # Unit specification for chemical cell
    unit_map              :: UnitMap                         # Bidirectional atom ↔ unit mapping
    bare_dipole_operators :: Vector{NTuple{3, HermitianC64}} # Product-space spin ops per physical atom
end

# True if the system carries entanglement metadata.
is_entangled(sys::System) = !isnothing(sys.entanglement)

# Needed because uncontracted.dipoles and related state are mutably updated
function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.uncontracted), ent.groupings,
                        ent.unit_map, ent.bare_dipole_operators)
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
    # Identify which atoms are grouped into units vs. standalone
    entangled_atoms = Set{Int}(atom for unit in groupings for (atom, _) in unit)
    unentangled_atoms = filter(!in(entangled_atoms), 1:natoms(crystal))

    @assert length(entangled_atoms) + length(unentangled_atoms) == natoms(crystal) "Invalid entangled unit specification."

    # Build mapping structures directly (no intermediate dicts)
    atom_to_unit = Vector{UnitPart}(undef, natoms(crystal))
    unit_to_members = Vector{Vector{MemberAtom}}()
    new_positions = Vec3[]

    # Process unentangled atoms first (each becomes a singleton unit)
    for atom in unentangled_atoms
        unit_idx = length(unit_to_members) + 1
        push!(unit_to_members, [MemberAtom(atom, zero(Vec3), zero(SVector{3,Int}))])
        push!(new_positions, crystal.positions[atom])
        atom_to_unit[atom] = UnitPart(unit_idx, 1)
    end

    # Process multi-atom entangled units
    for unit_spec in groupings
        # Compute unwrapped positions (applying user-specified cell offsets)
        unwrapped_positions = [
            crystal.positions[atom] + Vec3(cell_offset...)
            for (atom, cell_offset) in unit_spec
        ]
        center_unwrapped = sum(unwrapped_positions) / length(unwrapped_positions)
        center_cell = floor.(Int, center_unwrapped)

        # Store wrapped unit center in the contracted crystal
        push!(new_positions, center_unwrapped - Vec3(center_cell...))

        # Build members list with offsets relative to the unit center
        unit_idx = length(unit_to_members) + 1
        members = map(enumerate(zip(unit_spec, unwrapped_positions))) do (k, ((atom, cell_offset), pos))
            offset = pos - center_unwrapped
            relative_cell_offset = cell_offset - center_cell
            atom_to_unit[atom] = UnitPart(unit_idx, k)
            MemberAtom(atom, offset, relative_cell_offset)
        end
        push!(unit_to_members, members)
    end

    # Generate contracted crystal and mapping
    # Space group = 1 allows "invalid" anisotropies from pair → onsite conversions
    new_crystal = Crystal(crystal.latvecs, new_positions, 1)
    unit_map = UnitMap(atom_to_unit, unit_to_members)

    return new_crystal, unit_map
end


# Returns a list of length equal to the number of "units" in the a contracted
# crystal. Each list element is itself a list of integers, each of which
# corresponds to the N of the corresponding site of the original system. The
# order is consistent with that given by the `inverse` field of a
# `UnitMap`.
function Ns_in_units(sys_original, unit_map)
    return [[sys_original.Ns[member.atom] for member in members]
            for members in unit_map.unit_to_members]
end

# Pull out atom indices comprising a unit
atoms_in_unit(unit_map, i) = [member.atom for member in unit_map.unit_to_members[i]]

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
# Cell of a chemical atom relative to the cell of the contracted unit that
# contains it.
@inline function unit_member_cell_offset(um::UnitMap, atom::Int)
    (; unit, part) = um.atom_to_unit[atom]
    return um.unit_to_members[unit][part].cell_offset
end

# Convert a bond between physical atoms into the corresponding bond between
# contracted units. If `bond.n` points from atom i in one cell to atom j in a
# neighboring cell, then the unit-to-unit cell displacement must be corrected by
# the two atoms' cell offsets inside their units.
function contracted_bond(bond::Bond, um::UnitMap)
    ui = um.atom_to_unit[bond.i]
    uj = um.atom_to_unit[bond.j]
    n = SVector{3,Int}(bond.n)
    n_new = n + unit_member_cell_offset(um, bond.i) - unit_member_cell_offset(um, bond.j)
    return Bond(ui.unit, uj.unit, collect(n_new))
end

# Checks whether a bond given in terms of the original crystal lies inside a
# single unit of the contracted system.
function bond_is_in_unit(bond::Bond, um::UnitMap)
    newbond = contracted_bond(bond, um)
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
function accum_pair_coupling_into_bond_operator_in_unit!(op, pc, sys, contracted_site, unit_map)
    pc.isculled && return

    Ns_unit = Ns_in_units(sys, unit_map)[contracted_site]
    i, j = pc.bond.i, pc.bond.j
    ui = unit_map.atom_to_unit[i]
    uj = unit_map.atom_to_unit[j]
    Ni = sys.Ns[1, 1, 1, i]
    Nj = sys.Ns[1, 1, 1, j]

    # Both sites lie in the same unit, so both embed into the same product space.
    embed_i = A -> local_op_to_product_space(A, ui.part, Ns_unit)
    embed_j = B -> local_op_to_product_space(B, uj.part, Ns_unit)
    accum_bond_operator!(op, pc, embed_i, embed_j, Ni, Nj)
    return
end

# Converts a pair coupling between two distinct units in the original system into
# a pair coupling between the corresponding units of the contracted system.
function pair_coupling_into_bond_operator_between_units(pc, sys, unit_map)
    (; i, j) = pc.bond
    ui = unit_map.atom_to_unit[i]
    uj = unit_map.atom_to_unit[j]

    Ns_local = Ns_in_units(sys, unit_map)
    Ns1 = Ns_local[ui.unit]
    Ns2 = Ns_local[uj.unit]
    N1 = sys.Ns[1, 1, 1, i]
    N2 = sys.Ns[1, 1, 1, j]
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
# `bare` system supplies only Hilbert-space dimensions (value-independent).
function contract_param(param::ModelParam, bare::System, unit_map, Ns_unit)
    nunits = length(unit_map.unit_to_members)

    # Intra-unit terms become onsite (on-bond) operators, one per touched unit.
    onsites = Tuple{Int, OnsiteCoupling}[]
    for u in 1:nunits
        Ns = Ns_unit[u]
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
                accum_pair_coupling_into_bond_operator_in_unit!(unit_operator, pc, bare, u, unit_map)
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
        (; newbond, bond_operator) = pair_coupling_into_bond_operator_between_units(pc, bare, unit_map)
        Ni = prod(Ns_unit[newbond.i])
        Nj = prod(Ns_unit[newbond.j])
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

The input `sys` should be in `:SUN` mode and have 1×1×1 dimensions. All
interactions in `sys` will be transferred to the entangled system. Subsequent
reshapings of this entangled system are then allowed.
"""
function entangle_system(sys::System{M}, groupings) where M
    isnothing(sys.origin) || error("Entangle a single-cell system first, then reshape")

    if eltype(eltype(groupings)) <: Integer
        groupings = [[(x, SA[0, 0, 0]) for x in g] for g in groupings]
    end

    groupings = [[(x, SVector{3, Int}(y)) for (x, y) in g] for g in groupings]

    # Construct contracted (P1) crystal and the atom ↔ unit mapping
    contracted_crystal, unit_map = contract_crystal(sys.crystal, groupings)

    # The (uniform) external field couples to each unit's total moment
    @assert allequal(@view sys.extfield[:,:,:,:]) "Entangled units require a uniform applied field."
    B = sys.extfield[1,1,1,1]

    # Local Hilbert space dimensions per unit (all must be equal). (TODO:
    # Determine if alternative behavior preferable in mixed case.)
    Ns_unit = Ns_in_units(sys, unit_map)
    Ns_contracted = map(prod, Ns_unit)
    @assert allequal(Ns_contracted) "After contraction, the dimensions of the local Hilbert spaces on each generalized site must all be equal."
    N = first(Ns_contracted)

    # Build the contracted `System{N}` fields directly (cf. `reshape_supercell_aux`).
    # An entangled system is an ordinary single-cell system that carries
    # entanglement metadata: unit g-factors are the identity, κ = 1, and dipoles
    # are undefined (NaN) because a unit's moment is a derived quantity. The
    # couplings are installed below through the ordinary params backbone.
    nunits = natoms(contracted_crystal)
    dims = sys.dims
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
                           copy(sys.rng), nothing)

    # Contract each labeled parameter independently (contraction is linear in
    # coupling strength), then let the ordinary machinery fill the interactions.
    sys_entangled.params = [contract_param(p, sys, unit_map, Ns_unit) for p in sys.params]
    repopulate_couplings_from_params!(sys_entangled)

    # Clone the physical (bare) system and derive all entanglement metadata from
    # the immutable `groupings` truth. Both `sys_entangled` and the physical
    # `uncontracted` start as single chemical cells (`origin = nothing`) and
    # reshape in tandem through the ordinary backbone.
    uncontracted = clone_system(sys)
    rebuild_entanglement!(sys_entangled, uncontracted, groupings)

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
# geometry and the immutable `groupings` truth (a grouping of chemical-cell
# atoms of the *unreshaped* physical crystal). The derived `unit_map` and the
# cached operator tables are rebuilt here geometrically. This is the single
# entry point used by both `entangle_system` and every reshaping variant.
#
# The mapping is derived from positions in the shared lattice frame (contracted
# and physical crystals share `latvecs`); each member's `cell_offset` records
# the cell of its atom relative to the unit's cell, so a unit's atoms may lie in
# different cells of the reshaped system (a "straddling" unit).
function rebuild_entanglement!(sys::System, uncontracted, groupings::Vector{Groupings})
    bare_origin = something(uncontracted.origin, uncontracted)

    # Chemical-cell contraction: per chemical atom, its within-unit `part` and
    # global displacement from the unit center. Reshape-invariant.
    _, um_chem = contract_crystal(bare_origin.crystal, groupings)
    latvecs = bare_origin.crystal.latvecs
    goff = map(um_chem.atom_to_unit) do (; unit, part)
        latvecs * um_chem.unit_to_members[unit][part].offset
    end
    parts = [part for (; part) in um_chem.atom_to_unit]

    # Build the reshaped, homogeneous `unit_map` sized to the reshaped
    # physical crystal. For each atom `a` of the reshaped chemical cell, locate
    # the contracted site sitting at its unit center (which may fall in a
    # neighboring cell ⇒ straddling). `inverse[u]` is ordered by part to align
    # with the product-space factor order.
    nunits = natoms(sys.crystal)
    atoms_per_unit = length(first(groupings))
    natoms_bare = natoms(uncontracted.crystal)
    @assert natoms_bare == nunits * atoms_per_unit
    atom_to_unit = UnitPart[]
    unit_to_members = [Vector{MemberAtom}(undef, atoms_per_unit) for _ in 1:nunits]
    for a in 1:natoms_bare
        a_chem = map_atom_to_other_crystal(uncontracted.crystal, a, bare_origin.crystal)
        part = parts[a_chem]
        center = global_position_at(uncontracted, CartesianIndex(1, 1, 1, a)) - goff[a_chem]
        unit_site = position_to_site(sys, orig_crystal(sys).latvecs \ center)
        u = to_atom(unit_site)
        push!(atom_to_unit, UnitPart(u, part))
        # Physical offset from the unit center (reshaped fractional frame; used
        # by `entangled_measure`) and the integer cell of atom `a` relative to
        # the unit's cell (nonzero ⇒ straddling; used by `member_site`).
        offset = uncontracted.crystal.latvecs \ goff[a_chem]
        cell_offset = SVector(1, 1, 1) .- SVector(to_cell(unit_site))
        unit_to_members[u][part] = MemberAtom(a, offset, cell_offset)
    end
    unit_map = UnitMap(atom_to_unit, unit_to_members)

    Ns_unit = Ns_in_units(uncontracted, unit_map)
    bare_dipole_operators =  map(1:natoms_bare) do atom
        S_local = spin_matrices_of_dim(; N=uncontracted.Ns[1, 1, 1, atom])
        (; unit, part) = unit_map.atom_to_unit[atom]
        ntuple(α -> local_op_to_product_space(S_local[α], part, Ns_unit[unit]), 3)
    end

    sys.entanglement = Entanglement(uncontracted, groupings, unit_map, bare_dipole_operators)
    return sys
end


################################################################################
# Dynamics: syncing physical dipoles with coherent states
################################################################################

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
            S = bare_dipole_operators[member.atom]
            S_bare[bs] = Vec3(real(dot(Zu, S[1], Zu)), real(dot(Zu, S[2], Zu)), real(dot(Zu, S[3], Zu)))
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
        S = bare_dipole_operators[member.atom]
        uncontracted.dipoles[bs] = Vec3(real(dot(Z, S[1], Z)), real(dot(Z, S[2], Z)), real(dot(Z, S[3], Z)))
    end
    return
end


################################################################################
# Coherent states and measurements
################################################################################
# Find the unique coherent state corresponding to a set of fully-polarized
# dipoles on each site inside a specified entangled unit. FIXME: This function
# is dead. Fix `set_dipole!` instead.
function coherent_state_from_dipoles(sys::System, dipoles, unit)
    (; uncontracted, unit_map) = get_entanglement(sys)

    # Atom indices (of the physical system) that lie in the specified unit.
    atoms = [member.atom for member in unit_map.unit_to_members[unit]]
    @assert length(dipoles) == length(atoms) "Invalid number of dipoles for specified unit."

    # Local Hilbert space dimensions for those atoms.
    Ns = Ns_in_units(uncontracted, unit_map)[unit]

    # Coherent state per atom, in each local Hilbert space.
    coherents = []
    for (dipole, N) in zip(dipoles, Ns)
        S = spin_matrices((N-1)/2)
        coherent = eigvecs(S' * dipole)[:,N] # Highest-weight eigenvector
        push!(coherents, coherent)
    end

    # Tensor product gives the unit's coherent state.
    return kron(coherents...)
end

# Build a unit-level MeasureSpec for an entangled `System`. The measure is first
# constructed at the atom level on the physical `uncontracted` (reusing the
# ordinary `ssf_custom` for g-tensor, form factors, correlation pairs, and
# combiner), then transformed into a unit-level measure via `entangled_measure`.
# `ssf_custom(f, sys::System)` dispatches here when `sys` is entangled.
function ssf_custom_entangled(f, sys::System; apply_g, formfactors)
    (; uncontracted) = get_entanglement(sys)
    measure_atom = ssf_custom(f, uncontracted; apply_g, formfactors)
    return entangled_measure(measure_atom, sys)
end

# Transform an atom-level MeasureSpec (indexed to `uncontracted`) into a
# unit-level MeasureSpec indexed to the contracted `sys`. For each unit and each
# of its `atoms_per_unit` subsites, the atom operator is embedded into the
# product-space Hilbert space via `local_op_to_product_space`. Position offsets
# and form factors come from `unit_map`. Observables are uniform across
# cells (g-factors are uniform), so the per-cell operator is broadcast.
function entangled_measure(measure, sys::System)
    (; uncontracted, unit_map) = get_entanglement(sys)
    Ns_unit = Ns_in_units(uncontracted, unit_map)

    nobs = num_observables(measure)
    dims = sys.dims

    nunits = length(unit_map.unit_to_members)
    atoms_per_unit = length(unit_map.unit_to_members[1])  # uniform by construction

    Op = eltype(measure.operators)
    new_ops     = Array{Op, 6}(undef, atoms_per_unit, nobs, dims..., nunits)
    new_offsets = zeros(Vec3, atoms_per_unit, nunits)
    new_ff      = Array{FormFactor, 3}(undef, atoms_per_unit, nobs, nunits)

    for i in 1:nunits
        for (k, member) in enumerate(unit_map.unit_to_members[i])
            atom = member.atom  # atom index within a chemical cell of uncontracted
            new_offsets[k, i] = member.offset
            for μ in 1:nobs
                new_ff[k, μ, i] = measure.formfactors[1, μ, atom]
                for c in CartesianIndices(dims)
                    A = measure.operators[1, μ, c, atom]
                    A_product = local_op_to_product_space(A, k, Ns_unit[i])
                    new_ops[k, μ, c, i] = Hermitian(A_product)
                end
            end
        end
    end

    return MeasureSpec(new_ops, measure.corr_pairs, measure.combiner, new_ff; offsets=new_offsets)
end
