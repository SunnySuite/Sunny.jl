################################################################################
# Types
################################################################################
# A member atom within an entangled unit.
struct UnitMember
    atom        :: Int            # Atom index in the physical crystal
    offset      :: Vec3           # Position offset from unit center (fractional coords)
    cell_offset :: SVector{3,Int} # Cell of this atom relative to the unit's cell (for straddling)
end

# User-facing specification of an entangled unit. The atom index identifies a
# chemical-cell atom, and the cell offset chooses which periodic image of that
# atom participates in the unwrapped unit definition.
const UnitSpec = Vector{Tuple{Int, SVector{3,Int}}}

@inline floor_cell(v) = SVector{3,Int}(floor(Int, v[1]), floor(Int, v[2]), floor(Int, v[3]))

# Bidirectional mapping between physical atoms and entangled units.
# - atom_to_unit[i] = (unit_index, slot_in_unit) for physical atom i
# - unit_members[j] = ordered list of UnitMembers comprising unit j
struct UnitMap
    atom_to_unit :: Vector{Tuple{Int, Int}}
    unit_members :: Vector{Vector{UnitMember}}
end

function Base.copy(um::UnitMap)
    return UnitMap(copy(um.atom_to_unit), copy(um.unit_members))
end

# Metadata for a System with entangled units. Hilbert dimension N is _not_
# carried in the bare system's type. This avoids unecessary dynamic lookups.
#
# `units` is the immutable "truth": the grouping of chemical-cell atoms into
# entangled units, expressed against the *unreshaped* physical crystal
# (`something(bare_system.origin, bare_system)`). Everything else is derived and
# rebuilt by `rebuild_entanglement!` whenever the system is reshaped —
# `unit_map` and the operator tables track the current (possibly reshaped)
# `bare_system`. The physical sites comprising a unit are recovered by integer
# arithmetic from `unit_map` (see `member_site`); each member carries a
# `cell_offset` so a unit may straddle cell boundaries of a reshaped system.
struct Entanglement <: AbstractEntanglement
    bare_system           :: System                                    # Physical system
    units                 :: Vector{UnitSpec}                          # Unit specification for chemical cell
    unit_map              :: UnitMap                                   # Bidirectional atom ↔ unit mapping
    bare_dipole_operators :: Vector{NTuple{3, HermitianC64}}           # Product-space spin ops per physical atom (no g)
end

function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.bare_system), deepcopy(ent.units),
                        copy(ent.unit_map), ent.bare_dipole_operators)
end

# Concrete-typed view of the entanglement metadata. The field is declared
# `Union{Nothing, AbstractEntanglement}` on `System` (the abstract supertype
# breaks a recursive type dependency), so accessing it directly poisons
# inference with an abstract type. Mirroring the `interactions_union ::
# Vector{Interactions}` pattern, assert the concrete `Entanglement` here so JET
# and downstream inference stay clean. The `Entanglement{N}` parameter is the
# bare system's local dimension (generally distinct from the contracted `sys`'s
# N), so it is left free here; asserting the UnionAll still concretizes the
# `bare_system::System{N}` field per instance. Only call when `sys` is entangled.
@inline get_entanglement(sys::System) = sys.entanglement :: Entanglement


################################################################################
# Crystal contraction
################################################################################
# Takes a crystal and a list of entangled units. Each member of each unit is
# `(atom, cell_offset)`, where `atom` is a chemical-cell atom index and
# `cell_offset` chooses the periodic image of that atom in the unwrapped unit
# definition. Nonzero offsets allow a unit to straddle a periodic boundary.
function contract_crystal(crystal, units::Vector{UnitSpec})
    # Identify which atoms are grouped into units vs. standalone
    entangled_atoms = Set{Int}(atom for unit in units for (atom, _) in unit)
    unentangled_atoms = filter(∉(entangled_atoms), 1:natoms(crystal))

    @assert length(entangled_atoms) + length(unentangled_atoms) == natoms(crystal) "Invalid entangled unit specification."

    # Build mapping structures directly (no intermediate dicts)
    atom_to_unit = Vector{Tuple{Int, Int}}(undef, natoms(crystal))
    unit_members = Vector{Vector{UnitMember}}()
    new_positions = Vec3[]

    # Process unentangled atoms first (each becomes a singleton unit)
    for atom in unentangled_atoms
        unit_idx = length(unit_members) + 1
        push!(unit_members, [UnitMember(atom, zero(Vec3), zero(SVector{3,Int}))])
        push!(new_positions, crystal.positions[atom])
        atom_to_unit[atom] = (unit_idx, 1)
    end

    # Process multi-atom entangled units
    for unit_spec in units
        # Compute unwrapped positions (applying user-specified cell offsets)
        unwrapped_positions = [
            crystal.positions[atom] + Vec3(cell_offset...)
            for (atom, cell_offset) in unit_spec
        ]
        center_unwrapped = sum(unwrapped_positions) / length(unwrapped_positions)
        center_cell = floor_cell(center_unwrapped)

        # Store wrapped unit center in the contracted crystal
        push!(new_positions, center_unwrapped - Vec3(center_cell...))

        # Build members list with offsets relative to the unit center
        unit_idx = length(unit_members) + 1
        members = map(enumerate(zip(unit_spec, unwrapped_positions))) do (k, ((atom, cell_offset), pos))
            offset = pos - center_unwrapped
            relative_cell_offset = cell_offset - center_cell
            atom_to_unit[atom] = (unit_idx, k)
            UnitMember(atom, offset, relative_cell_offset)
        end
        push!(unit_members, members)
    end

    # Generate contracted crystal and mapping
    # Space group = 1 allows "invalid" anisotropies from pair → onsite conversions
    new_crystal = Crystal(crystal.latvecs, new_positions, 1)
    unit_map = UnitMap(atom_to_unit, unit_members)

    return new_crystal, unit_map
end


# Returns a list of length equal to the number of "units" in the a contracted
# crystal. Each list element is itself a list of integers, each of which
# corresponds to the N of the corresponding site of the original system. The
# order is consistent with that given by the `inverse` field of a
# `UnitMap`.
function Ns_in_units(sys_original, unit_map)
    return [[sys_original.Ns[member.atom] for member in members]
            for members in unit_map.unit_members]
end

# Pull out atom indices comprising a unit
atoms_in_unit(unit_map, i) = [member.atom for member in unit_map.unit_members[i]]

# Physical (bare-system) site of a unit member, given the containing unit's site
# and the member's `UnitMember`. The member atom sits in cell `cell_offset`
# relative to the unit's cell (periodically wrapped), so a unit may straddle
# cell boundaries. Mirrors `bonded_site`.
@inline function member_site(unit_site, member::UnitMember, dims)
    CartesianIndex(altmod1.(to_cell(unit_site) .+ Tuple(member.cell_offset), dims)..., member.atom)
end

# Physical (bare-system) sites comprising the entangled unit at contracted
# `unit_site` of `sys`. Accepts CartesianIndex, tuple, or integer atom index.
function entangled_unit_members(sys::System, unit_site::Union{CartesianIndex,NTuple{4,Int}})
    (; bare_system, unit_map) = get_entanglement(sys)
    unit_site = to_cartesian(unit_site)
    return [member_site(unit_site, member, bare_system.dims) for member in unit_map.unit_members[to_atom(unit_site)]]
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
    unit, k = um.atom_to_unit[atom]
    return um.unit_members[unit][k].cell_offset
end

# Convert a bond between physical atoms into the corresponding bond between
# contracted units. If `bond.n` points from atom i in one cell to atom j in a
# neighboring cell, then the unit-to-unit cell displacement must be corrected by
# the two atoms' cell offsets inside their units.
function contracted_bond(bond::Bond, um::UnitMap)
    unit_i, _ = um.atom_to_unit[bond.i]
    unit_j, _ = um.atom_to_unit[bond.j]
    n = SVector{3,Int}(bond.n)
    n_new = n + unit_member_cell_offset(um, bond.i) - unit_member_cell_offset(um, bond.j)
    return Bond(unit_i, unit_j, collect(n_new))
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
    i_unit = unit_map.atom_to_unit[i][2]
    j_unit = unit_map.atom_to_unit[j][2]
    Ni = sys.Ns[1, 1, 1, i]
    Nj = sys.Ns[1, 1, 1, j]

    # Both sites lie in the same unit, so both embed into the same product space.
    embed_i = A -> local_op_to_product_space(A, i_unit, Ns_unit)
    embed_j = B -> local_op_to_product_space(B, j_unit, Ns_unit)
    accum_bond_operator!(op, pc, embed_i, embed_j, Ni, Nj)
    return
end

# Converts a pair coupling between two distinct units in the original system into
# a pair coupling between the corresponding units of the contracted system.
function pair_coupling_into_bond_operator_between_units(pc, sys, unit_map)
    (; i, j) = pc.bond
    unit1, slot1 = unit_map.atom_to_unit[i]
    unit2, slot2 = unit_map.atom_to_unit[j]

    Ns_local = Ns_in_units(sys, unit_map)
    Ns1 = Ns_local[unit1]
    Ns2 = Ns_local[unit2]
    N1 = sys.Ns[1, 1, 1, i]
    N2 = sys.Ns[1, 1, 1, j]
    dim1 = prod(Ns1)
    dim2 = prod(Ns2)
    N = dim1 * dim2

    # Site i acts on unit1 (first tensor factor), site j on unit2 (second).
    embed_i = A -> kron(local_op_to_product_space(A, slot1, Ns1), I(dim2))
    embed_j = B -> kron(I(dim1), local_op_to_product_space(B, slot2, Ns2))

    bond_operator = zeros(ComplexF64, N, N)
    accum_bond_operator!(bond_operator, pc, embed_i, embed_j, N1, N2)
    newbond = contracted_bond(pc.bond, unit_map)
    return (; newbond, bond_operator)
end


################################################################################
# Public API
################################################################################
"""
    entangle_system(sys::System{N}, units)

Create a new [`System`](@ref), where `units` are collected into "entangled
units". Each unit member is `(atom, cell_offset)`, where `atom` is a
chemical-cell atom index and `cell_offset::SVector{3,Int}` chooses the periodic
image of that atom in the unwrapped unit definition. Sunny will model each such
unit within a tensor-product Hilbert space to allow for local quantum
entanglement. This feature currently requires a single unit type (all dimers,
all trimers, etc).

For example, `[(1, [0, 0, 0]), (2, [1, 0, 0])]` defines a dimer whose
second atom is taken from the neighboring +a1 cell.

The input `sys` should be in `:SUN` mode and have 1×1×1 dimensions. All
interactions in `sys` will be transferred to the entangled system. Subsequent
reshapings of this entangled system are then allowed.
"""
function entangle_system(sys::System{M}, units) where M
    isnothing(sys.origin) || error("Entangle a single-cell system first, then reshape")

    if eltype(eltype(units)) <: Integer
        units = [[(x, SA[0, 0, 0]) for x in g] for g in units]
    end

    units = [[(x, SVector{3, Int}(y)) for (x, y) in g] for g in units]

    # Construct contracted crystal
    contracted_crystal, unit_map = contract_crystal(sys.crystal, units)

    # Make sure we have a uniform external field
    @assert allequal(@view sys.extfield[:,:,:,:]) "Entangled units require a uniform applied field."
    B = sys.extfield[1,1,1,1]

    # Determine Ns for local Hilbert spaces (all must be equal). (TODO:
    # Determine if alternative behavior preferable in mixed case.)
    Ns_unit = Ns_in_units(sys, unit_map)
    Ns_contracted = map(Ns -> prod(Ns), Ns_unit)
    @assert allequal(Ns_contracted) "After contraction, the dimensions of the local Hilbert spaces on each generalized site must all be equal."

    # Construct empty contracted system
    dims = size(sys.dipoles)[1:3]
    spin_infos = [i => Moment(s=(N-1)/2, g=1.0) for (i, N) in enumerate(Ns_contracted)]  # TODO: Decisions about g-factor
    sys_entangled = System(contracted_crystal, spin_infos, :SUN; dims)

    # The (uniform) external field couples to the total moment of each unit. The
    # unit g-factor is the identity, so the field is stored directly per unit.
    # Non-uniform fields may be applied later via `set_field!`/`set_field_at!`.
    fill!(sys_entangled.extfield, B)

    # Transfer rng from origin system to entangled system
    copy!(sys_entangled.rng, sys.rng)

    # TODO: Extend to inhomogenous systems
    # For each contracted site, scan original interactions and reconstruct as necessary.
    new_pair_data = Tuple{Bond, Matrix{ComplexF64}}[]
    for (contracted_site, N) in zip(1:natoms(contracted_crystal), Ns_contracted)
        Ns = Ns_unit[contracted_site]

        ## Onsite portion of interaction
        relevant_sites = atoms_in_unit(unit_map, contracted_site)
        unit_operator = zeros(ComplexF64, N, N)

        # Pair couplings that become onsite (intra-unit) couplings
        original_interactions = sys.interactions_union[relevant_sites]
        for (site, interaction) in zip(relevant_sites, original_interactions)
            onsite_original = interaction.onsite
            unit_index = unit_map.atom_to_unit[site][2]
            unit_operator += local_op_to_product_space(onsite_original, unit_index, Ns)
        end

        # Sort all PairCouplings in couplings that will be within a unit and couplings that will be between units
        pcs_intra = PairCoupling[]
        pcs_inter = PairCoupling[]
        for interaction in original_interactions, pc in interaction.pair
            (; bond) = pc
            if bond_is_in_unit(bond, unit_map)
                push!(pcs_intra, pc)
            else
                push!(pcs_inter, pc)
            end
        end

        # Convert intra-unit PairCouplings to onsite couplings
        for pc in pcs_intra
            accum_pair_coupling_into_bond_operator_in_unit!(unit_operator, pc, sys, contracted_site, unit_map)
        end
        set_onsite_coupling!(sys_entangled, unit_operator, contracted_site)

        ## Convert inter-unit PairCouplings into new pair couplings
        for pc in pcs_inter
            (; newbond, bond_operator) = pair_coupling_into_bond_operator_between_units(pc, sys, unit_map)
            push!(new_pair_data, (newbond, bond_operator))
        end
    end

    # Now have list of bonds and bond operators. First we must find individual
    # exemplars of each symmetry class of bonds in terms of the *new* crystal.
    all_bonds_with_interactions = [data[1] for data in new_pair_data]
    exemplars = Bond[]
    while length(all_bonds_with_interactions) > 0
        exemplar = all_bonds_with_interactions[1]
        all_bonds_with_interactions = filter(all_bonds_with_interactions) do bond
            !is_related_by_symmetry(contracted_crystal, bond, exemplar)
        end
        push!(exemplars, exemplar)
    end

    # Sum bond operators associated with exemplar and set the interaction. Use
    # `extract_parts=false` because dipoles are NaNs for entangled systems.
    for bond in exemplars
        relevant_interactions = filter(data -> data[1] == bond, new_pair_data)
        bond_operator = sum(data[2] for data in relevant_interactions)
        set_pair_coupling!(sys_entangled, bond_operator, bond; extract_parts=false)
    end

    # Clone the physical (bare) system and derive all entanglement metadata from
    # the immutable `units` truth. The contracted `sys_entangled` carries `origin`
    # = its chemical cell (set by the `System` constructor), matching the physical
    # `bare_system`, so both reshape in tandem through the ordinary backbone.
    bare_system = clone_system(sys)
    rebuild_entanglement!(sys_entangled, bare_system, units)

    # Initialize coherent states of each entangled unit to the tensor product of
    # bare coherent states.
    for unit in eachsite(sys_entangled)
        members = entangled_unit_members(sys_entangled, unit)
        Z = kron((bare_system.coherents[bs] for bs in members)...)
        set_coherent!(sys_entangled, Z, unit)
    end

    return sys_entangled
end

function entangle_system(::System{0}, _)
    error("Cannot entangle units of a `:dipole`-mode `System`. Use `:SUN` mode.")
end


################################################################################
# Per-unit dipole operators
################################################################################
# Build the product-space spin operators for each atom of the physical
# `bare_system`, embedded into the local Hilbert space of the containing
# entangled unit. No g-tensor is applied. An entangled system is homogeneous, so
# these depend only on the atom index. These per-atom operators are used to
# populate `bare_system.dipoles`.
function build_bare_dipole_operators(bare_system, unit_map)
    Ns_unit = Ns_in_units(bare_system, unit_map)
    natoms = length(unit_map.atom_to_unit)
    return map(1:natoms) do atom
        S_local = spin_matrices_of_dim(; N=bare_system.Ns[1, 1, 1, atom])
        unit, k = unit_map.atom_to_unit[atom]
        ntuple(α -> local_op_to_product_space(S_local[α], k, Ns_unit[unit]), 3)
    end
end

# Populate the entanglement metadata of a contracted `sys` (which may be
# reshaped), given the physical `bare_system` reshaped to the same geometry and
# the immutable `units` truth (a grouping of chemical-cell atoms of the
# *unreshaped* physical crystal). The derived `unit_map` and the cached
# operator tables are rebuilt here geometrically. This is the single entry point
# used by both `entangle_system` and every reshaping variant.
#
# The mapping is derived from positions in the shared lattice frame (contracted
# and physical crystals share `latvecs`); each member's `cell_offset` records the
# cell of its atom relative to the unit's cell, so a unit's atoms may lie in
# different cells of the reshaped system (a "straddling" unit).
function rebuild_entanglement!(sys::System, bare_system, units::Vector{UnitSpec})
    bare_origin = something(bare_system.origin, bare_system)

    # Chemical-cell contraction: per chemical atom, its within-unit slot `k` and
    # global displacement from the unit center. Reshape-invariant.
    _, um_chem = contract_crystal(bare_origin.crystal, units)
    latvecs = bare_origin.crystal.latvecs
    goff = map(um_chem.atom_to_unit) do (unit, k)
        latvecs * um_chem.unit_members[unit][k].offset
    end
    slot = [k for (_, k) in um_chem.atom_to_unit]

    # Build the reshaped, homogeneous `unit_map` sized to the reshaped
    # physical crystal. For each atom `a` of the reshaped chemical cell, locate
    # the contracted site sitting at its unit center (which may fall in a
    # neighboring cell ⇒ straddling). `inverse[u]` is ordered by slot to align
    # with the product-space factor order.
    nunits = natoms(sys.crystal)
    atoms_per_unit = length(first(units))
    forward = Vector{Tuple{Int, Int}}(undef, natoms(bare_system.crystal))
    inverse = [Vector{UnitMember}(undef, atoms_per_unit) for _ in 1:nunits]
    for a in 1:natoms(bare_system.crystal)
        a_chem = map_atom_to_other_crystal(bare_system.crystal, a, bare_origin.crystal)
        k = slot[a_chem]
        center = global_position_at(bare_system, CartesianIndex(1, 1, 1, a)) - goff[a_chem]
        unit_site = position_to_site(sys, orig_crystal(sys).latvecs \ center)
        u = to_atom(unit_site)
        forward[a] = (u, k)
        # Physical offset from the unit center (reshaped fractional frame; used
        # by `entangled_measure`) and the integer cell of atom `a` relative to
        # the unit's cell (nonzero ⇒ straddling; used by `member_site`).
        offset = bare_system.crystal.latvecs \ goff[a_chem]
        cell_offset = SVector(1, 1, 1) .- SVector(to_cell(unit_site))
        inverse[u][k] = UnitMember(a, offset, cell_offset)
    end
    unit_map = UnitMap(forward, inverse)

    bare_dipole_operators = build_bare_dipole_operators(bare_system, unit_map)
    sys.entanglement = Entanglement(bare_system, units, unit_map, bare_dipole_operators)
    return sys
end


# Entangled path of `set_params!`. The labeled parameters live on the physical
# (bare) system; the contracted couplings are a derived quantity that must be
# regenerated. The entanglement *metadata* (contraction mapping, per-unit dipole
# operators) is purely geometric and unaffected by coupling values, so only the
# contracted interactions are rebuilt here.
function set_params_entangled!(sys::System, labels, vals)
    bare_system = get_entanglement(sys).bare_system
    units = get_entanglement(sys).units

    # 1. Update the labels on the bare system (the source of truth). This takes
    #    the ordinary path and also syncs `bare_system.origin`.
    set_params!(bare_system, labels, vals)

    # 2. Regenerate the contracted chemical-cell couplings from the updated bare
    #    chemical cell. `sys_chem` shares the crystal/structure of `sys`'s own
    #    contracted chemical cell (`something(sys.origin, sys)`); only coupling
    #    values differ.
    bare_origin = something(bare_system.origin, bare_system)
    sys_chem = entangle_system(bare_origin, units)

    # 3. Install the regenerated contracted couplings. For an unreshaped `sys`,
    #    write directly; for a reshaped `sys`, update its chemical-cell `origin`
    #    and re-transfer to the reshaped geometry (the ordinary backbone).
    target = something(sys.origin, sys)
    target.params = sys_chem.params
    target.interactions_union = sys_chem.interactions_union
    isnothing(sys.origin) || transfer_params_from_origin!(sys)
    return
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
# (not of the `bare_system`).
function accum_bare_field_grad_coherents!(HZ, Z::Array{CVec{N}, 4}, ent::Entanglement) where N
    (; bare_system, unit_map, bare_dipole_operators) = ent
    dims = bare_system.dims

    # Reuse the bare system's persistent scratch buffers.
    S_bare, dE_dS_bare = get_dipole_buffers(bare_system, 2)
    fill!(dE_dS_bare, zero(Vec3))

    # Physical dipoles ⟨Sₐ⟩ evaluated from the *passed* `Z` (which may be a trial
    # state not yet synced into `bare_system.dipoles`), so the dipole-dependent
    # Ewald field is consistent with the state whose gradient is requested.
    for u in CartesianIndices(Z)
        Zu = Z[u]
        for member in unit_map.unit_members[to_atom(u)]
            bs = member_site(u, member, dims)
            S = bare_dipole_operators[member.atom]
            S_bare[bs] = Vec3(real(dot(Zu, S[1], Zu)), real(dot(Zu, S[2], Zu)), real(dot(Zu, S[3], Zu)))
        end
    end

    # Per-atom dE/dSₐ on the bare system: Zeeman gₐ' B plus, if enabled, the
    # long-range dipole-dipole field (which itself includes the g-tensor).
    for bare_site in eachsite(bare_system)
        dE_dS_bare[bare_site] += bare_system.gs[bare_site]' * bare_system.extfield[bare_site]
    end
    if !isnothing(bare_system.ewald)
        accum_ewald_grad!(dE_dS_bare, S_bare, bare_system)
    end

    # Sum each physical atom's contribution into its containing unit's coherent
    # gradient. A unit's members may straddle cells; `member_site` routes them.
    for u in CartesianIndices(Z)
        Zu = Z[u]
        for member in unit_map.unit_members[to_atom(u)]
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
    (; bare_system, unit_map, bare_dipole_operators) = get_entanglement(sys)
    for member in unit_map.unit_members[to_atom(site)]
        bs = member_site(site, member, bare_system.dims)
        S = bare_dipole_operators[member.atom]
        bare_system.dipoles[bs] = Vec3(real(dot(Z, S[1], Z)), real(dot(Z, S[2], Z)), real(dot(Z, S[3], Z)))
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
    (; bare_system, unit_map) = get_entanglement(sys)

    # Atom indices (of the physical system) that lie in the specified unit.
    atoms = [member.atom for member in unit_map.unit_members[unit]]
    @assert length(dipoles) == length(atoms) "Invalid number of dipoles for specified unit."

    # Local Hilbert space dimensions for those atoms.
    Ns = Ns_in_units(bare_system, unit_map)[unit]

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
# constructed at the atom level on the physical `bare_system` (reusing the
# ordinary `ssf_custom` for g-tensor, form factors, correlation pairs, and
# combiner), then transformed into a unit-level measure via `entangled_measure`.
# `ssf_custom(f, sys::System)` dispatches here when `sys` is entangled.
function ssf_custom_entangled(f, sys::System; apply_g, formfactors)
    (; bare_system) = get_entanglement(sys)
    measure_atom = ssf_custom(f, bare_system; apply_g, formfactors)
    return entangled_measure(measure_atom, sys)
end

# Transform an atom-level MeasureSpec (indexed to `bare_system`) into a
# unit-level MeasureSpec indexed to the contracted `sys`. For each unit and each
# of its `atoms_per_unit` subsites, the atom operator is embedded into the
# product-space Hilbert space via `local_op_to_product_space`. Position offsets
# and form factors come from `unit_map`. Observables are uniform across
# cells (g-factors are uniform), so the per-cell operator is broadcast.
function entangled_measure(measure, sys::System)
    (; bare_system, unit_map) = get_entanglement(sys)
    Ns_unit = Ns_in_units(bare_system, unit_map)

    nobs = num_observables(measure)
    dims = sys.dims

    nunits = length(unit_map.unit_members)
    atoms_per_unit = length(unit_map.unit_members[1])  # uniform by construction

    Op = eltype(measure.operators)
    new_ops     = Array{Op, 6}(undef, atoms_per_unit, nobs, dims..., nunits)
    new_offsets = zeros(Vec3, atoms_per_unit, nunits)
    new_ff      = Array{FormFactor, 3}(undef, atoms_per_unit, nobs, nunits)

    for i in 1:nunits
        for (k, member) in enumerate(unit_map.unit_members[i])
            atom = member.atom  # atom index within a chemical cell of bare_system
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
