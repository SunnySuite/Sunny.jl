################################################################################
# Types
################################################################################
# Data for mapping one site inside a unit back to the site of the original
# system.
struct InverseData
    site        :: Int64          # Atom index of original, uncontracted crystal
    offset      :: Vec3           # Position offset of original atom relative to center of unit
    cell_offset :: SVector{3,Int} # Cell of this member relative to the unit's cell (nonzero for a straddling unit)
end

# `forward` contains a list from sites of the original crystal to a site of the
# contracted crystal, including an extra index to keep track entangled units:
# `(contracted_crystal_site_index, intra_unit_site_index)`. If the first index
# refers to a site in the new crystal that does not contain multiple units, than
# the second index will always be 1.
#
# `inverse` contains a list of length equal to the number of sites in the
# contracted crystal (corresponding to `contracted_crystal_site_index` above).
# Each element of this list is another list of tuples,
# `(site_of_original_crystal, position_offset)`. The position offset is applied
# to the position of the contracted crystal to recover the corresponding
# location in the original crystal. The length of these sublists corresponds to
# the number of sites within the entangled unit.
struct CrystalContractionInfo
    forward :: Vector{Tuple{Int64, Int64}}  # Original site index -> full unit index (contracted crystal site index and unit subindex)
    inverse :: Vector{Vector{InverseData}}  # List ordered according to contracted crystal sites. Each element is itself a list containing original crystal site indices and corresponding offset information
end

function Base.copy(cci::CrystalContractionInfo)
    return CrystalContractionInfo(copy(cci.forward), copy(cci.inverse))
end

# Metadata for a System with entangled units. Hilbert dimension N is _not_
# carried in the bare system's type. This avoids unecessary dynamic lookups.
#
# `units` is the immutable "truth": the grouping of chemical-cell atoms into
# entangled units, expressed against the *unreshaped* physical crystal
# (`something(bare_system.origin, bare_system)`). Everything else is derived and
# rebuilt by `rebuild_entanglement!` whenever the system is reshaped —
# `contraction_info` and the operator tables track the current (possibly
# reshaped) `bare_system`. The physical sites comprising a unit are recovered by
# integer arithmetic from `contraction_info` (see `member_site`); each member
# carries a `cell_offset` so a unit may straddle cell boundaries of a reshaped
# system.
struct Entanglement <: AbstractEntanglement
    bare_system           :: System                                    # Physical system
    units                 :: Vector{Vector{Int}}                       # TRUTH: chemical-cell atom grouping
    contraction_info      :: CrystalContractionInfo                    # Forward/inverse mapping (matches bare_system)
    bare_dipole_operators :: Vector{NTuple{3, HermitianC64}}           # Product-space spin ops per physical atom (no g)
end

function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.bare_system), deepcopy(ent.units),
                        copy(ent.contraction_info), ent.bare_dipole_operators)
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
# Takes a crystal and a list of integer tuples. The tuples indicate which sites
# in the original crystal are to be grouped, i.e., contracted into a single site
# of a new crystal.
function contract_crystal(crystal, units)

    # Determine which sites are entangled and which are not.
    unentangled_sites = 1:natoms(crystal)
    entangled_sites = Int64[]
    for unit in units, site in unit
        push!(entangled_sites, site)
        unentangled_sites = filter(!=(site), unentangled_sites)
    end

    # Check that unit definitions are valid.
    @assert length(entangled_sites) + length(unentangled_sites) == natoms(crystal) "Invalid entangled unit specification."
    @assert isempty(intersect(entangled_sites, unentangled_sites)) "Invalid entangled unit specification." # Sanity check. Should be true by construction. Remove later

    # Prepare mapping dictionaries and initialize iteration variable.
    forward_map = Dict{Int64, Tuple{Int64, Int64}}()
    inverse_map = Dict{Tuple{Int64, Int64}, InverseData}()
    new_site_current = 1

    # Add sites that are *not* encapsulated in any unit as individual sites in
    # the new crystal with no associated displacement.
    new_positions = []
    for site in unentangled_sites
        push!(new_positions, crystal.positions[site])
        new_pair = (new_site_current, 1)
        forward_map[site] = new_pair
        inverse_map[new_pair] = InverseData(site, Vec3(0, 0, 0), SVector(0, 0, 0))
        new_site_current += 1
    end

    # Assign entangled units to single site in new crystal and record mapping
    # information.
    for unit in units
        # Find new position by averaging location of entangled positions.
        old_positions = [crystal.positions[i] for i in unit]
        new_position = sum(old_positions) / length(old_positions)
        push!(new_positions, new_position)

        # Record forward and inverse mapping information, including
        # the displacement data from the new unit position.
        for (j, site) in enumerate(unit)
            new_pair = (new_site_current, j)
            offset = crystal.positions[site] - new_position

            forward_map[site] = new_pair
            inverse_map[new_pair] = InverseData(site, offset, SVector(0, 0, 0))
        end

        new_site_current += 1
    end

    nsites_new = new_site_current - 1

    # Sorted list version of forward map
    forward_keys = collect(keys(forward_map))
    forward_vals = collect(values(forward_map))
    idcs = sortperm(forward_keys)
    forward_list = [forward_vals[i] for i in idcs]

    # Sorted list version of inverse map
    inverse_keys = collect(keys(inverse_map))
    inverse_vals = collect(values(inverse_map))
    idcs = sortperm(inverse_keys)
    inverse_keys = inverse_keys[idcs]
    inverse_vals = inverse_vals[idcs]

    inverse_list = [InverseData[] for _ in 1:nsites_new]
    for (n, key) in enumerate(inverse_keys)
        new_site, _ = key  # `key` is a tuple (global_site, local_site_in_unit)
        push!(inverse_list[new_site], inverse_vals[n])
    end

    # Generate a new contracted crystal and information to invert contraction.
    # Space group must be set to 1 to allow "invalid" anisotropies -- these are
    # generated by PairCouplings that become onsite couplings. NB: Add option to
    # override to unconventional cell warnings and errors? (Probably yes)
    new_crystal = Crystal(crystal.latvecs, new_positions, 1)
    contraction_info = CrystalContractionInfo(forward_list, inverse_list)

    return new_crystal, contraction_info
end


# Returns a list of length equal to the number of "units" in the a contracted
# crystal. Each list element is itself a list of integers, each of which
# corresponds to the N of the corresponding site of the original system. The
# order is consistent with that given by the `inverse` field of a
# `CrystalContractionInfo`.
function Ns_in_units(sys_original, contraction_info)
    Ns = [Int64[] for _ in 1:length(contraction_info.inverse)]
    for (n, contracted_sites) in enumerate(contraction_info.inverse)
        for (; site) in contracted_sites
            push!(Ns[n], sys_original.Ns[site])
        end
    end
    Ns
end


# Pull out original indices of sites in entangled unit
atoms_in_unit(contraction_info, i) = [inverse_data.site for inverse_data in contraction_info.inverse[i]]

# Physical (bare-system) site of a unit member, given the containing unit's site
# and the member's `InverseData`. The member atom sits in cell `cell_offset`
# relative to the unit's cell (periodically wrapped), so a unit may straddle
# cell boundaries. Mirrors `bonded_site`.
@inline function member_site(unit_site, id::InverseData, dims)
    CartesianIndex(altmod1.(to_cell(unit_site) .+ Tuple(id.cell_offset), dims)..., id.site)
end

# Physical (bare-system) sites comprising the entangled unit at contracted
# `unit_site` of `sys`.
function entangled_unit_members(sys::System, unit_site)
    (; bare_system, contraction_info) = get_entanglement(sys)
    unit_site = to_cartesian(unit_site)
    return [member_site(unit_site, id, bare_system.dims) for id in contraction_info.inverse[to_atom(unit_site)]]
end


################################################################################
# Pair-coupling contraction
################################################################################
# Checks whether a bond given in terms of the original crystal lies inside a
# single unit of the contracted system.
function bond_is_in_unit(bond::Bond, ci::CrystalContractionInfo)
    (ci.forward[bond.i][1] == ci.forward[bond.j][1]) && (bond.n == [0, 0, 0])
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
function accum_pair_coupling_into_bond_operator_in_unit!(op, pc, sys, contracted_site, contraction_info)
    pc.isculled && return

    Ns_unit = Ns_in_units(sys, contraction_info)[contracted_site]
    i, j = pc.bond.i, pc.bond.j
    i_unit = contraction_info.forward[i][2]
    j_unit = contraction_info.forward[j][2]
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
function pair_coupling_into_bond_operator_between_units(pc, sys, contraction_info)
    (; i, j, n) = pc.bond
    unit1, unitsub1 = contraction_info.forward[i]
    unit2, unitsub2 = contraction_info.forward[j]

    Ns_local = Ns_in_units(sys, contraction_info)
    Ns1 = Ns_local[unit1]
    Ns2 = Ns_local[unit2]
    N1 = sys.Ns[1, 1, 1, i]
    N2 = sys.Ns[1, 1, 1, j]
    dim1 = prod(Ns1)
    dim2 = prod(Ns2)
    N = dim1 * dim2

    # Site i acts on unit1 (first tensor factor), site j on unit2 (second).
    embed_i = A -> kron(local_op_to_product_space(A, unitsub1, Ns1), I(dim2))
    embed_j = B -> kron(I(dim1), local_op_to_product_space(B, unitsub2, Ns2))

    bond_operator = zeros(ComplexF64, N, N)
    accum_bond_operator!(bond_operator, pc, embed_i, embed_j, N1, N2)
    return (; newbond=Bond(unit1, unit2, n), bond_operator)
end


################################################################################
# Public API
################################################################################
"""
    entangle_system(sys::System{N}, groups::Vector{Vector{Int}})

Create a new [`System`](@ref), where `groups` of atom indices are collected into
"entangled units". Sunny will model each such unit within a tensor-product
Hilbert space to allow for local quantum entanglement. This feature currently
requires a single unit type (all dimers, all trimers, etc).

The input `sys` should be in `:SUN` mode and have 1×1×1 dimensions. All
interactions in `sys` will be transferred to the entangled system. Subsequent
reshapings of this entangled system are then allowed.
"""
function entangle_system(sys::System{M}, groups) where M
    isnothing(sys.origin) || error("Entangle a single-cell system first, then reshape")

    # Construct contracted crystal
    contracted_crystal, contraction_info = contract_crystal(sys.crystal, groups)

    # Make sure we have a uniform external field
    @assert allequal(@view sys.extfield[:,:,:,:]) "Entangled units require a uniform applied field."
    B = sys.extfield[1,1,1,1]

    # Determine Ns for local Hilbert spaces (all must be equal). (TODO:
    # Determine if alternative behavior preferable in mixed case.)
    Ns_unit = Ns_in_units(sys, contraction_info)
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
        relevant_sites = atoms_in_unit(contraction_info, contracted_site)
        unit_operator = zeros(ComplexF64, N, N)

        # Pair couplings that become onsite (intra-unit) couplings
        original_interactions = sys.interactions_union[relevant_sites]
        for (site, interaction) in zip(relevant_sites, original_interactions)
            onsite_original = interaction.onsite
            unit_index = contraction_info.forward[site][2]
            unit_operator += local_op_to_product_space(onsite_original, unit_index, Ns)
        end

        # Sort all PairCouplings in couplings that will be within a unit and couplings that will be between units
        pcs_intra = PairCoupling[]
        pcs_inter = PairCoupling[]
        for interaction in original_interactions, pc in interaction.pair
            (; bond) = pc
            if bond_is_in_unit(bond, contraction_info)
                push!(pcs_intra, pc)
            else
                push!(pcs_inter, pc)
            end
        end

        # Convert intra-unit PairCouplings to onsite couplings
        for pc in pcs_intra
            accum_pair_coupling_into_bond_operator_in_unit!(unit_operator, pc, sys, contracted_site, contraction_info)
        end
        set_onsite_coupling!(sys_entangled, unit_operator, contracted_site)

        ## Convert inter-unit PairCouplings into new pair couplings
        for pc in pcs_inter
            (; newbond, bond_operator) = pair_coupling_into_bond_operator_between_units(pc, sys, contraction_info)
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
    units_truth = [collect(Int, u) for u in groups]
    rebuild_entanglement!(sys_entangled, bare_system, units_truth)

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
function build_bare_dipole_operators(bare_system, contraction_info)
    Ns_unit = Ns_in_units(bare_system, contraction_info)
    natoms = length(contraction_info.forward)
    return map(1:natoms) do atom
        S_local = spin_matrices_of_dim(; N=bare_system.Ns[1, 1, 1, atom])
        unit, k = contraction_info.forward[atom]
        ntuple(α -> local_op_to_product_space(S_local[α], k, Ns_unit[unit]), 3)
    end
end

# Populate the entanglement metadata of a contracted `sys` (which may be
# reshaped), given the physical `bare_system` reshaped to the same geometry and
# the immutable `units` truth (a grouping of chemical-cell atoms of the
# *unreshaped* physical crystal). The derived `contraction_info` and the cached
# operator tables are rebuilt here geometrically. This is the single entry point
# used by both `entangle_system` and every reshaping variant.
#
# The mapping is derived from positions in the shared lattice frame (contracted
# and physical crystals share `latvecs`); each member's `cell_offset` records the
# cell of its atom relative to the unit's cell, so a unit's atoms may lie in
# different cells of the reshaped system (a "straddling" unit).
function rebuild_entanglement!(sys::System, bare_system, units)
    bare_origin = something(bare_system.origin, bare_system)

    # Chemical-cell contraction: per chemical atom, its within-unit slot `k` and
    # global displacement from the unit center. Reshape-invariant.
    _, ci_chem = contract_crystal(bare_origin.crystal, units)
    latvecs = bare_origin.crystal.latvecs
    goff = map(ci_chem.forward) do (unit, k)
        latvecs * ci_chem.inverse[unit][k].offset
    end
    slot = [k for (_, k) in ci_chem.forward]

    # Build the reshaped, homogeneous `contraction_info` sized to the reshaped
    # physical crystal. For each atom `a` of the reshaped chemical cell, locate
    # the contracted site sitting at its unit center (which may fall in a
    # neighboring cell ⇒ straddling). `inverse[u]` is ordered by slot to align
    # with the product-space factor order.
    nunits = natoms(sys.crystal)
    atoms_per_unit = length(first(units))
    forward = Vector{Tuple{Int, Int}}(undef, natoms(bare_system.crystal))
    inverse = [Vector{InverseData}(undef, atoms_per_unit) for _ in 1:nunits]
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
        inverse[u][k] = InverseData(a, offset, cell_offset)
    end
    contraction_info = CrystalContractionInfo(forward, inverse)

    bare_dipole_operators = build_bare_dipole_operators(bare_system, contraction_info)
    sys.entanglement = Entanglement(bare_system, units, contraction_info, bare_dipole_operators)
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
    (; bare_system, contraction_info, bare_dipole_operators) = ent
    dims = bare_system.dims

    # Reuse the bare system's persistent scratch buffers.
    S_bare, dE_dS_bare = get_dipole_buffers(bare_system, 2)
    fill!(dE_dS_bare, zero(Vec3))

    # Physical dipoles ⟨Sₐ⟩ evaluated from the *passed* `Z` (which may be a trial
    # state not yet synced into `bare_system.dipoles`), so the dipole-dependent
    # Ewald field is consistent with the state whose gradient is requested.
    for u in CartesianIndices(Z)
        Zu = Z[u]
        for id in contraction_info.inverse[to_atom(u)]
            bs = member_site(u, id, dims)
            S = bare_dipole_operators[id.site]
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
        for id in contraction_info.inverse[to_atom(u)]
            bs = member_site(u, id, dims)
            S = bare_dipole_operators[id.site]
            h = dE_dS_bare[bs]
            for β in 1:3
                Sβ = SMatrix{N, N}(S[β])
                HZ[u] += h[β] * (Sβ * Zu)
            end
        end
    end
    return
end

# Sync the bare dipoles with the coherent state of an unit at `site`
function sync_bare_dipoles_at!(sys::System, site)
    Z = sys.coherents[site]
    (; bare_system, contraction_info, bare_dipole_operators) = get_entanglement(sys)
    for id in contraction_info.inverse[to_atom(site)]
        bs = member_site(site, id, bare_system.dims)
        S = bare_dipole_operators[id.site]
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
    (; bare_system, contraction_info) = get_entanglement(sys)

    # Atom indices (of the physical system) that lie in the specified unit.
    atoms = [id.site for id in contraction_info.inverse[unit]]
    @assert length(dipoles) == length(atoms) "Invalid number of dipoles for specified unit."

    # Local Hilbert space dimensions for those atoms.
    Ns = Ns_in_units(bare_system, contraction_info)[unit]

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
# and form factors come from `contraction_info`. Observables are uniform across
# cells (g-factors are uniform), so the per-cell operator is broadcast.
function entangled_measure(measure, sys::System)
    (; bare_system, contraction_info) = get_entanglement(sys)
    Ns_unit = Ns_in_units(bare_system, contraction_info)

    nobs = num_observables(measure)
    dims = sys.dims

    nunits = length(contraction_info.inverse)
    atoms_per_unit = length(contraction_info.inverse[1])  # uniform by construction

    Op = eltype(measure.operators)
    new_ops     = Array{Op, 6}(undef, atoms_per_unit, nobs, dims..., nunits)
    new_offsets = zeros(Vec3, atoms_per_unit, nunits)
    new_ff      = Array{FormFactor, 3}(undef, atoms_per_unit, nobs, nunits)

    for i in 1:nunits
        for (k, inverse_info) in enumerate(contraction_info.inverse[i])
            atom = inverse_info.site  # atom index within a chemical cell of bare_system
            new_offsets[k, i] = inverse_info.offset
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
