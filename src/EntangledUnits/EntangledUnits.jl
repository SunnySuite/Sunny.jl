################################################################################
# Types
################################################################################
# Data for mapping one site inside a unit back to the site of the original
# system.
struct InverseData
    site   :: Int64  # Atom index of original, uncontracted crystal
    offset :: Vec3   # Position offset of original atom relative to center of unit
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
# `contraction_info`, `unit_members`, and the operator tables all track the
# current (possibly reshaped) `bare_system`. Because the mapping is rebuilt
# geometrically from `units`, an entangled unit may straddle cell boundaries of a
# reshaped system.
struct Entanglement <: AbstractEntanglement
    bare_system           :: System                                    # Physical (reshaped) system
    units                 :: Vector{Vector{Int}}                       # TRUTH: chemical-cell atom grouping
    contraction_info      :: CrystalContractionInfo                    # Forward/inverse mapping (matches bare_system)
    unit_members          :: Array{Vector{CartesianIndex{4}}, 4}       # Physical sites comprising each unit site
    bare_dipole_operators :: Vector{NTuple{3, HermitianC64}}           # Product-space spin ops per physical atom (no g)
end

function clone_entanglement(ent::Entanglement)
    return Entanglement(clone_system(ent.bare_system), deepcopy(ent.units),
                        copy(ent.contraction_info), deepcopy(ent.unit_members),
                        ent.bare_dipole_operators)
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
        inverse_map[new_pair] = InverseData(site, Vec3(0, 0, 0))
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
            inverse_map[new_pair] = InverseData(site, offset)
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


# Reconstruct original crystal from contracted Crystal and a CrystalContractionInfo
function expand_crystal(contracted_crystal, contraction_info)
    (; forward, inverse) = contraction_info
    contracted_positions = contracted_crystal.positions
    nsites_expanded = length(forward)
    expanded_positions = [Vec3(0, 0, 0) for _ in 1:nsites_expanded]
    for (contracted_idx, original_site_data) in enumerate(inverse)
        for (; site, offset) in original_site_data
            expanded_positions[site] = contracted_positions[contracted_idx] + offset
        end
    end
    Crystal(contracted_crystal.latvecs, expanded_positions)
end


# Returns a list of length equal to the number of "units" in the a contracted
# crystal. Each list element is itself a list of integers, each of which
# corresponds to the N of the corresponding site of the original system. The
# order is consistent with that given by the `inverse` field of a
# `CyrstalContractionInfo`.
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

# Physical (bare-system) sites comprising the entangled unit at contracted
# `unit_site` of `sys`. Precomputed in `rebuild_entanglement!` (a plain lookup
# here) so that a unit may straddle cell boundaries — member sites need not
# share the unit's cell.
entangled_unit_members(sys::System, unit_site) = get_entanglement(sys).unit_members[to_cartesian(unit_site)]


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
# System contraction
################################################################################
function entangle_system(::System{0}, _)
    error("Cannot contract a dipole system. Use :SUN mode.")
end

function entangle_system(sys::System{M}, units) where M
    # Construct contracted crystal
    contracted_crystal, contraction_info = contract_crystal(sys.crystal, units)

    # Make sure we have a uniform external field
    @assert allequal(@view sys.extfield[:,:,:,:]) "Entangled units require a uniform applied field."
    B = sys.extfield[1,1,1,1]

    # Determine Ns for local Hilbert spaces (all must be equal). (TODO: Determine if alternative behavior preferable in mixed case.)
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

        # Pair interactions that become within-unit interactions
        original_interactions = sys.interactions_union[relevant_sites]
        for (site, interaction) in zip(relevant_sites, original_interactions)
            # Onsite anisotropy portion. The Zeeman term is *not* folded in here;
            # it is handled as a first-class term on the contracted system via
            # `sys.extfield` and the cached per-unit `dipole_operators` (see
            # `set_field_entangled!` and `set_energy_grad_coherents!`).
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

    # We collected all bond_operators associated with a particular exemplar, sum
    # them, and set the interaction. Use `extract_parts=false` to keep the whole
    # operator in the `general` (tensor-decomposed) channel. Extracting bilinear
    # and biquadratic parts would express them in terms of the unit's spin-(N-1)/2
    # operators, but the contracted `sys.dipoles` holds a unit's *total moment*
    # (Σₖ gₖ Sₖ), not ⟨spin_matrices_of_dim(N)⟩. Those channels would then be
    # evaluated against the wrong vector in `energy_aux`/`set_energy_grad_*`. The
    # general channel is evaluated directly from coherent states, so it is exact.
    for bond in exemplars
        relevant_interactions = filter(data -> data[1] == bond, new_pair_data)
        bond_operator = sum(data[2] for data in relevant_interactions)
        set_pair_coupling!(sys_entangled, bond_operator, bond; extract_parts=false)
    end

    return (; sys_entangled, contraction_info)
end


################################################################################
# Per-unit dipole operators
################################################################################
# Build the product-space spin operators for each atom of the physical
# `bare_system`, embedded into the local Hilbert space of the containing
# entangled unit. No g-tensor is applied. An entangled system is homogeneous, so
# these depend only on the atom index. These per-atom operators feed
# `sync_unit_dipoles!` (populating `bare_system.dipoles`) and are summed with
# g-weights to form the per-unit `sys.dipole_operators` below.
function build_bare_dipole_operators(bare_system, contraction_info)
    Ns_unit = Ns_in_units(bare_system, contraction_info)
    natoms = length(contraction_info.forward)
    return map(1:natoms) do atom
        S_local = spin_matrices_of_dim(; N=bare_system.Ns[1, 1, 1, atom])
        unit, k = contraction_info.forward[atom]
        ntuple(α -> local_op_to_product_space(S_local[α], k, Ns_unit[unit]), 3)
    end
end

# Build the per-unit product-space operators `sys.dipole_operators[u]` for an
# entangled system: the g-weighted total magnetic moment of each unit,
#
#   T^α = Σ_{k ∈ u} Σ_β (g_k)_{αβ} embed(Sₖ^β),
#
# so that ⟨Z|T^α|Z⟩ is the α-component of the unit's total moment M = Σ_k g_k Sₖ.
# Derived from the per-atom `bare_dipole_operators` (one source of truth). The
# unit g-factor `sys.gs[u]` is the identity, so the per-atom g_k live *inside*
# these operators. Used by the Zeeman energy/gradient and SWT (see
# `set_energy_grad_coherents!` and `swt_data`).
function build_unit_dipole_operators(bare_system, contraction_info, bare_dipole_operators)
    Ns_unit = Ns_in_units(bare_system, contraction_info)
    nunits = length(contraction_info.inverse)
    return map(1:nunits) do unit
        Nprod = prod(Ns_unit[unit])
        T = ntuple(_ -> zeros(ComplexF64, Nprod, Nprod), 3)
        for id in contraction_info.inverse[unit]
            atom = id.site
            g = bare_system.gs[1, 1, 1, atom]
            S_embed = bare_dipole_operators[atom]
            for α in 1:3, β in 1:3
                T[α] .+= g[α, β] .* S_embed[β]
            end
        end
        # Real g-weighted sum of Hermitian embedded spins ⇒ Hermitian.
        return ntuple(α -> Hermitian(T[α]), 3)
    end
end

# Populate the entanglement metadata of a contracted `sys` (which may be
# reshaped), given the physical `bare_system` reshaped to the same geometry and
# the immutable `units` truth (a grouping of chemical-cell atoms of the
# *unreshaped* physical crystal). All derived data — the site mapping
# (`unit_members`), the contraction bookkeeping (`contraction_info`), and the
# cached operator tables — is rebuilt here geometrically. This is the single
# entry point used by both `entangle_units` and every reshaping variant.
#
# The mapping is derived purely from positions in the shared lattice frame
# (contracted and physical crystals share `latvecs`), so a unit's member atoms
# may lie in different cells of the reshaped system (a "straddling" unit).
function rebuild_entanglement!(sys::System, bare_system, units)
    bare_origin = something(bare_system.origin, bare_system)

    # Chemical-cell contraction: gives, per chemical atom, its slot `k` within
    # its unit and the physical offset from the unit center. Reshape-invariant.
    _, ci_chem = contract_crystal(bare_origin.crystal, units)
    latvecs = bare_origin.crystal.latvecs
    natoms_chem = length(ci_chem.forward)
    # Global displacement of each chemical atom from its unit center, and slot.
    goff = map(1:natoms_chem) do a
        (unit, k) = ci_chem.forward[a]
        latvecs * ci_chem.inverse[unit][k].offset
    end
    slot = [ci_chem.forward[a][2] for a in 1:natoms_chem]

    # Map each physical site of the (reshaped) `bare_system` to the full unit
    # site of the (reshaped) contracted `sys` that owns it. Straddle-safe: the
    # unit center is found by subtracting the atom's physical offset and locating
    # the contracted site there, with cell wrapping handled by `position_to_site`.
    source_idcs = Array{CartesianIndex{4}}(undef, size(eachsite(bare_system)))
    for bs in eachsite(bare_system)
        a_chem = map_atom_to_other_crystal(bare_system.crystal, to_atom(bs), bare_origin.crystal)
        center = global_position_at(bare_system, bs) - goff[a_chem]
        source_idcs[bs] = position_to_site(sys, orig_crystal(sys).latvecs \ center)
    end

    # Invert the site map into `unit_members[u] = [physical sites in unit u]`,
    # ordered by within-unit slot. This is the concrete lookup used by the sync
    # and field-mirroring hot paths (`sync_entangled_unit!`, `set_field_at!`).
    nunits = natoms(sys.crystal)
    atoms_per_unit = length(first(units))
    unit_members = Array{Vector{CartesianIndex{4}}}(undef, size(eachsite(sys)))
    for u in eachsite(sys)
        unit_members[u] = Vector{CartesianIndex{4}}(undef, atoms_per_unit)
    end
    for bs in eachsite(bare_system)
        a_chem = map_atom_to_other_crystal(bare_system.crystal, to_atom(bs), bare_origin.crystal)
        unit_members[source_idcs[bs]][slot[a_chem]] = bs
    end

    # Build the reshaped `contraction_info` (homogeneous; sized to the reshaped
    # physical crystal's atoms). `forward[a] = (unit atom, slot)`; `inverse[u]` is
    # ordered by slot so it aligns with the product-space factor order.
    forward = Vector{Tuple{Int, Int}}(undef, natoms(bare_system.crystal))
    inverse = [Vector{InverseData}(undef, atoms_per_unit) for _ in 1:nunits]
    for a in 1:natoms(bare_system.crystal)
        a_chem = map_atom_to_other_crystal(bare_system.crystal, a, bare_origin.crystal)
        k = slot[a_chem]
        unit_atom = to_atom(source_idcs[CartesianIndex(1, 1, 1, a)])
        forward[a] = (unit_atom, k)
        # Physical offset from the unit center, in the reshaped lattice frame.
        # (For a straddling atom this may fall outside the [0,1) cell — that is
        # the true physical displacement and is what `expand_crystal` /
        # `entangled_measure` expect.)
        offset = bare_system.crystal.latvecs \ goff[a_chem]
        inverse[unit_atom][k] = InverseData(a, offset)
    end
    contraction_info = CrystalContractionInfo(forward, inverse)

    bare_dipole_operators = build_bare_dipole_operators(bare_system, contraction_info)
    sys.dipole_operators = build_unit_dipole_operators(bare_system, contraction_info, bare_dipole_operators)
    sys.entanglement = Entanglement(bare_system, units, contraction_info, unit_members, bare_dipole_operators)
    return sys
end


################################################################################
# Public API
################################################################################
"""
    entangle_units(sys::System{N}, units)

Create a new [`System`](@ref) of "entangled units" from an existing `System`.
`units` is a list of tuples specifying the atoms inside each unit cell that will
be grouped into a single "entangled unit." This feature is only supported for
systems that can be viewed as a regular lattice of a single unit type (all
dimers, all trimers, etc). Sunny will use the SU(_N_) formalism to model each
one of these units as a distinct Hilbert space in which the full quantum
mechanical structure is locally preserved.

Interactions must be specified for the original `System`. Sunny will
automatically reconstruct the appropriate interactions for the entangled system.
The returned `System` reports physical geometry (positions, dipoles) against the
original crystal, while its dynamical variables are the coherent states of the
entangled units.
"""
function entangle_units(sys::System{N}, units) where {N}
    # External field is folded into the onsite interactions of the contracted
    # system, which is homogeneous (indexed by unit, not site). So g-factors must
    # be uniform across unit cells.
    for atom in axes(sys.coherents, 4)
        @assert allequal(@view sys.gs[:,:,:,atom]) "Entangled units require g-factors be uniform across unit cells"
    end

    # Generate the contracted system and the physical (bare) system. The
    # contracted system carries `origin` = its chemical cell (set by the `System`
    # constructor), matching the physical `bare_system`, so both reshape in
    # tandem through the ordinary backbone.
    sys_entangled = entangle_system(sys, units).sys_entangled
    bare_system = clone_system(sys)

    # Store the immutable `units` truth and derive all mapping metadata from it.
    units_truth = [collect(Int, u) for u in units]
    rebuild_entanglement!(sys_entangled, bare_system, units_truth)

    # Coordinate the contracted coherent states and the physical dipoles.
    set_expected_dipoles!(sys_entangled)

    return sys_entangled
end

function entangle_units(::System{0}, _)
    error("Cannot entangle units of a `:dipole`-mode `System`. Use `:SUN` mode.")
end


################################################################################
# Dynamics: syncing physical dipoles with coherent states
################################################################################
# An entangled `System` *is* the contracted system: `eachsite`, `sys.dipoles`,
# `sys.crystal` are unit-level, while physical data lives in
# `sys.entanglement.bare_system`. Dynamics (`step!`, `randomize_spins!`,
# `minimize_energy!`, `set_coherent!`, `suggest_timestep`, `clone_system`) work
# on the contracted system directly and self-sync physical dipoles (see
# `setspin!`/`set_expected_dipoles!`), so no entangled-specific overrides are
# needed for them. The methods below cover the cases that need physical data or
# entangled-specific behavior.
#
# NB: for a unified entangled `System`, `eachsite(sys)` iterates the entangled
# *units* (the native sites of the contracted system). To iterate the physical
# atoms, iterate `eachsite(sys.entanglement.bare_system)`.

# Dipole-dipole energy for an entangled system, evaluated on the physical bare
# system (whose dipoles are synced with the coherent states).
function entangled_ewald_energy(ent::Entanglement)
    bare = ent.bare_system
    return isnothing(bare.ewald) ? 0.0 : ewald_energy(bare)
end

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
    (; bare_system, unit_members, bare_dipole_operators) = ent

    # Reuse the bare system's persistent scratch buffers.
    S_bare, dE_dS_bare = get_dipole_buffers(bare_system, 2)
    fill!(dE_dS_bare, zero(Vec3))

    # Physical dipoles ⟨Sₐ⟩ evaluated from the *passed* `Z` (which may be a trial
    # state not yet synced into `bare_system.dipoles`), so the dipole-dependent
    # Ewald field is consistent with the state whose gradient is requested.
    for u in CartesianIndices(unit_members)
        Zu = Z[u]
        for bs in unit_members[u]
            S = bare_dipole_operators[to_atom(bs)]
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
    # gradient. A unit's members may straddle cells; `unit_members` routes them.
    for u in CartesianIndices(unit_members)
        Zu = Z[u]
        for bs in unit_members[u]
            S = bare_dipole_operators[to_atom(bs)]
            h = dE_dS_bare[bs]
            for β in 1:3
                Sβ = SMatrix{N, N}(S[β])
                HZ[u] += h[β] * (Sβ * Zu)
            end
        end
    end
    return
end

# Sync the physical dipoles of the unit at contracted site `unit_site` from its
# coherent state `Zu`, writing into `bare_system.dipoles`, and return the unit's
# total magnetic moment M = Σₖ gₖ Sₖ over its physical atoms.
function sync_unit_dipoles!(ent::Entanglement, unit_site, Zu)
    (; bare_system, unit_members, bare_dipole_operators) = ent
    M = zero(Vec3)
    for bs in unit_members[unit_site]
        S = bare_dipole_operators[to_atom(bs)]  # per-atom product-space spin operators (no g)
        d = Vec3(real(dot(Zu, S[1], Zu)), real(dot(Zu, S[2], Zu)), real(dot(Zu, S[3], Zu)))
        bare_system.dipoles[bs] = d
        M += bare_system.gs[bs] * d
    end
    return M
end

# Keep the cached dipoles of a single unit coherent with its coherent state,
# given a unit `site` of the contracted `sys`. Two things are updated:
#
#   1. The physical dipoles of every atom in the unit (`bare_system.dipoles`).
#   2. The contracted `sys.dipoles[site]`, redefined to carry a meaningful value:
#      the unit's *total physical magnetic moment* expressed in the unit's own
#      g-convention. `expected_spin` is nonsense here (it would assume the
#      product-space dimension N=∏Nₖ is a single spin-(N-1)/2 irrep), so instead
#      we form M = Σₖ gₖ Sₖ over the physical atoms and store `sys.gs[site] \ M`.
#      With the conventional `sys.gs[site] = I`, this is just M. Then
#      `sys.gs[site] * sys.dipoles[site] = M` recovers the total moment, which
#      will facilitate a Zeeman coupling on the entangled unit.
#
# Called after a targeted coherent-state update (`set_coherent!`/`setspin!`).
function sync_entangled_unit!(sys::System, site)
    ent = get_entanglement(sys)
    site = to_cartesian(site)
    M = sync_unit_dipoles!(ent, site, sys.coherents[site])
    sys.dipoles[site] = sys.gs[site] \ M
    return
end

# Entangled path of `set_expected_dipoles!`: refresh every unit's cached dipoles
# (physical + contracted total moment) from the coherent states.
function set_expected_dipoles_entangled!(sys::System)
    ent = get_entanglement(sys)
    for site in eachsite(sys)
        M = sync_unit_dipoles!(ent, site, sys.coherents[site])
        sys.dipoles[site] = sys.gs[site] \ M
    end
    return
end


################################################################################
# Coherent states and measurements
################################################################################
# Find the unique coherent state corresponding to a set of fully-polarized
# dipoles on each site inside a specified entangled unit.
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
