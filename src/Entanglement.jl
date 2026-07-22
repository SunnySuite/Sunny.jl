################################################################################
# Types
################################################################################

struct UnitAndPartIndex
    unit :: Int             # Unit index in contracted crystal
    part :: Int             # Part index within unit
end

# One part of an entangled unit. In principle, Δcell is redundant -- it could be
# recalculated from Δpos and the unit centers.
struct UnitPart
    atom  :: Int            # Atom index in the uncontracted (bare) crystal
    Δpos  :: Vec3           # Fractional position relative to unit center
    Δcell :: SVector{3,Int} # Cell of this atom relative to the unit's cell
end

# Bidirectional mapping between entangled units and their constituent parts
struct UnitLayout
    atom_to_unit  :: Vector{UnitAndPartIndex} # bare atom → (unit, part)
    unit_to_parts :: Vector{Vector{UnitPart}} # unit → bare atom data
end

# Metadata for a System with entangled units. Hilbert dimension N of the
# uncontracted System{N} is intentionally omitted to avoid dynamic lookup.
struct Entanglement <: AbstractEntanglement
    uncontracted    :: System                           # System prior to entanglement
    layout          :: UnitLayout                       # Membership of atoms within units
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

# Bare atom indices comprising a unit, ordered to match the unit_to_parts
# field of a UnitLayout.
function atoms_in_unit(layout::UnitLayout, u)
    [atom for (; atom) in layout.unit_to_parts[u]]
end
function atoms_in_unit(sys::System, u)
    is_entangled(sys) ? atoms_in_unit(get_entanglement(sys).layout, u) : [u]
end

# Given a site representing an entangled unit, return an iterator over the
# corresponding "bare sites" of the uncontracted system.
function bare_sites_for_unit(ent::Entanglement, unit_site)
    (; uncontracted, layout) = ent
    (; dims) = uncontracted
    parts = layout.unit_to_parts[to_atom(unit_site)]
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
    return Entanglement(clone_system(ent.uncontracted), ent.layout,
                        ent.lifted_spin_ops)
end
clone_entanglement(::Nothing) = nothing


################################################################################
# Pair-coupling contraction
################################################################################

# UnitPart data (Δpos and Δcell) for a bare atom.
function unit_offsets_for_atom(layout::UnitLayout, atom)
    (; unit, part) = layout.atom_to_unit[atom]
    return layout.unit_to_parts[unit][part]
end

# Convert a bond between physical atoms into the corresponding bond between
# contracted units. If `bond.n` points from atom i in one cell to atom j in a
# neighboring cell, then the unit-to-unit cell displacement must be corrected by
# the two atoms' cell offsets inside their units.
function contracted_bond(bond::Bond, layout::UnitLayout)
    ui = layout.atom_to_unit[bond.i]
    uj = layout.atom_to_unit[bond.j]
    n = SVector{3,Int}(bond.n)
    Δcell_i = unit_offsets_for_atom(layout, bond.i).Δcell
    Δcell_j = unit_offsets_for_atom(layout, bond.j).Δcell
    n_new = n + Δcell_i - Δcell_j
    return Bond(ui.unit, uj.unit, n_new)
end

# Checks whether a bond given in terms of the original crystal lies inside a
# single unit of the contracted system.
function bond_is_in_unit(bond::Bond, layout::UnitLayout)
    newbond = contracted_bond(bond, layout)
    return newbond.i == newbond.j && all(iszero, newbond.n)
end

# Converts what was a pair coupling between two different sites of a single unit
# in the original system into an "onsite coupling" for the entangled unit.
# Accumulates to `op` in place.
function accum_intra_unit_pair_coupling!(op, pc, uncontracted_Ns, layout::UnitLayout)
    pc.isculled && return

    ui = layout.atom_to_unit[pc.bond.i]
    uj = layout.atom_to_unit[pc.bond.j]
    Ns_unit = uncontracted_Ns[atoms_in_unit(layout, ui.unit)]
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
function promote_inter_unit_pair_coupling(pc, uncontracted_Ns, layout::UnitLayout)
    (; i, j) = pc.bond
    ui = layout.atom_to_unit[i]
    uj = layout.atom_to_unit[j]

    Ns1 = uncontracted_Ns[atoms_in_unit(layout, ui.unit)]
    Ns2 = uncontracted_Ns[atoms_in_unit(layout, uj.unit)]
    Ni = uncontracted_Ns[1, 1, 1, i]
    Nj = uncontracted_Ns[1, 1, 1, j]

    # Decompose the pair operator as ∑ₖ aₖ ⊗ bₖ. Then promote each factor to the
    # tensor product space, aₖ → Aₖ.
    data = map(pair_coupling_tensor_data(pc, Ni, Nj)) do (a, b)
        A = Hermitian(local_op_to_product_space(Matrix(a), ui.part, Ns1))
        B = Hermitian(local_op_to_product_space(Matrix(b), uj.part, Ns2))
        (A, B)
    end

    newbond = contracted_bond(pc.bond, layout)
    gen1 = spin_matrices_of_dim(; N=prod(Ns1))
    gen2 = spin_matrices_of_dim(; N=prod(Ns2))
    return (; newbond, pc.scalar, tensordec=TensorDecomposition(gen1, gen2, data))
end

# Map ModelParam of an uncontracted system to that of an entangled system.
function contract_param(param::ModelParam, uncontracted_Ns, layout::UnitLayout)
    nunits = length(layout.unit_to_parts)

    # Intra-unit terms become onsite (on-bond) operators, one per touched unit.
    onsites = Tuple{Int, OnsiteCoupling}[]
    for u in 1:nunits
        Ns = uncontracted_Ns[atoms_in_unit(layout, u)]
        unit_operator = zeros(ComplexF64, prod(Ns), prod(Ns))
        relevant = atoms_in_unit(layout, u)

        # Bare onsite couplings embedded into the unit's product space.
        for (atom, oc) in param.onsites
            atom in relevant || continue
            (; part) = layout.atom_to_unit[atom]
            unit_operator += local_op_to_product_space(oc, part, Ns)
        end

        # Intra-unit pair couplings fold into the same onsite operator.
        for pc in param.pairs
            if bond_is_in_unit(pc.bond, layout) && layout.atom_to_unit[pc.bond.i].unit == u
                accum_intra_unit_pair_coupling!(unit_operator, pc, uncontracted_Ns, layout)
            end
        end

        iszero(unit_operator) || push!(onsites, (u, Hermitian(unit_operator)))
    end

    # Inter-unit pair couplings become pair couplings between contracted units.
    pairs = PairCoupling[]
    for pc in param.pairs
        bond_is_in_unit(pc.bond, layout) && continue
        (; newbond, scalar, tensordec) = promote_inter_unit_pair_coupling(pc, uncontracted_Ns, layout)
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
    unit_to_parts = Vector{UnitPart}[]
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
        parts = map(enumerate(zip(atoms, Δcells, parts_pos))) do (part, (atom, Δcell, pos))
            atom_to_unit[atom] = UnitAndPartIndex(unit_idx, part)
            UnitPart(atom, pos - unit_pos, Δcell - unit_cell)
        end
        push!(unit_to_parts, parts)
    end

    # Reuse the spacegroup data, but drop the symops to mark unavailability of
    # symmetry analysis (in analogy to `reshape_crystal`).
    new_types = fill("", length(new_positions))
    new_classes = 1:length(unit_atoms)
    new_sg = Spacegroup(SymOp[], crystal.sg.label, crystal.sg.number, crystal.sg.setting)
    new_crystal = Crystal(nothing, crystal.latvecs, crystal.recipvecs, new_positions,
                          new_types, new_classes, new_sg)

    layout = UnitLayout(atom_to_unit, unit_to_parts)

    return new_crystal, layout
end

# Populate the entanglement metadata of a (possibly reshaped) contracted `sys`,
# given the corresponding `uncontracted` system. The `orig_layout` defines the
# entanged units in the context of the "original" chemical crystal.
function rebuild_entanglement!(sys::System, uncontracted, orig_layout::UnitLayout)
    nunits = natoms(sys.crystal)
    nparts = length(first(orig_layout.unit_to_parts))
    natoms_bare = natoms(uncontracted.crystal)
    @assert natoms_bare == nunits * nparts

    # Build a UnitLayout for the possibly reshaped sys.crystal
    atom_to_unit = UnitAndPartIndex[]
    unit_to_parts = [Vector{UnitPart}(undef, nparts) for _ in 1:nunits]
    for a in 1:natoms_bare
        orig_a = map_atom_to_other_crystal(uncontracted.crystal, a, orig_crystal(uncontracted))
        p = orig_layout.atom_to_unit[orig_a].part

        # Fractional displacement from the unit center to the part's atom
        orig_Δpos = unit_offsets_for_atom(orig_layout, orig_a).Δpos
        global_Δpos = orig_crystal(sys).latvecs * orig_Δpos
        Δpos = uncontracted.crystal.latvecs \ global_Δpos

        # Cell offset from the unit center to the part's atom at cell [0, 0, 0]
        unit_center = global_position_at(uncontracted, (1, 1, 1, a)) - global_Δpos
        u, unit_cell = position_to_atom_and_offset(sys.crystal, sys.crystal.latvecs \ unit_center)
        Δcell =  Vec3(0, 0, 0) - unit_cell

        # Store offset data
        unit_to_parts[u][p] = UnitPart(a, Δpos, Δcell)

        # Store unit part data
        push!(atom_to_unit, UnitAndPartIndex(u, p))
    end
    layout = UnitLayout(atom_to_unit, unit_to_parts)

    lifted_spin_ops = map(1:natoms_bare) do atom
        S_local = spin_matrices_of_dim(; N=uncontracted.Ns[1, 1, 1, atom])
        (; unit, part) = layout.atom_to_unit[atom]
        Ns_unit = uncontracted.Ns[atoms_in_unit(layout, unit)]
        ntuple(α -> local_op_to_product_space(S_local[α], part, Ns_unit), 3)
    end

    sys.entanglement = Entanglement(uncontracted, layout, lifted_spin_ops)
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
the cell boundary, however, then all atoms must be paired with cell offsets. For
example, replacing `[1, 2]` with `[(1, [0, 0, 0]), (2, [-1, 0, 0])]` would
indicate that atom `2` participates in the dimer via its periodic image in the
neighboring cell ``-𝐚_1``.

The input system should be in `:SUN` mode with interactions already populated.
These interactions will be transferred to the entangled system in tensor-product
form.
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
    contracted_crystal, layout = contract_crystal(uncontracted.crystal, unit_atoms, unit_Δcells)

    # Local Hilbert space dimension per unit. TODO: Relax uniformity constraint.
    nunits = natoms(contracted_crystal)
    Ns_contracted = [prod(uncontracted.Ns[atoms_in_unit(layout, u)]) for u in 1:nunits]
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
    sys.params = [contract_param(p, uncontracted.Ns, layout) for p in uncontracted.params]
    repopulate_couplings_from_params!(sys)

    # Set entanglement metadata in sys. Keep a clone of the uncontracted system,
    # e.g., as a helper for dipole-level physics.
    uncontracted = clone_system(uncontracted)
    rebuild_entanglement!(sys, uncontracted, layout)

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

    (; uncontracted, layout) = get_entanglement(sys)

    nobs = num_observables(measure)
    dims = sys.dims
    nunits = length(layout.unit_to_parts)
    nparts = length(layout.unit_to_parts[1])  # uniform by construction

    Op = eltype(measure.observables)
    new_obs     = Array{Op, 6}(undef, nobs, dims..., nunits, nparts)
    new_offsets = zeros(Vec3, nunits, nparts)
    new_ff      = Array{FormFactor, 3}(undef, nobs, nunits, nparts)

    for u in 1:nunits
        Ns_unit = uncontracted.Ns[atoms_in_unit(layout, u)]
        for (p, (; atom, Δpos, Δcell)) in enumerate(layout.unit_to_parts[u])
            new_offsets[u, p] = Δpos
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
