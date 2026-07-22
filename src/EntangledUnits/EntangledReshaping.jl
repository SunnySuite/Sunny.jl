#TODO: Test this rigorously -- this is key to reshaping entangled systems and
# making spin wave theories from them.

# An entangled system can be specified with a list of
# tuples, e.g. [(1,2), (3,4)], which will group the original sites 1 and 2 into
# one unit and 3 and 4 into another. Given an entangled `System` and a reshaped
# copy of its physical (bare) system, this function returns a list of tuples for
# specifying the corresponding entanglement of the reshaped system.
function units_for_reshaped_system(reshaped_sys_origin, sys)
    sys_origin = get_entanglement(sys).bare_system
    units = original_unit_spec(sys)
    new_crystal = reshaped_sys_origin.crystal
    new_atoms = collect(1:natoms(new_crystal))
    new_units = []

    # Take a list of all the new atoms in the reshaped system. Pick the first.
    # Map it back to the original system to determine what unit it belongs to.
    # Then map all members of the unit forward to define the unit in terms of
    # the atoms of the reshaped system. Remove these atoms from the list of
    # sites left to be "entangled" and repeat until list of new atoms is
    # exhausted.
    while length(new_atoms) > 0
        # Pick any site from list of new sites
        new_atom = new_atoms[1]
        new_site = (1, 1, 1, new_atom) # Only need to define for a single unit cell, may as well be the first 
        new_position = position_at(reshaped_sys_origin, new_site)

        # Find corresponding original atom number.
        site = position_to_site(sys_origin, position_at(reshaped_sys_origin, new_site))
        original_atom = site[4]
        position_of_corresponding_atom = position_at(sys_origin, (1, 1, 1, original_atom))
        offset = new_position - position_of_corresponding_atom

        # Find the unit to which this original atom belongs.
        unit = units[findfirst(unit -> original_atom in unit, units)]

        # Find positions of all atoms in the unit, find corresponding sites in reshape system, and define unit for reshaped system.
        unit_positions = [position_at(sys_origin, (1, 1, 1, atom)) for atom in unit]
        new_unit_sites = [position_to_site(reshaped_sys_origin, position + offset) for position in unit_positions]
        new_unit = Int64[]
        for new_site in new_unit_sites
            i, j, k, a = new_site.I
            if !(i == j == k == 1)
                error("Specified reshaping incompatible with specified entangled units. (Unit split between crystalographic unit cells.)")
            end
            push!(new_unit, a)
        end
        push!(new_units, Tuple(new_unit))

        # Delete members of newly defined unit from list of all atoms in the
        # reshaped system.
        idcs = findall(atom -> atom in new_unit, new_atoms)
        deleteat!(new_atoms, idcs)
    end

    return new_units
end

# Reshaping of an entangled `System`. The contracted system and the physical
# `bare_system` are reshaped together, and fresh entanglement metadata is
# attached to the reshaped contracted system. These are called from the ordinary
# `reshape_supercell`/`repeat_periodically`/`resize_supercell` (Reshaping.jl)
# when `sys` is entangled.
function reshape_supercell_entangled(sys::System, shape)
    bare_system = get_entanglement(sys).bare_system

    # Reshape the physical (bare) system.
    bare_system_new = reshape_supercell(bare_system, shape)

    # Determine the entangled units of the reshaped system, then reshape the
    # contracted system and attach fresh entanglement metadata.
    units_new = units_for_reshaped_system(bare_system_new, sys)
    _, contraction_info = contract_crystal(bare_system_new.crystal, units_new)
    sys_new = reshape_contracted(sys, s -> reshape_supercell(s, shape))

    return attach_entanglement!(sys_new, bare_system_new, contraction_info)
end

function repeat_periodically_entangled(sys::System, counts)
    (; bare_system, contraction_info) = get_entanglement(sys)

    # Repeat both the contracted and physical systems periodically. The unit
    # structure is unchanged, so reuse the contraction info (copied).
    sys_new = reshape_contracted(sys, s -> repeat_periodically(s, counts))
    bare_system_new = repeat_periodically(bare_system, counts)

    return attach_entanglement!(sys_new, bare_system_new, copy(contraction_info))
end

function resize_supercell_entangled(sys::System, dims::NTuple{3,Int})
    return reshape_supercell_entangled(sys, diagm(Vec3(dims)))
end

# Flatten an entangled `sys` into a single `(1,1,1)` supercell, preserving
# entanglement metadata. This mirrors the plain flatten performed in the
# `SpinWaveTheory` constructor (`reshape_supercell_aux` onto
# `resize_and_flatten_crystal`), but reshapes the physical `bare_system` in
# tandem and reattaches fresh metadata so that the flattened contracted system
# still carries `moment_operators` (needed for the Zeeman term in `swt_data`).
function flatten_supercell_entangled(sys::System)
    bare = get_entanglement(sys).bare_system

    # Flatten the contracted system via the ordinary (non-entangled) path.
    new_cryst = resize_and_flatten_crystal(sys.crystal, sys.dims)
    sys_flat = reshape_contracted(sys, s -> reshape_supercell_aux(s, new_cryst, (1, 1, 1)))

    # Flatten the physical system the same way.
    bare_cryst = resize_and_flatten_crystal(bare.crystal, bare.dims)
    bare_flat = reshape_supercell_aux(bare, bare_cryst, (1, 1, 1))

    # Recompute contraction info for the flattened system and reattach metadata.
    units_new = units_for_reshaped_system(bare_flat, sys)
    _, contraction_info = contract_crystal(bare_flat.crystal, units_new)
    return attach_entanglement!(sys_flat, bare_flat, contraction_info)
end

# Apply an ordinary reshaping operation `f` to the *contracted* part of an
# entangled `sys`, treating it as a plain system. The `entanglement` field is
# temporarily detached so `f` (which dispatches on `is_entangled`) takes the
# ordinary path and does not recurse. The reshaped system returned by `f`
# already carries `entanglement = nothing` (built via the positional `System`
# constructor); the caller re-attaches fresh metadata.
function reshape_contracted(sys::System, f)
    ent = sys.entanglement
    sys.entanglement = nothing
    try
        return f(sys)
    finally
        sys.entanglement = ent
    end
end