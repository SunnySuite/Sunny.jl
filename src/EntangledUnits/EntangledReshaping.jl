#TODO: Test this rigorously -- this is key to reshaping EntangledSystems and
# making EntangledSpinWaveTheorys.

# An entangled system can be specified with a list of
# tuples, e.g. [(1,2), (3,4)], which will group the original sites 1 and 2 into
# one unit and 3 and 4 into another. Given an EntangledSystem constructed from a
# sys_uncontracted and a reshaped sys_uncontracted, this function returns a list of
# tuples for specifying the corresponding entanglement of the reshaped system.
function units_for_reshaped_system(reshaped_sys_uncontracted, esys)
    (; sys_uncontracted) = esys
    units = original_unit_spec(esys)
    new_crystal = reshaped_sys_uncontracted.crystal
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
        new_position = position_at(reshaped_sys_uncontracted, new_site)

        # Find corresponding original atom number.
        site = position_to_site(sys_uncontracted, position_at(reshaped_sys_uncontracted, new_site))
        original_atom = site[4]
        position_of_corresponding_atom = position_at(sys_uncontracted, (1, 1, 1, original_atom))
        offset = new_position - position_of_corresponding_atom

        # Find the unit to which this original atom belongs.
        unit = units[findfirst(unit -> original_atom in unit, units)]

        # Find positions of all atoms in the unit, map into the reshaped
        # crystal, and reject any mapping that requires a nonzero lattice
        # offset (indicating a split across crystallographic cells).
        unit_positions = [position_at(sys_uncontracted, (1, 1, 1, atom)) for atom in unit]
        new_unit_sites = map(unit_positions) do position
            position_to_atom_and_offset(
                reshaped_sys_uncontracted.crystal,
                reshaped_sys_uncontracted.crystal.latvecs \ orig_crystal(reshaped_sys_uncontracted).latvecs * (position + offset),
            )
        end
        new_unit = Int64[]
        for (a, n) in new_unit_sites
            if !all(==(0), n)
                error("Given shape incompatible with entangled unit structure. Unit split between crystallographic cells.")
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

function reshape_supercell(esys::EntangledSystem, shape)
    (; sys_contracted, sys_uncontracted) = esys

    # Reshape the origin System.
    sys_uncontracted_new = reshape_supercell(sys_uncontracted, shape)

    # Reshape the the underlying "entangled" System.
    units_new = units_for_reshaped_system(sys_uncontracted_new, esys)
    contracted_crystal, contraction_info = contract_crystal(sys_uncontracted_new.crystal, units_new)
    sys_contracted_new = reshape_supercell(sys_contracted, shape)

    # The two reshape paths must agree on crystal positions
    contracted_crystal.positions ≈ sys_contracted_new.crystal.positions || error("Given shape incompatible with entangled unit structure.")

    # Construct dipole operator field for reshaped EntangledSystem
    dipole_operators_origin = all_dipole_observables(sys_uncontracted_new; apply_g=false)
    (; observables, source_idcs) = observables_to_product_space(dipole_operators_origin, sys_uncontracted_new, contraction_info)

    return EntangledSystem(sys_contracted_new, sys_uncontracted_new, contraction_info, observables, source_idcs)
end

function repeat_periodically(esys::EntangledSystem, counts)
    (; sys_contracted) = esys
    all(>=(1), counts) || error("Require at least one count in each direction.")
    shape = cell_shape(sys_contracted) * diagm(Vec3(sys_contracted.dims .* counts))
    return reshape_supercell(esys, shape)
end

function resize_supercell(esys::EntangledSystem, dims::NTuple{3,Int})
    return reshape_supercell(esys, diagm(Vec3(dims)))
end
