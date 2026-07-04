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
    # the atoms of the reshaped system. In so doing ensure that the unit has not
    # been split across multiple reshaped cells. Remove these atoms from the
    # list of sites left to be "entangled" and repeat until list of new atoms is
    # exhausted.
    while length(new_atoms) > 0
        # Pick any site from list of new sites
        new_atom = new_atoms[1]
        new_site = (1, 1, 1, new_atom) # Any representative cell works.
        new_site_global = global_position_at(reshaped_sys_uncontracted, new_site)

        # Map the reshaped site back to an atom in the original crystal.
        new_site_in_original_basis = sys_uncontracted.crystal.latvecs \ new_site_global
        original_site = position_to_site(sys_uncontracted, new_site_in_original_basis)
        original_atom = original_site[4]

        # Find the original unit containing this atom.
        unit = units[findfirst(unit -> original_atom in unit, units)]

        # Build displacements from the reference atom to every atom in the
        # unit, then apply those displacements to the reshaped reference site.
        original_reference_position = position_at(sys_uncontracted, (1, 1, 1, original_atom))
        unit_displacements = [position_at(sys_uncontracted, (1, 1, 1, atom)) - original_reference_position for atom in unit]
        unit_displacements_global = [sys_uncontracted.crystal.latvecs * displacement for displacement in unit_displacements]

        # Map each target position back to a reshaped-crystal atom + lattice
        # offset. Any nonzero offset means the unit was split across cells.
        new_unit_sites = map(unit_displacements_global) do displacement_global
            target_global = new_site_global + displacement_global
            target_in_reshaped_basis = reshaped_sys_uncontracted.crystal.latvecs \ target_global
            position_to_atom_and_offset(reshaped_sys_uncontracted.crystal, target_in_reshaped_basis)
        end
        new_unit = Int64[]
        for (a, lattice_offset) in new_unit_sites
            if !allequal(lattice_offset)
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
