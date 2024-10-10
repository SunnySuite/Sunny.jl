#TODO: Test this rigorously -- this is key to reshaping EntangledSystems and
# making EntangledSpinWaveTheorys. 

# An entangled system can be specified with a list of
# tuples, e.g. [(1,2), (3,4)], which will group the original sites 1 and 2 into
# one unit and 3 and 4 into another. Given an EntangledSystem constructed from a
# sys_origin and and a reshaped sys_origin, this function returns a list of
# tuples for specifying the corresponding entanglement of the reshaped system.
function units_for_reshaped_system(reshaped_sys_origin, esys)
    (; sys_origin) = esys                 
    units = original_unit_spec(esys)
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
        new_site = CartesianIndex(1, 1, 1, new_atom) # Only need to define for a single unit cell, may as well be the first 
        new_position = position(reshaped_sys_origin, new_site)

        # Find corresponding original atom number.
        site = position_to_site(sys_origin, position(reshaped_sys_origin, new_site))
        original_atom = site[4]
        position_of_corresponding_atom = position(sys_origin, (1, 1, 1, original_atom))
        offset = new_position - position_of_corresponding_atom

        # Find the unit to which this original atom belongs.
        unit = units[findfirst(unit -> original_atom in unit, units)]

        # Find positions of all atoms in the unit, find corresponding sites in reshape system, and define unit for reshaped system.
        unit_positions = [position(sys_origin, CartesianIndex(1, 1, 1, atom)) for atom in unit]
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

function reshape_supercell(esys::EntangledSystem, shape)
    (; sys, sys_origin) = esys

    # Reshape the origin System.
    sys_origin_new = reshape_supercell(sys_origin, shape)

    # Reshape the the underlying "entangled" System.
    units_new = units_for_reshaped_system(sys_origin_new, esys)
    _, contraction_info = contract_crystal(sys_origin_new.crystal, units_new)
    sys_new = reshape_supercell(sys, shape)

    # Construct dipole operator field for reshaped EntangledSystem
    dipole_operators_origin = all_dipole_observables(sys_origin_new; apply_g=false) 
    (; observables, source_idcs) = observables_to_product_space(dipole_operators_origin, sys_origin_new, contraction_info)

    return EntangledSystem(sys_new, sys_origin_new, contraction_info, observables, source_idcs)
end

function repeat_periodically(esys, counts)
    (; sys, sys_origin, contraction_info) = esys

    # Repeat both entangled and original system periodically
    sys_new = repeat_periodically(sys, counts)
    sys_origin_new = repeat_periodically(sys_origin, counts)

    # Construct dipole operator field for reshaped EntangledSystem
    dipole_operators_origin = all_dipole_observables(sys_origin_new; apply_g=false) 
    (; observables, source_idcs) = observables_to_product_space(dipole_operators_origin, sys_origin_new, contraction_info)

    return EntangledSystem(sys_new, sys_origin_new, contraction_info, observables, source_idcs)
end

function resize_supercell(esys::EntangledSystem, dims::NTuple{3,Int})
    (; sys, sys_origin, contraction_info) = esys

    # Resize both entangled and original system periodically
    sys_new = reshape_supercell(sys, diagm(Vec3(dims)))
    sys_origin_new = reshape_supercell(sys_origin, diagm(Vec3(dims)))

    # Construct dipole operator field for reshaped EntangledSystem
    dipole_operators_origin = all_dipole_observables(sys_origin_new; apply_g=false) 
    (; observables, source_idcs) = observables_to_product_space(dipole_operators_origin, sys_origin_new, contraction_info)

    return EntangledSystem(sys_new, sys_origin_new, contraction_info, observables, source_idcs)
end