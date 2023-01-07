"""
    SiteInfo(atom::Int; S, g=2)

Characterizes the spin at a given `atom` index within the crystal unit cell. `S`
is an integer multiple of 1/2 and gives the spin angular momentum in units of ħ.
`g` is the g-factor or tensor, such that an angular momentum dipole ``s``
produces a magnetic moment ``g s``.
"""
struct SiteInfo
    atom   :: Int     # Index of atom in unit cell
    S      :: Float64 # Spin magnitude in units of ħ
    g      :: Mat3    # Spin g-tensor

    function SiteInfo(atom::Int; S, g=2)
        if !isinteger(2S)
            error("Spin S=$S for atom $atom is not a multiple of 1/2")
        end
        # If g is scalar then convert to tensor
        g = typeof(g) <: Number ? Mat3(I*g) : Mat3(g)
        new(atom, Float64(S), g)
    end
end


# For each atom in the unit cell of cryst, return the corresponding reference
# atom in ref_atom that is symmetry equivalent.
function propagate_reference_atoms(cryst::Crystal, ref_atoms::Vector{Int})
    # Sort infos by site equivalence class
    ref_atoms = sort(ref_atoms; by = (a -> cryst.classes[a]))
    ref_classes = cryst.classes[ref_atoms]

    # Verify that none of the atoms belong to the same class
    for i = 1:length(ref_atoms)-1
        a1, a2 = ref_atoms[[i,i+1]]
        c1, c2 = ref_classes[[i,i+1]]
        if c1 == c2
            error("Atoms $a1 and $a2 are symmetry equivalent.")
        end
    end
    @assert allunique(ref_classes)

    # Verify that every class has been specified
    missing_classes = setdiff(cryst.classes, ref_classes)
    if !isempty(missing_classes)
        c = first(missing_classes)
        a = findfirst(==(c), cryst.classes)
        error("Not all sites are specified; consider including atom $a.")
    end
    @assert length(ref_atoms) == length(unique(cryst.classes))

    # Return a symmetry-equivalent reference atom for each atom in the unit cell
    return map(1:nbasis(cryst)) do a
        c = cryst.classes[a]
        ref_atoms[only(findall(==(c), ref_classes))]
    end
end

"""
    propagate_site_info(cryst::Crystal, site_infos::Vector{SiteInfo})

Propagates spin magnitudes and symmetry-transformed g-tensors to all
symmetry-equivalent atoms. Throws an error if two symmetry-equivalent atoms are
provided in `site_infos`, or if some atoms remain unspecified.
"""
function propagate_site_info(cryst::Crystal, infos::Vector{SiteInfo})
    # Verify that all g tensors are consistent with the the site symmetries
    for info in infos
        if !is_coupling_valid(cryst, Bond(info.atom, info.atom, (0,0,0)), info.g)
            error("g-tensor $(info.g) is inconsistent with the site symmetry of atom $(info.atom).")
        end
    end

    ref_atoms = [info.atom for info in infos]
    atoms = propagate_reference_atoms(cryst, ref_atoms)

    return map(enumerate(infos[atoms])) do (a, info)
        a′ = info.atom
        g′ = info.g
        g = transform_coupling_for_bonds(cryst, Bond(a,a,(0,0,0)), Bond(a′,a′,(0,0,0)), g′)
        SiteInfo(a; info.S, g)
    end
end
