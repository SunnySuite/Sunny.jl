"""
    SpinInfo(atom::Int; S, g=2)

Characterizes the spin at a given `atom` index within the crystal unit cell. `S`
is an integer multiple of 1/2 and gives the spin angular momentum in units of ħ.
`g` is the g-factor or tensor, such that an angular momentum dipole ``s``
produces a magnetic moment ``g s`` in units of the Bohr magneton.
"""
struct SpinInfo
    atom   :: Int     # Index of atom in unit cell
    S      :: Float64 # Spin magnitude in units of ħ
    g      :: Mat3    # Spin g-tensor

    function SpinInfo(atom::Int; S, g)
        S > 0 || error("Spin S must be positive. Use `subcrystal` to discard non-magnetic ions.")
        isinteger(2S) || error("Spin S must be an exact multiple of 1/2")
        g = typeof(g) <: Number ? Mat3(I*g) : Mat3(g)
        new(atom, S, g)
    end
end

# Propagates spin magnitudes and symmetry-transformed g-tensors to all
# symmetry-equivalent atoms. Throws an error if two symmetry-equivalent atoms
# are provided in `infos`, or if some atoms remain unspecified.
function propagate_site_info(cryst::Crystal, infos::Vector{SpinInfo})
    # Verify that all g tensors are consistent with the the site symmetries
    for info in infos
        if !is_coupling_valid(cryst, Bond(info.atom, info.atom, (0,0,0)), info.g)
            error("g-tensor $(info.g) is inconsistent with the site symmetry of atom $(info.atom).")
        end
    end

    ref_atoms = [info.atom for info in infos]
    atom_to_ref_atom = propagate_reference_atoms(cryst, ref_atoms)

    return map(enumerate(atom_to_ref_atom)) do (a, a′)
        info = infos[findfirst(==(a′), ref_atoms)]
        S = info.S
        g = transform_coupling_for_bonds(cryst, Bond(a,a,(0,0,0)), Bond(a′,a′,(0,0,0)), info.g)
        SpinInfo(a; S, g)
    end
end
