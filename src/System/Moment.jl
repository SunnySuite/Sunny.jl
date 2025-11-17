"""
    Moment(; s, g)

Characterizes a effective spin magnetic moment on an atom. Quantum spin-`s` is a
multiple of 1/2 in units of ħ. The `g`-factor or tensor defines the
[`magnetic_moment`](@ref) ``μ = - g 𝐒`` in units of the Bohr magneton.

# Example
```julia
Moment(s=3/2, g=2)
```
"""
struct Moment
    s :: Float64 # quantum spin
    g :: Mat3    # g-tensor

    function Moment(; s, g)
        s > 0 || error("Spin s must be positive. Use `subcrystal` to discard non-magnetic ions.")
        isinteger(2s) || error("Spin s must be an exact multiple of 1/2")
        g = typeof(g) <: Number ? Mat3(I*g) : Mat3(g)
        new(s, g)
    end
end


# Given `pairs` that map reference atoms to values in `cryst`, return a full
# vector of values for all atoms in `new_cryst`.
function propagate_atom_data(cryst, new_cryst, pairs::Vector{Pair{Int, T}}; check_consistent=(i, v)->true, transform=(i, j, v)->v) where T
    for (i, v) in pairs
        1 <= i <= natoms(cryst) || error("Atom $i outside the valid range 1:$(natoms(cryst))")
        check_consistent(i, v)
    end

    # Reference data with respect to original crystal
    ref_atoms = [i for (i, _) in pairs]
    ref_classes = cryst.classes[ref_atoms]

    # One value for each atom in the original crystal
    data = map(enumerate(cryst.classes)) do (i, c)
        idxs = findall(==(c), ref_classes)
        isempty(idxs) && error("Not all sites are specified; consider including atom $i.")
        length(idxs) > 1 && error("Atoms $(ref_atoms[idxs]) are symmetry equivalent.")
        (j, v) = pairs[only(idxs)]
        transform(i, j, v)
    end

    # Map this data to the new crystal if distinct
    if new_cryst == cryst
        return data
    else
        return map(new_cryst.positions) do new_r
            r = cryst.latvecs \ new_cryst.latvecs * new_r
            data[position_to_atom(cryst, r)]
        end
    end
end


# Propagates each atom-moment pair to every symmetry-equivalent atom in the
# crystal. Throws an error if two symmetry-equivalent atoms are provided in
# `moments`, or if some atoms remain unspecified.
function propagate_moments(cryst::Crystal, moments::Vector{Pair{Int, Moment}})
    function check_consistent(i, v)
        if !is_coupling_valid(cryst, Bond(i, i, [0,0,0]), v.g)
            error("g-tensor on site $i is symmetry inconsistent; see `print_site(cryst, $i)` for more information.")
        end
    end
    function transform(i, j, v)
        g = transform_coupling_for_bonds(cryst, Bond(i, i, [0,0,0]), Bond(j, j, [0,0,0]), v.g)
        Moment(; v.s, g)
    end
    return propagate_atom_data(cryst, cryst, moments; check_consistent, transform)
end
