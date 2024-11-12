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
        s > 0 || error("Spin S must be positive. Use `subcrystal` to discard non-magnetic ions.")
        isinteger(2s) || error("Spin S must be an exact multiple of 1/2")
        g = typeof(g) <: Number ? Mat3(I*g) : Mat3(g)
        new(s, g)
    end
end


# Propagates each atom-moment pair to every symmetry-equivalent atom in the
# crystal. Throws an error if two symmetry-equivalent atoms are provided in
# `moments`, or if some atoms remain unspecified.
function propagate_moments(cryst::Crystal, moments::Vector{Pair{Int, Moment}})
    # Verify that all g tensors are consistent with the the site symmetries
    for (i, m) in moments
        1 <= i <= natoms(cryst) || error("Atom $i outside the valid range 1:$(natoms(cryst))")
        if !is_coupling_valid(cryst, Bond(i, i, [0,0,0]), m.g)
            error("g-tensor on site $i is symmetry inconsistent; see `print_site(cryst, $i)` for more information.")
        end
    end

    ref_atoms = [i for (i, _) in moments]
    ref_classes = cryst.classes[ref_atoms]

    return map(enumerate(cryst.classes)) do (i, c)
        js = findall(==(c), ref_classes)
        isempty(js) && error("Not all sites are specified; consider including atom $i.")
        length(js) > 1 && error("Atoms $(ref_atoms[js]) are symmetry equivalent.")
        (j, m) = moments[only(js)]
        g = transform_coupling_for_bonds(cryst, Bond(i, i, [0,0,0]), Bond(j, j, [0,0,0]), m.g)
        Moment(; m.s, g)
    end
end
