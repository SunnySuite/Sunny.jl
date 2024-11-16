# Basis vectors for primitive reciprocal lattice.
function prim_recipvecs(cryst::Crystal)
    prim_latvecs = cryst.latvecs * primitive_cell(cryst)
    return 2π * inv(prim_latvecs)'
end

function brillouin_zone_path(cryst::Crystal)
    # Lattice vectors in ITA standard setting, [aₛ bₛ cₛ] = [a, b, c] P⁻¹.
    std_latvecs = cryst.latvecs / cryst.sg.setting.R

    # High symmetry path in Brillouin zone
    bzpath = Brillouin.irrfbz_path(cryst.sg.number, eachcol(std_latvecs))

    # Verify primitive basis by calculating two ways
    @assert Brillouin.basis(bzpath) ≈ eachcol(prim_recipvecs(cryst))

    # Map from "primitive cell" coordinates to conventional RLU
    A = inv(primitive_cell(cryst)')
    @assert A ≈ cryst.recipvecs \ prim_recipvecs(cryst)
    map!(q -> A * q, values(bzpath.points))
    return (; bzpath.paths, bzpath.points)
end

"""
    special_points(cryst::Crystal)

The high-symmetry points of the irreducible Brillouin zone for `cryst`. Returns
a dictionary that maps symbol names to positions in conventional reciprocal
lattice units (RLU). Data is obtained from the
[Brillouin.jl](https://github.com/thchr/Brillouin.jl) package and originates
from SeeK-path [1, 2].

## Examples

```jl
latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], 227)
pts = special_points(cryst)
@assert keys(pts) == Set([:U, :W, :K, :Γ, :L, :X])
@assert pts[:Γ] == [0, 0, 0]
@assert pts[:U] == [1/4, 1, 1/4]
```

## References

1. [Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, _Band structure diagram
   paths based on crystallography_, Comp. Mat. Sci. **128**, 140
   (2017)](https://doi.org/10.1016/j.commatsci.2016.10.015).
2. SeeK-path software, [available
   online](https://www.materialscloud.org/work/tools/seekpath).
"""
function special_points(cryst::Crystal)
    return brillouin_zone_path(cryst).points
end

"""
    special_paths(cryst::Crystal)

Standardized paths that connect high-symmetry points within the irreducible
Brillouin zone for `cryst`. Data is obtained from the
[Brillouin.jl](https://github.com/thchr/Brillouin.jl) package and originates
from SeeK-path [1, 2].

## Examples

```jl
latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
cryst = Crystal(latvecs, [[0, 0, 0]], 227)
paths = special_paths(cryst)
@assert paths == [[:Γ, :X, :U], [:K, :Γ, :L, :W, :X]]
```

See also [`special_points`](@ref) which maps from symbols to coordinates in
reciprocal lattice units (RLU).

## References

1. [Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, _Band structure diagram
   paths based on crystallography_, Comp. Mat. Sci. **128**, 140
   (2017)](https://doi.org/10.1016/j.commatsci.2016.10.015).
2. SeeK-path software, [available
   online](https://www.materialscloud.org/work/tools/seekpath).
"""
function special_paths(cryst::Crystal)
    return brillouin_zone_path(cryst).paths
end
