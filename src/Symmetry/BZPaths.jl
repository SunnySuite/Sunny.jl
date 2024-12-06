# Basis vectors for primitive reciprocal lattice.
function prim_recipvecs(cryst::Crystal)
    prim_latvecs = cryst.latvecs * primitive_cell(cryst)
    return 2π * inv(prim_latvecs)'
end

# Equivalent to Bravais.bravaistype(sgnum, 3, normalize=false) with the Bravais
# package.
function standard_bravais_type(sgnum)
    hall = standard_setting[sgnum]
    letters = Dict(triclinic => 'a', monoclinic => 'm', orthorhombic => 'o',
                   tetragonal => 't', hexagonal => 'h', cubic => 'c')
    return letters[cell_type(hall)] * centering_symbol(hall)
end

function extended_bravais_symbol(cryst)
    sgnum = cryst.sg.number
    std_latvecs = cryst.latvecs / cryst.sg.setting.R
    bt = standard_bravais_type(sgnum)
    return Brillouin.KPaths.extended_bravais(sgnum, bt, SVector{3}(eachcol(std_latvecs)), Val{3}())
end


"""
    print_irreducible_bz_paths(cryst::Crystal)

Prints certain high-symmetry points with suggested paths between them. The
points lie in the irreducible Brillouin, which is reduced from the first
Brillouin zone by the point group symmetries of the crystal. Coordinates are
printed in reciprocal lattice units (RLU), i.e., as multiples of the
conventional lattice vectors of `cryst`. These high-symmetry paths were
originally formulated in Ref. [1] and implemented in
[SeeK-path](https://github.com/giovannipizzi/seekpath). Sunny obtains this
functionality from the [Brillouin.jl](https://github.com/thchr/Brillouin.jl)
package.

See also [`view_bz`](@ref) for an interactive visualization of these
high-symmetry paths.

## References

1. [Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, _Band structure diagram
   paths based on crystallography_, Comp. Mat. Sci. **128**, 140
   (2017)](https://doi.org/10.1016/j.commatsci.2016.10.015).
"""
function print_irreducible_bz_paths(cryst::Crystal)
    (; points, paths) = try 
        irreducible_bz_paths(cryst)
    catch e
        if startswith(e.msg, "Triclinic")
            rethrow(ErrorException("""Triclinic lattice angles must currently be all-acute or all-obtuse.
                                          See https://github.com/thchr/Brillouin.jl/issues/34
                                   """))
        else
            rethrow()
        end
    end

    pt_strs = map(collect(keys(points))) do k
        "$k = " * fractional_vec3_to_string(points[k])
    end
    path_strs = [join(String.(path), ", ") for path in paths]

    println("High-symmetry points in irreducible BZ:")
    println("    ", join(pt_strs, "\n    "))
    println("High-symmetry paths:")
    println("    [" * join(path_strs, "]\n    [") * "]")
    println("Extended Bravais symbol for SeeK-path: '$(extended_bravais_symbol(cryst))'")
end

function irreducible_bz_paths(cryst::Crystal)
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

    return (; bzpath.points, bzpath.paths)
end
