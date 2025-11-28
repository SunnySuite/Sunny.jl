
"""
An object describing a crystallographic unit cell and its space group symmetry.
Constructors are as follows:

    Crystal(filename; symprec=nothing)

Reads the crystal from a CIF or mCIF at `filename`. Lattice vectors will be in
units of angstrom. In most cases, `symprec` can be omitted. If provided, it will
specify the precision of the dimensionless site position data (commonly between
1e-2 and 1e-5).

In case of an mCIF, returns the standardized chemical cell. Use this crystal to
build an appropriately shaped [`System`](@ref) and then load the magnetic order
with [`set_dipoles_from_mcif!`](@ref).

    Crystal(latvecs, positions; types=nothing, symprec=1e-5)

Constructs a crystal from the complete list of atom positions `positions`, with
coordinates (between 0 and 1) in units of lattice vectors `latvecs`. Spacegroup
symmetry information is automatically inferred using the [Spglib
package](https://github.com/spglib/spglib) [1] at the dimensionless symmetry
precision `symprec`. The optional parameter `types` is a list of strings, one
for each atom, and can be used to break symmetry-equivalence between atoms.

    Crystal(latvecs, positions, spacegroup; types=nothing, choice=nothing, symprec=1e-5)

Builds a crystal using the symmetries of a spacegroup. One representative atom
must be specified for each occupied Wyckoff. The `spacegroup` may be passed as a
number 1..230 or as a string. A spacegroup number will imply the ITA standard
setting [2], consistent with conventions from the [Bilbao crystallographic
server](https://www.cryst.ehu.es). Alternatively, a spacegroup symbol can select
among all possible ITA settings (origin `choice` and orientation) for the
spacegroup. In case of ambiguity, an error message will list the options.


# Examples

```julia
# Read a Crystal from a .cif or .mcif file
Crystal("filename.cif")

# Build a BCC crystal in the conventional cubic unit cell by specifying both
# atoms. The spacegroup 229 is inferred.
latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
positions = [[0, 0, 0], [1/2, 1/2, 1/2]]
Crystal(latvecs, positions)

# Build a CsCl crystal (two simple cubic sublattices). Because of the distinct
# atom types, the spacegroup number 221 is now inferred.
types = ["Na", "Cl"]
cryst = Crystal(latvecs, positions; types)

# Build a diamond cubic crystal from its spacegroup number 227 and a single
# atom position belonging to Wyckoff 8a.
positions = [[1/8, 1/8, 1/8]]
cryst = Crystal(latvecs, positions, 227)

# Build the same crystal for 227 (Fd-3m), but with non-standard origin choice.
positions = [[0, 0, 0]]
cryst = Crystal(latvecs, positions, "Fd-3m"; choice="1")

# Spacegroup 5 has the short name "C2", which admits multiple ITA settings.
Crystal(latvecs, positions, 5)         # Standard ITA setting, "C 1 2 1"
Crystal(latvecs, positions, "C2")      # Error due to ambiguous setting
Crystal(latvecs, positions, "B 1 1 2") # Specific non-standard ITA setting
```

See also [`lattice_vectors`](@ref).

## References

1. [A. Togo, K. Shinohara, I. Tanaka, _Spglib: a software library for crystal
   symmetry search_ (2018)
   [arXiv:1808.01590]](https://arxiv.org/abs/1808.01590).
2. [International Tables of Crystallography, Volume A
   (2016)](https://doi.org/10.1107/97809553602060000114).
"""
struct Crystal
    root      :: Union{Nothing, Crystal} # Root crystal (invariant under `subcrystal` and reshaping)
    latvecs   :: Mat3                    # Lattice vectors as columns
    recipvecs :: Mat3                    # Reciprocal lattice vectors (conventional)
    positions :: Vector{Vec3}            # Positions in fractional coords
    types     :: Vector{String}          # Types
    classes   :: Vector{Int}             # Symmetry-equivalent class indices
    sg        :: Spacegroup              # Spacegroup symmetries and setting
end

# Constructs a crystal from the complete list of atom positions `positions`,
# representing fractions (between 0 and 1) of the lattice vectors `latvecs`.
# All symmetry information is automatically inferred.
function Crystal(latvecs, positions; types::Union{Nothing, Vector{String}}=nothing, symprec=1e-5)
    if length(positions) >= 100
        @warn """Very large crystal cells are not recommended. To model chemical inhomogeneity:
                 1. Create the `Crystal` with an idealized chemical cell.
                 2. Create the `System` with large volume dimensions, `dims`.
                 3. Use `to_inhomogeneous` and related functions to set inhomogeneities."""
    end
    latvecs = convert(Mat3, latvecs)
    positions = [convert(Vec3, p) for p in positions]
    if isnothing(types)
        types = fill("", length(positions))
    end
    return crystal_from_inferred_symmetry(latvecs, positions, types; symprec)
end


# Builds a crystal by applying the symmetry operators for a given spacegroup
# (number or symbol).
function Crystal(latvecs, positions, symbol::Union{Int, String}; types::Union{Nothing, Vector{String}}=nothing,
                 setting=nothing, choice::Union{Nothing, String}=nothing, symprec=1e-5)
    if !isnothing(setting)
        if symbol isa Integer
            @warn "`setting` argument is deprecated! Omit to get the ITA standard setting or use `choice`."
        else
            @warn "`setting` argument is deprecated! Use `choice` instead."
        end
        choice = setting
    end

    latvecs = Mat3(latvecs)
    positions = Vec3.(positions)
    if isnothing(types)
        types = fill("", length(positions))
    end

    # Enforce right-handedness convention for aesthetics
    det(latvecs) < 0 && error("Lattice vectors are not right-handed.")

    # Will report an informative error if latvecs not consistent with symbol
    sgt = unique_spacegroup_type(symbol, latvecs; choice)
    sg = Spacegroup(Int(sgt.hall_number))
    return crystal_from_spacegroup(latvecs, positions, types, sg; symprec)
end


# Sunny crystal specification is agnostic to length units. Spglib, however,
# expects lattice vector magnitudes to be order one. The wrappers below do the
# following: (1) Non-dimensionalize the lattice vectors using a natural length
# scale, (2) Perform Spglib symmetry analysis, (3) Re-introduce length
# dimensions in the appropriate return values.
function dimensionless_spg_cell(cell)
    (; lattice, positions, atoms, magmoms) = cell
    # See https://github.com/spglib/spglib/issues/595 for discussion of a0
    a0 = minimum(svdvals(Mat3(lattice)))
    cell = Spglib.Cell(lattice/a0, positions, atoms, magmoms)
    return (; cell, a0)
end
function spg_get_dataset_scaled(cell, symprec)
    (; cell, a0) = dimensionless_spg_cell(cell)
    ret = Spglib.get_dataset(cell, symprec)
    ret.std_lattice .*= a0
    ret.primitive_lattice .*= a0
    return ret
end


"""
    natoms(cryst::Crystal)

Number of atoms in the unit cell, i.e., number of Bravais sublattices.
"""
@inline natoms(cryst::Crystal) = length(cryst.positions)

"""
    cell_volume(cryst::Crystal)

Volume of the crystal unit cell.
"""
cell_volume(cryst::Crystal) = abs(det(cryst.latvecs))


# Human-readable spacegroup.
function spacegroup_label(hall_number::Int)
    sgt = Spglib.get_spacegroup_type(hall_number)
    # Note that sgt.international is more verbose than sgt.international_full
    return "'$(sgt.international)' ($(sgt.number))"
end


# Crystal is immutable in its public API. Therefore, this is an internal
# function should not overload Base.permute!
function permute_sites!(cryst::Crystal, p)
    cryst.positions .= cryst.positions[p]
    cryst.classes .= cryst.classes[p]
    cryst.types .= cryst.types[p]
    return cryst
end

# Sort the sites according to class and fractional coordinates. Any changes here
# would likely break indexing of user scripts. This is an internal function.
function sort_sites!(cryst::Crystal)
    function less_than(i, j)
        ci = cryst.classes[i]
        cj = cryst.classes[j]
        if ci != cj
            return ci < cj
        end

        # Sort in order of (a3, a2, a1).
        ri = cryst.positions[i]
        rj = cryst.positions[j]
        for k = 3:-1:1
            if !is_integer(ri[k]-rj[k]; tol=1e-12)
                return ri[k] < rj[k]
            end
        end

        # This should be impossible after the `validate_positions` or
        # `validate_orbits` check.
        @assert false "Symmetry-equivalent positions $(pos_to_string(ri)) and $(pos_to_string(rj))"
    end
    p = sort(eachindex(cryst.positions), lt=less_than)
    return permute_sites!(cryst, p)
end


"""
    standardize(cryst::Crystal)

Converts the provided crystal to an equivalent one that is standardized
following ITA conventions [1].

1. [International Tables of Crystallography, Volume A
   (2016)](https://doi.org/10.1107/97809553602060000114).
"""
function standardize(cryst::Crystal)
    # Lattice vectors of ITA standard cell
    latvecs = lattice_vectors(lattice_params(cryst.latvecs / cryst.sg.setting.R)...)

    # Map symmetry-distinct atoms to standard cell
    inds = unique(i -> cryst.classes[i], eachindex(cryst.classes))
    positions = transform.(Ref(cryst.sg.setting), cryst.positions[inds])
    types = cryst.types[inds]

    # Construct ITA standard cell for spacegroup number
    ret = Crystal(latvecs, positions, cryst.sg.number; types)

    # Verify consistency of Wyckoffs
    foreach(inds, positions) do i, r
        i′ = position_to_atom(ret, r)
        @assert get_wyckoff(cryst, i).letter == get_wyckoff(ret, i′).letter
    end

    return ret
end


function score_setting(setting::SymOp; latvecs::Mat3)
    # Lattice vectors in proposed setting
    A = latvecs / setting.R

    # Rotation with respect to idealized Cartesian frame
    R = A / lattice_vectors(lattice_params(A)...)
    (n, θ) = matrix_to_axis_angle(R * det(R))
    @assert 0 <= θ <= π

    # Quantize so that floats can be compared without precision artifacts
    θ = round(θ; digits=12)
    n = round.(n; digits=12)
    T = round.(setting.T; digits=12)

    # Lexicographic score in decreasing priority
    return (
        -sign(det(A)),              # no reflection
        -(norm(T) < 1e-12),         # no translation
        -(abs(θ) < 1e-12),          # no rotation
        abs(θ),                     # small rotation
        norm(T),                    # small translation
        -abs.(reverse(n) .+ 1e-12), # break ties by favoring z-axis rotation, clockwise
        +abs.(reverse(T) .- 1e-12), # break ties by favoring small a3 translation, positive
    )
end

function conventionalize_setting(latvecs::Mat3, setting::SymOp, sgnum::Int)
    symops = SymOp.(Spglib.get_symmetry_from_database(standard_setting[sgnum])...)
    if length(symops) == 1
        return setting
    end

    # The `setting` is ambiguous because `s * setting` is another valid setting
    # for each symop `s`. Build list of candidates. Without loss of generality,
    # wrap components of translation vector T to (-0.5, 0.5].
    candidate_settings = unique(map(symops) do s
        (; R, T) = s * setting
        T = mod.(T .+ (0.5 - 1e-12), 1) .- (0.5 - 1e-12)
        SymOp(R, T)
    end)

    scores = score_setting.(candidate_settings; latvecs)
    p = sortperm(scores)

    if !allunique(scores)
        @warn "Could not select standard cell unambiguously"
    end

    return candidate_settings[p[1]]
end

function crystal_from_inferred_symmetry(latvecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec, suppress_warnings=false)
    validate_positions(positions; symprec)

    d = spg_get_dataset_scaled(Spglib.Cell(latvecs, positions, types), symprec)
    classes = d.crystallographic_orbits
    # classes = d.equivalent_atoms
    symops = SymOp.(d.rotations, d.translations)
    label = spacegroup_label(Int(d.hall_number))
    sgnum = Int(d.spacegroup_number)
    cell_std = cell_type(standard_setting[sgnum])
    setting = mapping_to_standard_setting_from_spglib_dataset(d)
    setting = conventionalize_setting(latvecs, setting, sgnum)

    # Spglib-inferred symops are unreliable for large cells. Empty the list
    # to prevent faulty symmetry analysis.
    volume_ratio = length(positions) // d.n_std_atoms
    volume_ratio > 1 && empty!(symops)
    volume_ratio_str = number_to_math_string(volume_ratio)

    if !suppress_warnings
        if volume_ratio > 1
            volume_ratio_str = number_to_math_string(volume_ratio)
            types_str = allunique(types) ? "" : " or distinct `types=[...]`"
            @error "Symmetry analysis disabled! Cell is $volume_ratio_str times too large. Fix with `standardize`$types_str."
        elseif volume_ratio < 1
            @info "Cell is $volume_ratio_str the standard size for spacegroup $sgnum. Consider `standardize`."
        elseif cell_std in (tetragonal, hexagonal, cubic) && setting.R ≉ I
            # For a tetragonal, hexagonal, trigonal or cubic space group, the
            # provided cell should already be Spglib-idealized.
            @info "Nonstandard $cell_std cell for spacegroup $sgnum. Consider `standardize`."
        elseif !(cell_type(latvecs) in all_compatible_cells(cell_std))
            # For triclinic, monoclinic, and orthorhombic spacegroups, check
            # that the provided cell is symmetry-consistent (not necessarily
            # Spglib-idealized). For orthorhombic, arbitrary axis ordering is
            # acceptable. For monoclinic, arbitrary β magnitude is acceptable,
            # whereas Spglib finds a cell with 90 < β < 120. For triclinic,
            # Spglib's idealization algorithm is even more complex, involving
            # Niggli reduction.
            @info "Nonstandard $cell_std cell for spacegroup $sgnum. Consider `standardize`."
        elseif det(latvecs) < 0
            @info "Lattice vectors are not right-handed. Consider `standardize`."
        end
    end

    sg = Spacegroup(symops, label, sgnum, setting)

    # If the spacegroup setting matches a Hall number then use its tabulated
    # spacegroup data. TODO: A better implementation would snap `sg.setting.T`
    # to zero and recalculate symops from table data.
    sg = idealize_spacegroup(sg; tol=symprec)

    # Idealize the basis vectors for the lattice system
    latvecs = idealize_latvecs(sg, latvecs; tol=symprec)
    recipvecs = 2π*Mat3(inv(latvecs)')

    # Idealize positions using Wyckoff tables
    positions = idealize_position.(Ref(sg), d.wyckoffs, positions; tol=symprec)

    # Renumber class indices so that they are ascending, from 1..max_class.
    classes = [findfirst(==(c), unique(classes)) for c in classes]
    @assert unique(classes) == 1:maximum(classes)

    ret = Crystal(nothing, latvecs, recipvecs, positions, types, classes, sg)
    validate_crystal(ret)

    return ret
end

function validate_positions(positions::Vector{Vec3}; symprec)
    for i in eachindex(positions), j in i+1:length(positions)
        ri, rj = positions[[i, j]]
        overlapping = is_periodic_copy(ri, rj; tol=1.001symprec)
        too_close = is_periodic_copy(ri, rj; tol=4.001symprec)
        if overlapping || too_close
            descriptor = overlapping ? "Overlapping" : "Near-overlapping"
            ri_str, rj_str = pos_to_string.((ri, rj))
            symprec_str = number_to_simple_string(symprec; digits=2)
            error("$descriptor positions $ri_str and $rj_str at symprec=$symprec_str")
        end
    end
end

function validate_orbits(positions::Vector{Vec3}, orbits::Vector{Vector{Vec3}}; symprec, wyckoffs=nothing)
    @assert size(positions) == size(orbits)
    # Check that orbits are distinct
    for i in eachindex(positions), j in i+1:length(positions)
        ri, rj = positions[[i, j]]
        overlapping = any(is_periodic_copy.(Ref(ri), orbits[j]; tol=1.001symprec))
        too_close = any(is_periodic_copy.(Ref(ri), orbits[j]; tol=4.001symprec))
        if overlapping || too_close
            descriptor = overlapping ? "Equivalent" : "Near-equivalent"
            ri_str, rj_str = pos_to_string.((ri, rj))
            symprec_str = number_to_simple_string(symprec; digits=2)
            if isnothing(wyckoffs)
                error("$descriptor positions $ri_str and $rj_str at symprec=$symprec_str")
            else
                wyckstr = wyckoff_string(wyckoffs[i])
                error("$descriptor positions $ri_str and $rj_str in Wyckoff $wyckstr at symprec=$symprec_str")
            end
        end
    end

    return orbits
end

function validate_crystal(cryst::Crystal)
    (; latvecs, positions, types, classes, sg) = cryst

    for i in eachindex(positions)
        # Atoms of the same class must have the same type
        for j in eachindex(positions)
            if classes[i] == classes[j]
                @assert types[i] == types[j]
            end
        end

        # Wyckoffs must have the correct multiplicity
        w = get_wyckoff(cryst, i)
        cell_multiplicity = count(==(classes[i]), classes)
        @assert w.multiplicity ≈ cell_multiplicity / abs(det(sg.setting.R))
    end

    # Symop rotation/reflection matrices R must be orthogonal in Cartesian
    # coordinates.
    for s in sg.symops
        R = latvecs * s.R * inv(latvecs)
        # Due to possible imperfections in the lattice vectors, only require
        # that R is approximately orthogonal
        @assert norm(R*R' - I) < 1e-12 "Lattice vectors and symmetry operations are incompatible."
    end

    # TODO: Check that space group is closed and that symops have inverse?
end

function validate_symops(cryst::Crystal)
    if isempty(cryst.sg.symops)
        error("Crystal is missing symmetry information; consider `standardize`.")
    end
end


# repeat_multiple(["a", "b"], [2, 3]) == ["a", "a", "b", "b", "b"]
function repeat_multiple(vals, lens)
    return reduce(vcat, map(vals, lens) do val, len
        fill(val, len)
    end)
end

# Builds a crystal from an explicit set of symmetry operations and a minimal set
# of positions
function crystal_from_spacegroup(latvecs::Mat3, positions::Vector{Vec3}, types::Vector{String}, sg::Spacegroup; symprec)
    latvecs = idealize_latvecs(sg, latvecs; tol=symprec)
    recipvecs = 2π*inv(latvecs)'

    wyckoffs = find_wyckoff.(Ref(sg), positions; tol=symprec)
    orbits = map(wyckoffs) do w
        # Transform Wyckoff expression into custom setting
        expr0 = transform(inv(sg.setting), w.expr)
        # Map orbit Vector{WyckoffExpr} to Vector{Vec3}
        map(crystallographic_orbit(expr0; sg.symops)) do (; F, c)
            wrap_to_unit_cell(F * w.θ + c)
        end
    end
    validate_orbits(positions, orbits; symprec, wyckoffs)

    all_positions = reduce(vcat, orbits)
    all_types = repeat_multiple(types, length.(orbits))
    all_classes = repeat_multiple(eachindex(orbits), length.(orbits))

    ret = Crystal(nothing, latvecs, recipvecs, all_positions, all_types, all_classes, sg)
    sort_sites!(ret)
    validate_crystal(ret)

    return ret
end


function get_wyckoff(cryst::Crystal, i::Int)
    return find_wyckoff(cryst.sg, cryst.positions[i]; tol=1e-12)
end


"""
    primitive_cell(cryst)

The shape of a primitive cell in multiples of the lattice vectors for `cryst`.
Specifically, columns of `cryst.latvecs * primitive(cryst)` define the primitive
lattice vectors in global Cartesian coordinates. Returns `nothing` if spacegroup
setting information is missing.

This function may be useful for constructing inputs to
[`reshape_supercell`](@ref).
"""
function primitive_cell(cryst::Crystal)
    cryst = @something cryst.root cryst
    (; number, setting) = cryst.sg

    # Primitive lattice vectors in units of ITA standard lattice vectors
    centering = centering_symbol(standard_setting[number])
    prim_shape = standard_primitive_basis[centering]

    # Lattice vectors of ITA standard setting:
    #     std_latvecs = cryst.latvecs * inv(setting.R)
    # Primitive lattice vectors:
    #     prim_latvecs = std_latvecs * prim_shape
    # Primitive shape units of cryst.latvecs
    #     return inv(cryst.latvecs) * prim_latvecs

    # Shortcut for above
    return setting.R \ prim_shape
end


function check_shape_commensurate(cryst, shape)
    prim_cell = primitive_cell(cryst)
    shape_in_prim = prim_cell \ shape
    if !all_integer(shape_in_prim; tol=1e-12)
        if prim_cell ≈ I
            error("Elements of `shape` must be integer. Received $shape.")
        else
            error("Elements of `primitive_cell(cryst) \\ shape` must be integer. Calculated $shape_in_prim.")
        end
    end
end


# Indices of atoms that wrap into unique positions within the primitive cell.
# The choice is ambiguous because the un-wrapped positions are not restricted to
# `latvecs * primitive_cell()`.
function primitive_atoms_arbitrary(cryst::Crystal)
    P = primitive_cell(cryst)
    P ≈ I && return collect(1:natoms(cryst))
    invP = inv(P)

    atoms = Int[]
    for (i, ri) in enumerate(cryst.positions)
        if all(rj -> !is_periodic_copy(invP * ri, invP * rj), cryst.positions[atoms])
            push!(atoms, i)
        end
    end

    @assert length(atoms) / natoms(cryst) ≈ abs(det(P))
    return atoms
end

function reshape_crystal_aux(cryst::Crystal, new_latvecs)
    # Atom indices within the primitive cell and associated positions
    prim_inds = primitive_atoms_arbitrary(cryst)
    prim_global_positions = Ref(cryst.latvecs) .* cryst.positions[prim_inds]

    # new_shape_in_prim defines the new supercell in multiples of primitive
    # cells. Since the primitive cell is the smallest discrete unit, all
    # elements of the shape matrix are integer.
    prim_shape = primitive_cell(cryst)
    prim_latvecs = cryst.latvecs * prim_shape
    new_shape_in_prim = prim_latvecs \ new_latvecs
    @assert all_integer(new_shape_in_prim; tol=1e-12)
    new_shape_in_prim = round.(Int, new_shape_in_prim)

    # Factorize shape matrix as lower triangular H (column Hermite normal form)
    # times unimodular U.
    H = MatInt.col_hermite(new_shape_in_prim)
    @assert istril(H) && all(diag(H) .> 0)

    # Create a grid of shifted copies of the primitive cell and map them into
    # fractional coordinates for new_latvecs. The diagonal elements of H
    # determine the appropriate grid dimensions.
    new_positions = Vec3[]
    for n in Iterators.product(0:H[1,1]-1, 0:H[2,2]-1, 0:H[3,3]-1)
        x = [new_latvecs \ (r + prim_latvecs * collect(n)) for r in prim_global_positions]
        append!(new_positions, wrap_to_unit_cell.(x))
    end

    ncopies = round(Int, abs(det(new_shape_in_prim)))
    @assert length(new_positions) == length(prim_inds) * ncopies
    new_types = repeat(cryst.types[prim_inds], ncopies)
    new_classes = repeat(cryst.classes[prim_inds], ncopies)

    return (; new_positions, new_types, new_classes)
end

function reshape_crystal(cryst::Crystal, new_shape::Mat3)
    # Check that desired lattice vectors are commensurate with root cell
    check_shape_commensurate(cryst, new_shape)

    # root.latvecs defines the fractional coordinate system used by new_shape
    root = @something cryst.root cryst

    # Lattice vectors of the new unit cell in global coordinates
    new_latvecs = root.latvecs * new_shape

    # Return this crystal if possible
    new_latvecs ≈ cryst.latvecs && return cryst

    # Calculate atoms in new cell
    (; new_positions, new_types, new_classes) = reshape_crystal_aux(cryst, new_latvecs)

    # Reciprocal lattice vectors of the new unit cell
    new_recipvecs = 2π*inv(new_latvecs)'

    # Note that latvecs_std = latvecs * R⁻¹ = new_latvecs * new_R⁻¹ is an
    # invariant quantity. Using new_latvecs = latvecs * new_shape, this implies
    # that R⁻¹ = new_shape * new_R⁻¹, or new_R = R * new_shape.
    sg_setting = SymOp(root.sg.setting.R * new_shape, root.sg.setting.T)

    # Symops no longer available
    new_sg = Spacegroup(SymOp[], root.sg.label, root.sg.number, sg_setting)

    ret = Crystal(root, new_latvecs, new_recipvecs, new_positions, new_types, new_classes, new_sg)
    sort_sites!(ret)
    validate_crystal(ret)

    return ret
end


"""
    subcrystal(cryst, types) :: Crystal

Filters sublattices of a `Crystal` by atom `types`, keeping the space group
unchanged.

    subcrystal(cryst, classes) :: Crystal

Filters sublattices of `Crystal` by equivalence `classes`, keeping the space
group unchanged.
"""
function subcrystal(cryst::Crystal, types::Vararg{String, N}) where N
    for s in types
        if !(s in cryst.types)
            error("types string '$s' is not present in crystal.")
        end
    end
    atoms = findall(in(types), cryst.types)
    classes = unique(cryst.classes[atoms])
    return subcrystal(cryst, classes...)
end

function subcrystal(cryst::Crystal, classes::Vararg{Int, N}) where N
    root = @something cryst.root cryst
    for c in classes
        if !(c in cryst.classes)
            error("Class '$c' is not present in crystal.")
        end
    end

    atoms = findall(in(classes), cryst.classes)
    new_positions = cryst.positions[atoms]
    new_types = cryst.types[atoms]
    new_classes = cryst.classes[atoms]

    if atoms != 1:maximum(atoms)
        @info "Atoms have been renumbered in subcrystal."
    end

    ret = Crystal(root, cryst.latvecs, cryst.recipvecs, new_positions, new_types, new_classes, cryst.sg)
    return ret
end

# Avoids ambiguity error
subcrystal(cryst::Crystal) = cryst


function Base.show(io::IO, cryst::Crystal)
    spg = isempty(cryst.sg.label) ? "" : "$(cryst.sg.label), "
    println(io, "Crystal($spg$(natoms(cryst)) atoms)")
end

function Base.show(io::IO, ::MIME"text/plain", cryst::Crystal)
    printstyled(io, "Crystal\n"; bold=true, color=:underline)
    println(io, "Spacegroup $(cryst.sg.label)")

    if is_standard_form(cryst.latvecs)
        (a, b, c, α, β, γ) = lattice_params(cryst.latvecs)
        @printf io "Lattice params a=%.4g, b=%.4g, c=%.4g, α=%.4g°, β=%.4g°, γ=%.4g°\n" a b c α β γ
    else
        println(io, formatted_matrix(number_to_math_string.(cryst.latvecs); prefix="Lattice vectors: "))
    end

    @printf io "Cell volume %.4g\n" cell_volume(cryst)

    for c in unique(cryst.classes)
        i = findfirst(==(c), cryst.classes)
        descr = String[]
        if cryst.types[i] != ""
            push!(descr, "Type '$(cryst.types[i])'")
        end

        w = get_wyckoff(cryst, i)
        push!(descr, "Wyckoff $(wyckoff_string(w)) (site sym. '$(w.sitesym)')")
        if isempty(descr)
            push!(descr, "Class $c")
        end
        println(io, join(descr, ", "), ":")

        for i in findall(==(c), cryst.classes)
            pos = pos_to_string(cryst.positions[i])
            println(io, "   $i. $pos")
        end
    end
end


#= Definitions of common crystals =#

function square_crystal(; a=1.0, c=10a)
    latvecs = lattice_vectors(a, a, c, 90, 90, 90)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    return cryst
end

function triangular_crystal(; a=1.0, c=10a)
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    positions = [[0, 0, 0]]
    cryst = Crystal(latvecs, positions)
    return cryst
end

function kagome_crystal(; a=1.0, c=10a)
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    positions = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
    cryst = Crystal(latvecs, positions)
    sort_sites!(cryst)
    return cryst
end

function hexagonal_crystal(; a=1.0, c=10a)
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    positions = [[0, 0, 0], [1/3, 2/3, 0]]
    cryst = Crystal(latvecs, positions)
    sort_sites!(cryst)
    return cryst
end

function cubic_crystal(; a=1.0)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [[0, 0, 0]]
    return Crystal(latvecs, positions)
end

function fcc_crystal(; a=1.0)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [[0, 0, 0]/2,
                 [1, 1, 0]/2,
                 [1, 0, 1]/2,
                 [0, 1, 1]/2]
    cryst = Crystal(latvecs, positions)
    sort_sites!(cryst)
    return cryst
end

function bcc_crystal(; a=1.0)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [[0, 0, 0]/2, [1, 1, 1]/2]
    return Crystal(latvecs, positions)
end

function diamond_crystal(; a=1.0)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    positions = [
        [0, 0, 0]/4,
        [2, 2, 0]/4,
        [1, 1, 1]/4,
        [3, 3, 1]/4,
        [2, 0, 2]/4,
        [0, 2, 2]/4,
        [3, 1, 3]/4,
        [1, 3, 3]/4,
    ]
    cryst = Crystal(latvecs, positions)
    sort_sites!(cryst)
    return cryst
end

function pyrochlore_crystal(; a=1.0)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    return Crystal(latvecs, [[0, 0, 0]], 227)
end

function hyperkagome_crystal(; a=1.0)
    latvecs = lattice_vectors(a, a, a, 90, 90, 90)
    # This shift is appropriate to the Ir atoms of a conventional Na4Ir3O8 cell,
    # but any other shift would yield the same symmetries
    # https://materials.springer.com/isp/crystallographic/docs/sd_1723194
    x = 0.141
    cryst = Crystal(latvecs, [[1/8, x, x+1/4]], 213)
    return cryst
end

function hcp_crystal(; a=1.0, c=sqrt(8/3)*a)
    latvecs = lattice_vectors(a, a, c, 90, 90, 120)
    positions = [(1/3,2/3,1/4), (2/3,1/3,3/4)]
    return Crystal(latvecs, positions)
end
