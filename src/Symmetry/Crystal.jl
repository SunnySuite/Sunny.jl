
"""
An object describing a crystallographic unit cell and its space group symmetry.
Constructors are as follows:


    Crystal(filename; override_symmetry=false, symprec=nothing)

Reads the crystal from a `.cif` file located at the path `filename`. If
`override_symmetry=true`, the spacegroup will be inferred based on atom
positions and the returned unit cell may be reduced in size. For an mCIF file,
the return value is the magnetic supercell, unless `override_symmetry=true`. If
a precision for spacegroup symmetries cannot be inferred from the CIF file, it
must be specified with `symprec`. The `latvecs` field of the returned `Crystal`
will be in units of angstrom.

    Crystal(latvecs, positions; types=nothing, symprec=1e-5)

Constructs a crystal from the complete list of atom positions `positions`, with
coordinates (between 0 and 1) in units of lattice vectors `latvecs`. Spacegroup
symmetry information is automatically inferred using the [Spglib
package](https://github.com/spglib/spglib) [1]. The optional parameter `types`
is a list of strings, one for each atom, and can be used to break
symmetry-equivalence between atoms.

    Crystal(latvecs, positions, spacegroup; types=nothing, choice=nothing, symprec=1e-5)

Builds a crystal using the symmetries of a `spacegroup`. One representative atom
must be specified for each occupied Wyckoff. The `spacegroup` may be specified
as a number 1..230 or as a string, e.g., Hermann–Mauguin or Hall symbol. If only
a spacegroup number is provided, the ITA standard setting [2] will be employed,
consistent with conventions from the [Bilbao crystallographic
server](https://www.cryst.ehu.es). If a spacegroup symbol is provided that
allows for multiple ITA settings, the possible disambiguations will be included
in an error message, and may involve a `choice` string.


# Examples

```julia
# Read a Crystal from a .cif file
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

# Build an equivalent crystal using a different origin choice.
positions = [[0, 0, 0]]
cryst = Crystal(latvecs, positions, 227; choice="1")
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
    symprec   :: Float64                 # Tolerance to imperfections in symmetry
end

# Constructs a crystal from the complete list of atom positions `positions`,
# representing fractions (between 0 and 1) of the lattice vectors `latvecs`.
# All symmetry information is automatically inferred.
function Crystal(latvecs, positions; types::Union{Nothing, Vector{String}}=nothing, symprec=1e-5)
    print_crystal_warnings(latvecs, positions)
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
            @warn """`setting` argument is deprecated! Omit to get the ITA standard setting. 
                     Alternatively, use `choice` instead."""
        else
            @warn "`setting` argument is deprecated! Use `choice` instead."
        end
        choice = setting
    end

    print_crystal_warnings(latvecs, positions)
    latvecs = Mat3(latvecs)
    positions = Vec3.(positions)
    if isnothing(types)
        types = fill("", length(positions))
    end

    sgt = unique_spacegroup_type(symbol, latvecs; choice)
    sg = Spacegroup(Int(sgt.hall_number))
    return crystal_from_spacegroup(latvecs, positions, types, sg; symprec)
end


function print_crystal_warnings(latvecs, positions)
    det(latvecs) < 0 && @warn "Lattice vectors are not right-handed."
    if length(positions) >= 100
        @info """This a very large crystallographic cell, which Sunny does not handle well.
                 If the intention is to model chemical inhomogeneity, the recommended procedure is as
                 follows: First, create a small unit cell with an idealized structure. Next, create
                 a perfectly periodic `System` of the desired size. Finally, use `to_inhomogeneous`
                 and related functions to design a system with the desired inhomogeneities."""
    end
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
end

# Sort the sites according to class and fractional coordinates. This is an
# internal function.
function sort_sites!(cryst::Crystal)
    function less_than(i, j)
        ci = cryst.classes[i]
        cj = cryst.classes[j]
        if ci != cj
            return ci < cj
        end
        ri = cryst.positions[i]
        rj = cryst.positions[j]
        for k = 3:-1:1
            if !isapprox(ri[k], rj[k], atol=10cryst.symprec)
                return ri[k] < rj[k]
            end
        end
        str1, str2 = fractional_vec3_to_string.((ri, rj))
        error("""Detected two very close atoms ($str1 and $str2).
                 If positions inferred from spacegroup, try increasing `symprec` parameter to `Crystal`.""")
    end
    p = sort(eachindex(cryst.positions), lt=less_than)
    permute_sites!(cryst, p)
end


# Attempt to find a permutation matrix P (with sign flips) such that latvecs*P
# has a more standard form.
function permute_to_standardize_lattice_vectors(latvecs)
    P = Mat3(I) # Iteratively build permutation matrix

    # Clip small matrix elements to zero
    latvecs = [abs(x) < 1e-12 ? 0 : x for x in latvecs]

    # Conventionally, a1 should be aligned with x direction
    a1, a2, a3 = eachcol(latvecs)
    if iszero(a2[2]) && iszero(a2[3]) && norm(a2) ≈ norm(a1)
        P = P * SA[0 1 0; 1 0 0; 0 0 1] # Permute (a1, a2)
    elseif iszero(a3[2]) && iszero(a3[3]) && norm(a3) ≈ norm(a1)
        P = P * SA[0 0 1; 0 1 0; 1 0 0] # Permute (a1, a3)
    end

    # Conventionally, a2 should be in xy plane
    _, a2, a3 = eachcol(latvecs * P)
    if iszero(a3[3]) && norm(a3) ≈ norm(a2)
        P = P * SA[1 0 0; 0 0 1; 0 1 0] # Permute (a2, a3)
    end

    # Flip columns so that the diagonal elements are positive
    signs = [a < 0 ? -1 : 1 for a in diag(latvecs*P)]
    P = P * Diagonal(signs)

    # To preserve volume-orientation, we may need to flip one lattice
    # vector. Pick a2 arbitrarily.
    @assert det(P) ≈ 1 || det(P) ≈ -1
    if det(P) ≈ -1
        P = P * Diagonal(SA[1,-1,1])
    end

    @assert det(P) ≈ 1
    @assert P' ≈ inv(P)
    return P
end


"""
    standardize(cryst::Crystal; idealize=true)

Return the symmetry-inferred standardized crystal unit cell. If `idealize=true`,
then the lattice vectors and site positions will be adapted. See "definitions
and conventions" of the [spglib
documentation](https://spglib.readthedocs.io/en/stable/) for more information.
"""
function standardize(cryst::Crystal; idealize=true)
    !isnothing(cryst.root) && error("Call this function on root crystal instead")

    (; symprec) = cryst
    cell = Spglib.Cell(cryst.latvecs, cryst.positions, cryst.types)
    (; lattice, positions, atoms) = Spglib.standardize_cell(cell, symprec; no_idealize=!idealize)
    positions = Vec3.(positions)
    lattice = Mat3(lattice)

    if !idealize
        # Empirically, these lattice vectors from spglib are accurate to about 6
        # digits. However, spglib produces much higher accuracy with the
        # `idealize=true` option. Rotate the higher precision lattice vectors so
        # that they give the best match to the ones for the non-idealized cell.
        std_lattice = Mat3(Spglib.standardize_cell(cell, symprec; no_idealize=false).lattice)
        R = closest_unitary(lattice / std_lattice)
        isapprox(R*std_lattice, lattice; rtol=1e-5) || error("Failed to standardize the cell")
        lattice = R * std_lattice
        # In the non-idealized case, the spglib choice of lattice vectors can
        # sometimes be strange. For example, in a tetrahedral cell, (a1, a2)
        # might be pointing along (y, -x), whereas (x, y) would be a more
        # natural choice. Attempt to permute lattice vectors back to a standard
        # order, with sign-flips as needed.
        P = permute_to_standardize_lattice_vectors(lattice)
        # These transformations preserve global positions, `lattice * r`
        lattice = lattice * P
        positions = [P' * r for r in positions]
    end

    ret = crystal_from_inferred_symmetry(lattice, positions, atoms; symprec)
    sort_sites!(ret)

    # TODO: Make this case work by avoiding Spglib.standardize_cell and using
    # sg.setting data instead. It might come up if one is using an ITA setting
    # that is not standard.
    if ret.sg.number != cryst.sg.number
        error("Inferred spacegroup $(ret.sg.number) doesn't match $(cryst.sg.number); are any atoms missing from chemical cell?")
    end

    return ret
end

function crystal_from_inferred_symmetry(latvecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec, check_cell=true)
    # Print a warning if non-conventional lattice vectors are detected.
    try cell_type(latvecs) catch e @warn e.msg end

    for i in 1:length(positions)
        for j in i+1:length(positions)
            ri = positions[i]
            rj = positions[j]
            if all_integer(ri-rj; symprec)
                error("Positions $ri and $rj are equivalent.")
            end
        end
    end

    recipvecs = 2π*Mat3(inv(latvecs)')
    positions = wrap_to_unit_cell.(positions; symprec)

    cell = Spglib.Cell(latvecs, positions, types)
    d = Spglib.get_dataset(cell, symprec)

    if check_cell
        ratio = length(positions) / d.n_std_atoms
        if ratio > 1
            ratio_str = number_to_simple_string(ratio; digits=4)
            types_str = allunique(types) ? "" : " Alternatively, break site symmetry with `types=[...]`."
            @warn """The symmetry-inferred, conventional unit cell is $ratio_str times smaller. Obtain
                     it with `standardize`.$types_str"""
        end
    end

    classes = d.crystallographic_orbits
    # classes = d.equivalent_atoms
    symops = SymOp.(d.rotations, d.translations)
    label = spacegroup_label(Int(d.hall_number))
    number = d.spacegroup_number
    setting = mapping_to_standard_setting_from_spglib_dataset(d)
    sg = Spacegroup(symops, label, number, setting)

    # Renumber class indices so that they are ascending, from 1..max_class.
    classes = [findfirst(==(c), unique(classes)) for c in classes]
    @assert unique(classes) == 1:maximum(classes)

    ret = Crystal(nothing, latvecs, recipvecs, positions, types, classes, sg, symprec)
    validate(ret)
    for i in 1:natoms(ret)
        w = get_wyckoff(ret, i)
        @assert w.letter == d.wyckoffs[i]
        @assert w.sitesym == d.site_symmetry_symbols[i]
    end

    return ret
end


function is_spacegroup_type_consistent(sgt, latvecs)
    cell = cell_type(latvecs)
    hall_cell = cell_type(Int(sgt.hall_number))

    if hall_cell == Sunny.monoclinic
        # Special handling of monoclinic space groups. There are three possible
        # conventions for the unit cell, depending on which of α, β, or γ is
        # special.
        _, _, _, α, β, γ = lattice_params(latvecs)
        x = first(replace(sgt.choice, "-" => ""))
        if x == 'a'
            return β≈90 && γ≈90
        elseif x == 'b'
            return α≈90 && γ≈90
        elseif x == 'c'
            return α≈90 && β≈90
        end
    else
        return cell in all_compatible_cells(hall_cell)
    end
end


function crystallographic_orbit(symops::Vector{SymOp}, position::Vec3; symprec)
    orbit = Vec3[]
    for s = symops
        x = wrap_to_unit_cell(transform(s, position); symprec)
        if !any(y -> all_integer(x-y; symprec), orbit)
            push!(orbit, x)
        end
    end
    return orbit
end

function crystallographic_orbits_distinct(symops::Vector{SymOp}, positions::Vector{Vec3}; symprec, wyckoffs=nothing)
    orbits = crystallographic_orbit.(Ref(symops), positions; symprec)

    # Check orbits are distinct
    for (i, ri) in enumerate(positions), j in i+1:length(orbits)
        if any(rj -> all_integer(ri-rj; symprec), orbits[j])
            if isnothing(wyckoffs)
                error("Reference positions $(positions[i]) and $(positions[j]) are symmetry equivalent.")
            else
                @assert wyckoffs[i].letter == wyckoffs[j].letter
                (; multiplicity, letter) = wyckoffs[i]
                error("Reference positions $(positions[i]) and $(positions[j]) are symmetry equivalent in Wyckoff $multiplicity$letter.")
            end
        end
    end

    return orbits
end

function repeat_multiple(vals, lens)
    reduce(vcat, map(vals, lens) do val, len
        fill(val, len)
    end)
end

# Builds a crystal from an explicit set of symmetry operations and a minimal set
# of positions
function crystal_from_spacegroup(latvecs::Mat3, positions::Vector{Vec3}, types::Vector{String}, sg::Spacegroup; symprec)
    wyckoffs = find_wyckoff_for_position.(Ref(sg), positions; symprec)
    orbits = crystallographic_orbits_distinct(sg.symops, positions; symprec, wyckoffs)

    # Symmetry-propagated orbits must match known multiplicities of the Wyckoffs
    foreach(orbits, wyckoffs) do orbit, wyckoff
        @assert wyckoff.multiplicity ≈ length(orbit) / abs(det(sg.setting.R))
    end

    all_positions = reduce(vcat, orbits)
    all_types = repeat_multiple(types, length.(orbits))
    all_classes = repeat_multiple(eachindex(orbits), length.(orbits))

    recipvecs = 2π*Mat3(inv(latvecs)')
    ret = Crystal(nothing, latvecs, recipvecs, all_positions, all_types, all_classes, sg, symprec)
    sort_sites!(ret)
    validate(ret)

    return ret
end

function get_wyckoff(cryst::Crystal, i::Int)
    (; classes, positions, sg, symprec) = cryst
    wyckoff = find_wyckoff_for_position(sg, positions[i]; symprec)
    cell_multiplicity = count(==(classes[i]), classes)
    @assert wyckoff.multiplicity ≈ cell_multiplicity / abs(det(sg.setting.R))
    return wyckoff
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
    prim_cell = @something primitive_cell(cryst) Mat3(I)
    shape_in_prim = prim_cell \ shape
    if !all_integer(shape_in_prim; cryst.symprec)
        if prim_cell ≈ I
            error("Elements of `shape` must be integer. Received $shape.")
        else
            error("Elements of `primitive_cell(cryst) \\ shape` must be integer. Calculated $shape_in_prim.")
        end
    end
end

function reshape_crystal(cryst::Crystal, new_shape::Mat3)
    # Check that desired lattice vectors are commensurate with root cell
    check_shape_commensurate(cryst, new_shape)

    # The `root.latvecs` defines the fractional coordinate system, but note that
    # `cryst` may be formed as a subcrystal of `root`. For reshaping purposes,
    # therefore, we must use `cryst.positions` and not `root.positions`, etc.
    root = @something cryst.root cryst

    # Lattice vectors of the new unit cell in global coordinates
    new_latvecs = root.latvecs * new_shape

    # Return this crystal if possible
    new_latvecs ≈ cryst.latvecs && return cryst

    # Symmetry precision needs to be rescaled for the new unit cell. Ideally we
    # would have three separate rescalings (one per lattice vector), but we're
    # forced to pick just one. Scale according to volume change.
    new_symprec = root.symprec / cbrt(abs(det(new_shape)))

    # This matrix defines a mapping from fractional coordinates `x` in `cryst`
    # to fractional coordinates `y` in the new unit cell.
    B = new_latvecs \ cryst.latvecs

    # Our goal is to loop over `cryst` cells (n1, n2, n3) and make sure we
    # completely cover `new_latvecs` cell, {0 ≤ y₁ ≤ 1, 0 ≤ y₂ ≤ 1, 0 ≤ y₃ ≤ 1}.
    # Due to the linear relationship `y = Bx`, the required range of (n1, n2,
    # n3) is associated with inv(B). For example, `inv(B) * [1, 1, 1]` should be
    # covered in the space of `x` sampling.
    #
    # It is not clear to me what the right formula is, and the
    # `sum.(eachrow(...))` formula below is a conservative, heuristic guess.
    # Fortunately, it is well protected against mistakes via the check on atom
    # count. Any mistake will yield an assertion error: "Missing atoms in
    # reshaped unit cell".
    nmax = round.(Int, sum.(eachrow(abs.(inv(B))))) .+ 1

    new_positions = Vec3[]
    new_types     = String[]
    new_classes   = Int[]

    for i in 1:natoms(cryst)
        for n1 in -nmax[1]:nmax[1], n2 in -nmax[2]:nmax[2], n3 in -nmax[3]:nmax[3]
            x = cryst.positions[i] + Vec3(n1, n2, n3)
            y = B*x

            # Check whether the new position y (in fractional coordinates
            # associated with `new_latvecs`) is within the new unit cell.
            # Account for finite symprec ϵ by checking the bounds [-ϵ,1-ϵ). See
            # related comment in `wrap_to_unit_cell`.
            if all(-new_symprec .<= y .< 1 - new_symprec)
                push!(new_positions, y)
                push!(new_types, cryst.types[i])
                push!(new_classes, cryst.classes[i])
            end
        end
    end

    # Check that we have the right number of atoms
    @assert abs(det(B)) * length(new_positions) ≈ natoms(cryst) "Missing atoms in reshaped unit cell. Please report this bug!"

    # Reciprocal lattice vectors of the new unit cell
    new_recipvecs = 2π * Mat3(inv(new_latvecs)')

    # Note that latvecs_std = latvecs * R⁻¹ = new_latvecs * new_R⁻¹ is an
    # invariant quantity. Using new_latvecs = latvecs * new_shape, this implies
    # that R⁻¹ = new_shape * new_R⁻¹, or new_R = R * new_shape.
    sg_setting = SymOp(root.sg.setting.R * new_shape, root.sg.setting.T)

    # Empty symops list indicates that this information has been lost.
    new_sg = Spacegroup(SymOp[], root.sg.label, root.sg.number, sg_setting)

    return Crystal(root, new_latvecs, new_recipvecs, new_positions, new_types, new_classes, new_sg, new_symprec)
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

    ret = Crystal(root, cryst.latvecs, cryst.recipvecs, new_positions, new_types, new_classes, cryst.sg, cryst.symprec)
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
        (; multiplicity, letter, sitesym) = get_wyckoff(cryst, i)
        push!(descr, "Wyckoff $multiplicity$letter (site sym. '$sitesym')")
        if isempty(descr)
            push!(descr, "Class $c")
        end
        println(io, join(descr, ", "), ":")

        for i in findall(==(c), cryst.classes)
            pos = fractional_vec3_to_string(cryst.positions[i])
            println(io, "   $i. $pos")
        end
    end
end

function validate(cryst::Crystal)
    # Atoms of the same class must have the same type
    for i in eachindex(cryst.positions)
        for j in eachindex(cryst.positions)
            if cryst.classes[i] == cryst.classes[j]
                @assert cryst.types[i] == cryst.types[j]
            end
        end
    end

    # Rotation matrices in global coordinates must be orthogonal
    for s in cryst.sg.symops
        R = cryst.latvecs * s.R * inv(cryst.latvecs)
        # Due to possible imperfections in the lattice vectors, only require
        # that R is approximately orthogonal
        @assert norm(R*R' - I) < cryst.symprec "Lattice vectors and symmetry operations are incompatible."
    end

    # TODO: Check that space group is closed and that symops have inverse?
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
