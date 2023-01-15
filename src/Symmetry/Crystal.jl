struct SiteSymmetry
    symbol       :: String
    multiplicity :: Int
    wyckoff      :: Char
end


"""
A `Crystal` object describes a crystallographic unit cell, and contains
information about its space group symmetry. Constructors are as follows:


    Crystal(lat_vecs, positions; types=nothing, symprec=1e-5)

Constructs a crystal from the complete list of atom positions `positions`,
representing fractions (between 0 and 1) of the lattice vectors `lat_vecs`. All
symmetry information is automatically inferred.

    Crystal(lat_vecs, positions, symbol::String; types=nothing, setting=nothing, symprec=1e-5)

Builds a crystal by applying the symmetry operators for a given spacegroup
symbol.

    Crystal(lat_vecs, positions, spacegroup_number; types=nothing, setting=nothing, symprec=1e-5)

Builds a crystal by applying symmetry operators for a given international
spacegroup number.

    Crystal(filename::AbstractString; symprec=1e-5)

Reads the crystal from a `.cif` file located at the path `filename`.
"""
struct Crystal
    lat_vecs       :: Mat3                                 # Lattice vectors as columns
    positions      :: Vector{Vec3}                         # Positions in fractional coords
    types          :: Vector{String}                       # Types
    classes        :: Vector{Int}                          # Class indices
    sitesyms       :: Union{Nothing, Vector{SiteSymmetry}} # Optional site symmetries
    symops         :: Vector{SymOp}                        # Symmetry operations
    spacegroup     :: String                               # Description of space group
    symprec        :: Float64                              # Tolerance to imperfections in symmetry
end

"""
    nbasis(crystal::Crystal)

Number of basis positions (sublattices) in the unit cell.
"""
nbasis(cryst::Crystal) = length(cryst.positions)

"""
    cell_volume(crystal::Crystal)

Volume of the crystal unit cell.
"""
cell_volume(cryst::Crystal) = abs(det(cryst.lat_vecs))

"""
    position(crystal::Crystal, i::Int, cell=(0,0,0))

Position of an atom in global Cartesian coordinates. The optional `cell`
parameter denotes a displacement in unit cell indices.
"""
position(cryst::Crystal, i::Int, offset=(0,0,0)) = cryst.lat_vecs * (convert(Vec3, offset) + cryst.positions[i])


# Constructs a crystal from the complete list of atom positions `positions`,
# representing fractions (between 0 and 1) of the lattice vectors `lat_vecs`.
# All symmetry information is automatically inferred.
function Crystal(lat_vecs, positions; types::Union{Nothing,Vector{String}}=nothing, symprec=1e-5)
    lat_vecs = convert(Mat3, lat_vecs)
    positions = [convert(Vec3, p) for p in positions]
    if isnothing(types)
        types = fill("", length(positions))
    end
    return crystal_from_inferred_symmetry(lat_vecs, positions, types; symprec)
end


# Builds a crystal by applying the symmetry operators for a given spacegroup
# symbol.
function Crystal(lat_vecs, positions, symbol::String; types::Union{Nothing,Vector{String}}=nothing, setting=nothing, symprec=1e-5)
    lat_vecs = convert(Mat3, lat_vecs)
    positions = [convert(Vec3, p) for p in positions]
    if isnothing(types)
        types = fill("", length(positions))
    end
    crystal_from_symbol(lat_vecs, positions, types, symbol; setting, symprec)
end

# Builds a crystal by applying symmetry operators for a given international
# spacegroup number.
function Crystal(lat_vecs, positions, spacegroup_number::Int; types::Union{Nothing,Vector{String}}=nothing, setting=nothing, symprec=1e-5)
    lat_vecs = convert(Mat3, lat_vecs)
    positions = [convert(Vec3, p) for p in positions]
    if isnothing(types)
        types = fill("", length(positions))
    end
    symbol = string(spacegroup_number)
    crystal_from_symbol(lat_vecs, positions, types, symbol; setting, symprec)
end


function spacegroup_name(hall_number::Int)
    # String representation of space group
    sgt = Spglib.get_spacegroup_type(hall_number)
    return "HM symbol '$(sgt.international)' ($(sgt.number))"
end

function symops_from_spglib(rotations, translations)
    Rs = Mat3.(transpose.(eachslice(rotations, dims=3)))
    Ts = Vec3.(eachcol(translations))
    return SymOp.(Rs, Ts)
end


# Sort the sites according to class and fractional coordinates.
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
            if !isapprox(ri[k], rj[k], atol=cryst.symprec)
                return ri[k] < rj[k]
            end
        end
        error("Positions $i and $j cannot be distinguished.")
    end
    perm = sort(eachindex(cryst.positions), lt=less_than)
    cryst.positions .= cryst.positions[perm]
    cryst.classes .= cryst.classes[perm]
    cryst.types .= cryst.types[perm]
end

# Wrap each coordinate of position r into the range [0,1). To account for finite
# precision, wrap 1-ϵ to -ϵ, where ϵ=symprec is a tolerance parameter.
function wrap_to_unit_cell(r::Vec3; symprec)
    return @. mod(r+symprec, 1) - symprec
end

function all_integer(x; symprec)
    return norm(x - round.(x)) < symprec
end

function crystal_from_inferred_symmetry(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec=1e-5)
    for i in 1:length(positions)
        for j in i+1:length(positions)
            ri = positions[i]
            rj = positions[j]
            if all_integer(ri-rj; symprec)
                error("Positions $ri and $rj are symmetry equivalent.")
            end
        end
    end

    positions = wrap_to_unit_cell.(positions; symprec)

    cell = Spglib.Cell(lat_vecs, positions, types)
    d = Spglib.get_dataset(cell, symprec)
    classes = d.crystallographic_orbits
    # classes = d.equivalent_atoms
    symops = symops_from_spglib(d.rotations, d.translations)
    spacegroup = spacegroup_name(d.hall_number)

    # renumber class indices so that they go from 1:max_class
    classes = [findfirst(==(c), unique(classes)) for c in classes]
    @assert unique(classes) == 1:maximum(classes)

    # multiplicities for the equivalence classes
    multiplicities = map(classes) do c
        # indices that belong to class c
        idxs = findall(==(c), classes)
        # indices in the primitive cell that belong to class c
        prim_idxs = unique(d.mapping_to_primitive[idxs])
        # number of atoms in the standard cell that correspond to each primitive index for class c
        counts = [count(==(i), d.std_mapping_to_primitive) for i in prim_idxs]
        # sum over all equivalent atoms in the primitive cell
        sum(counts)
    end

    sitesyms = SiteSymmetry.(d.site_symmetry_symbols, multiplicities, d.wyckoffs)

    ret = Crystal(lat_vecs, positions, types, classes, sitesyms, symops, spacegroup, symprec)
    validate(ret)
    return ret
end


# Build Crystal using the space group denoted by a unique Hall number. The complete
# list is given at http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
function crystal_from_hall_number(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}, hall_number::Int; symprec=1e-5)
    cell = cell_type(lat_vecs)
    hall_cell = cell_type(hall_number)
    allowed_cells = all_compatible_cells(hall_cell)
    @assert cell in allowed_cells "Hall number $hall_number requires a $hall_cell cell, but found $cell."

    if hall_cell == Sunny.monoclinic
        is_compatible = is_compatible_monoclinic_cell(lat_vecs, hall_number)
        @assert is_compatible "Lattice vectors define a monoclinic cell that is incompatible with Hall number $hall_number."
    end

    symops = symops_from_spglib(Spglib.get_symmetry_from_database(hall_number)...)
    spacegroup = spacegroup_name(hall_number)

    return crystal_from_symops(lat_vecs, positions, types, symops, spacegroup; symprec)
end

function crystal_from_symbol(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}, symbol::String; setting=nothing, symprec=1e-5)
    hall_numbers = Int[]
    crysts = Crystal[]

    n_hall_numbers = 530
    for hall_number in 1:n_hall_numbers
        sgt = Spglib.get_spacegroup_type(hall_number)

        if (replace(symbol, " "=>"") == sgt.international_short || 
            symbol in [string(sgt.number), sgt.hall_symbol, sgt.international, sgt.international_full])

            # Some Hall numbers may be incompatible with unit cell of provided
            # lattice vectors; skip them.
            is_compatible = true

            cell = cell_type(lat_vecs)
            hall_cell = cell_type(hall_number)
            allowed_cells = all_compatible_cells(hall_cell)

            # Special handling of trigonal space groups
            if Sunny.is_trigonal_symmetry(hall_number)
                # Trigonal symmetry must have either hexagonal or rhombohedral
                # cell, according to the Hall number.
                is_latvecs_valid = cell in [Sunny.rhombohedral, Sunny.hexagonal]
                if !is_latvecs_valid
                    error("Symbol $symbol requires a rhomobohedral or hexagonal cell, but found $cell.")
                end
                is_compatible = cell in allowed_cells
            else
                # For all other symmetry types, there is a unique cell for each Hall number
                if !(cell in allowed_cells)
                    error("Symbol $symbol requires a $hall_cell cell, but found $cell.")
                end
            end

            if hall_cell == Sunny.monoclinic
                is_compatible = is_compatible_monoclinic_cell(lat_vecs, hall_number)
            end

            if is_compatible
                cryst = crystal_from_hall_number(lat_vecs, positions, types, hall_number; symprec)
                push!(hall_numbers, hall_number)
                push!(crysts, cryst)
            end
        end
    end

    if length(crysts) == 0
        error("Could not find symbol '$symbol' in database.")
    elseif length(crysts) == 1
        return first(crysts)
    else
        if !isnothing(setting)
            i = findfirst(hall_numbers) do hall_number
                sgt = Spglib.get_spacegroup_type(hall_number)
                setting == sgt.choice
            end
            if isnothing(i)
                error("The symbol '$symbol' is ambiguous, and the specified setting '$setting' is not valid.")
            else
                return crysts[i]
            end
        end

        println("The spacegroup '$symbol' allows for multiple settings!")
        println("Returning a list of the possible crystals:")
        for (i, (hall_number, c)) in enumerate(zip(hall_numbers, crysts))
            sgt = Spglib.get_spacegroup_type(hall_number)
            hm_symbol = sgt.international
            choice = sgt.choice
            n_atoms = length(c.positions)
            i_str = @sprintf "%2d" i
            natoms_str = @sprintf "%2d" n_atoms
            println("   $i_str. \"$hm_symbol\", setting=\"$choice\", with $natoms_str atoms")
        end
        println()
        println("Note: To disambiguate, you may pass a named parameter, setting=\"...\".")
        println()
        return crysts
    end
end

"Builds a crystal from an explicit set of symmetry operations and a minimal set of positions "
function crystal_from_symops(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}, symops::Vector{SymOp}, spacegroup::String; symprec=1e-5)
    all_positions = Vec3[]
    all_types = String[]
    classes = Int[]
    
    for i = eachindex(positions)
        for s = symops
            x = wrap_to_unit_cell(transform(s, positions[i]); symprec)

            idx = findfirst(y -> all_integer(x-y; symprec), all_positions)
            if isnothing(idx)
                push!(all_positions, x)
                push!(all_types, types[i])
                push!(classes, i)
            else
                j = classes[idx]
                if i != j
                    error("Reference positions $(positions[i]) and $(positions[j]) are symmetry equivalent.")
                end
            end
        end
    end

    # Atoms are sorted by contiguous equivalence classes: 1, 2, ..., n
    @assert unique(classes) == 1:maximum(classes)

    # Ask Spglib to infer the spacegroup for the given positions and types
    inferred = crystal_from_inferred_symmetry(lat_vecs, all_positions, all_types; symprec)

    # Compare the inferred symops to the provided ones
    is_subgroup = all(symops) do s
        any(inferred.symops) do s′
            isapprox(s, s′; atol=symprec)
        end
    end
    is_supergroup = all(inferred.symops) do s
        any(symops) do s′
            isapprox(s, s′; atol=symprec)
        end
    end

    if !is_subgroup
        println("""Warning: User provided symmetry operation could not be inferred by Spglib,
                   which likely indicates a non-conventional unit cell.""")
    end

    # If the inferred symops match the provided ones, then we use the inferred
    # Crystal. Otherwise we must construct a new Crystal and do not have site
    # symmetry information.
    ret = if is_subgroup && is_supergroup
        inferred
    else
        Crystal(lat_vecs, all_positions, all_types, classes, nothing, symops, spacegroup, symprec)
    end
    sort_sites!(ret)
    validate(ret)
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
    idxs = findall(in(types), cryst.types)
    classes = unique(cryst.classes[idxs])
    return subcrystal(cryst, classes...)
end

function subcrystal(cryst::Crystal, classes::Vararg{Int, N}) where N
    for c in classes
        if !(c in cryst.classes)
            error("Class '$c' is not present in crystal.")
        end
    end

    idxs = findall(in(classes), cryst.classes)
    new_positions = cryst.positions[idxs]
    new_types = cryst.types[idxs]
    new_classes = cryst.classes[idxs]
    new_sitesyms = cryst.sitesyms[idxs]

    if idxs != 1:maximum(idxs)
        println("Warning: atoms are being renumbered.")
    end

    ret = Crystal(cryst.lat_vecs, new_positions, new_types, new_classes, new_sitesyms,
                  cryst.symops, cryst.spacegroup, cryst.symprec)
    return ret
end


function Base.show(io::IO, ::MIME"text/plain", cryst::Crystal)
    printstyled(io, "Crystal\n"; bold=true, color=:underline)
    println(io, cryst.spacegroup)

    if is_standard_form(cryst.lat_vecs)
        (a, b, c, α, β, γ) = lattice_params(cryst.lat_vecs)
        @printf io "Lattice params a=%.4g, b=%.4g, c=%.4g, α=%.4g°, β=%.4g°, γ=%.4g°\n" a b c α β γ
    else
        println(io, "Lattice vectors:")
        for a in eachcol(cryst.lat_vecs)
            @printf io "   [%.4g %.4g %.4g]\n" a[1] a[2] a[3]
        end
    end

    @printf io "Cell volume %.4g\n" cell_volume(cryst)

    for c in unique(cryst.classes)
        i = findfirst(==(c), cryst.classes)
        descr = String[]
        if cryst.types[i] != ""
            push!(descr, "Type '$(cryst.types[i])'")
        end
        if !isnothing(cryst.sitesyms)
            symbol = cryst.sitesyms[i].symbol
            multiplicity = cryst.sitesyms[i].multiplicity
            wyckoff = cryst.sitesyms[i].wyckoff
            push!(descr, "Wyckoff $multiplicity$wyckoff (point group '$symbol')")
        end
        println(io, join(descr, ", "), ":")

        for i in findall(==(c), cryst.classes)
            pos = atom_pos_to_string(cryst.positions[i])
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
    for s in cryst.symops
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        # Due to possible imperfections in the lattice vectors, only require
        # that R is approximately orthogonal
        @assert norm(R*R' - I) < cryst.symprec "Lattice vectors and symmetry operations are incompatible."
    end

    # TODO: Check that space group is closed and that symops have inverse?
end

#= Definitions of common crystals =#

function cubic_crystal(; a=1.0)
    lat_vecs = lattice_vectors(a, a, a, 90, 90, 90)
    basis_vecs = [[0, 0, 0]]
    Crystal(lat_vecs, basis_vecs)
end

function fcc_crystal(; a=1.0)
    lat_vecs = lattice_vectors(a, a, a, 90, 90, 90)
    basis_vecs = [[0, 0, 0]/2,
                  [1, 1, 0]/2,
                  [1, 0, 1]/2,
                  [0, 1, 1]/2]
    cryst = Crystal(lat_vecs, basis_vecs)
    sort_sites!(cryst)
    cryst
end

function fcc_primitive_crystal(; a=1.0)
    lat_vecs = [1 1 0; 1 0 1; 0 1 1]' * a/2
    basis_vecs = [[0, 0, 0]]
    Crystal(lat_vecs, basis_vecs)
end

function bcc_crystal(; a=1.0)
    lat_vecs = lattice_vectors(a, a, a, 90, 90, 90)
    basis_vecs = [[0, 0, 0]/2,
                  [1, 1, 1]/2,]
    Crystal(lat_vecs, basis_vecs)
end

function bcc_primitive_crystal(; a=1.0)
    lat_vecs = [1 1 -1; 1 -1 1; -1 1 1]' * a/2
    basis_vecs = [[0, 0, 0]]
    Crystal(lat_vecs, basis_vecs)
end


function diamond_crystal(; a=1.0)
    lat_vecs = lattice_vectors(a, a, a, 90, 90, 90)
    basis_vecs = [
        [0, 0, 0]/4,
        [2, 2, 0]/4,
        [1, 1, 1]/4,
        [3, 3, 1]/4,
        [2, 0, 2]/4,
        [0, 2, 2]/4,
        [3, 1, 3]/4,
        [1, 3, 3]/4,
    ]
    cryst = Crystal(lat_vecs, basis_vecs)
    sort_sites!(cryst)
    cryst
end

function diamond_primitive_crystal(; a=1.0)
    lat_vecs = [1 1 0; 1 0 1; 0 1 1]' * a/2
    basis_vecs = [
        [0, 0, 0]/4,
        [1, 1, 1]/4,
    ]
    Crystal(lat_vecs, basis_vecs)
end

function pyrochlore_lattice(; a=1.0)
    lat_vecs = [1 1 0; 1 0 1; 0 1 1]' * a/2
    basis_vecs = [
        [5, 5, 5]/8,
        [1, 5, 5]/8,
        [5, 5, 1]/8,
        [5, 1, 5]/8
    ]
    cryst = Crystal(lat_vecs, basis_vecs)
    sort_sites!(cryst)
    cryst
end