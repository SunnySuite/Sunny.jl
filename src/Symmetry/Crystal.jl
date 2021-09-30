"""Defines the SymOp, Crystal types"""

"""
    SymOp

Defines a symmetry operation belonging to a 3D space group, operating on fractional coordinates.
"""
struct SymOp
    R::Mat3
    T::Vec3
end

function Base.isapprox(s1::SymOp, s2::SymOp; atol)
    return isapprox(s1.R, s2.R; atol) && isapprox(s1.T, s2.T; atol)
end

function transform(s::SymOp, r::Vec3)
    return s.R*r + s.T
end


# Wrap each coordinate of position r into the range [0,1). To account for finite
# precision, wrap 1-ϵ to -ϵ, where ϵ=symprec is a tolerance parameter.
function wrap_to_unit_cell(r::Vec3; symprec=1e-5)
    return @. mod(r+symprec, 1) - symprec
end

function is_same_position(x, y; symprec=1e-5)
    return norm(rem.(x-y, 1, RoundNearest)) < symprec
end

struct SiteSymmetry
    symbol       :: String
    multiplicity :: Int
    wyckoff      :: Char
end


"""
    Crystal

A type holding all geometry and symmetry information needed to represent
 a three-dimensional crystal.
"""
struct Crystal
    lat_vecs       :: Mat3                # Lattice vectors as columns
    positions      :: Vector{Vec3}        # Full set of atoms, fractional coords
    classes        :: Vector{Int}         # Index to equivalent atom type
    types          :: Vector{String}      # Types for each atom
    symops         :: Vector{SymOp}       # Symmetry operations
    spacegroup     :: String              # Description of space group
    sitesymmetries :: Union{Nothing, Vector{SiteSymmetry}} # Optional site symmetries, indexed by class
    symprec        :: Float64             # Tolerance to imperfections in symmetry
end

nbasis(cryst::Crystal) = length(cryst.positions)
cell_volume(cryst::Crystal) = abs(det(cryst.lat_vecs))
lattice_params(cryst::Crystal) = lattice_params(cryst.lat_vecs)
lattice_vectors(cryst::Crystal) = cryst.lat_vecs


"""
    Crystal(lat_vecs, positions; types=nothing, symprec=1e-5)

Construct a `Crystal` using explicit geometry information, with all symmetry
information automatically inferred. `positions` should be a complete list of
site positions (in fractional coordinates) within the unit cell defined by
lattice vectors `lat_vecs`.
"""
function Crystal(lat_vecs, positions; types::Union{Nothing,Vector{String}}=nothing, symprec=1e-5)
    lat_vecs = convert(Mat3, lat_vecs)
    positions = [convert(Vec3, p) for p in positions]
    if isnothing(types)
        types = fill("", length(positions))
    end
    return crystal_from_inferred_symmetry(lat_vecs, positions, types; symprec)
end

"""
    Crystal(lat_vecs, positions, symbol::String; types=nothing, symprec=1e-5)

Build `Crystal` by applying symmetry operators for a given spacegroup symbol.
"""
function Crystal(lat_vecs::Mat3, positions, symbol::String; types::Union{Nothing,Vector{String}}=nothing, setting=nothing, symprec=1e-5)
    positions = [convert(Vec3, p) for p in positions]
    if isnothing(types)
        types = fill("", length(positions))
    end
    crystal_from_symbol(lat_vecs, positions, types, symbol; setting, symprec)
end

"""
    Crystal(lat_vecs, positions, spacegroup_number; types=nothing, symprec=1e-5)

Build `Crystal` by applying symmetry operators for a given international spacegroup
number.
"""
function Crystal(lat_vecs::Mat3, positions, spacegroup_number::Int; types::Union{Nothing,Vector{String}}=nothing, setting=nothing, symprec=1e-5)
    positions = [convert(Vec3, p) for p in positions]
    if isnothing(types)
        types = fill("", length(positions))
    end
    symbol = international_short_names[spacegroup_number]
    crystal_from_symbol(lat_vecs, positions, types, symbol; setting, symprec)
end


const n_hall_numbers = 530
const n_space_groups = 230

const international_short_names = begin
    ret = fill("", n_space_groups)
    for n in 1:n_hall_numbers
        sgt = Spglib.get_spacegroup_type(n)
        ret[sgt.number] = sgt.international_short
    end
    ret
end

function spacegroup_name(hall_number::Int)
    # String representation of space group
    sgt = Spglib.get_spacegroup_type(hall_number)
    return "HM symbol '$(sgt.international)' ($(sgt.number))"
end

function symops_from_spglib(rotations, translations)
    Rs = Mat3.(transpose.(eachslice(rotations, dims=3)))
    Ts = Vec3.(eachcol(translations))
    symops = SymOp.(Rs, Ts)
end


# Sort the sites according to class and fractional coordinates.
function sort_sites!(cryst::Crystal)
    function less_than(i, j)
        eci = cryst.classes[i]
        ecj = cryst.classes[j]
        if eci != ecj
            return eci < ecj
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


function crystal_from_inferred_symmetry(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec=1e-5)
    for i in 1:length(positions)
        for j in i+1:length(positions)
            ri = positions[i]
            rj = positions[j]
            if is_same_position(ri, rj; symprec)
                error("Positions $ri and $rj are symmetry equivalent.")
            end
        end
    end

    positions = wrap_to_unit_cell.(positions; symprec)

    cell = Spglib.Cell(lat_vecs, hcat(positions...), types)
    d = Spglib.get_dataset(cell, symprec)
    classes = d.equivalent_atoms .+ 1
    symops = symops_from_spglib(d.rotations, d.translations)
    spacegroup = spacegroup_name(d.hall_number)

    multiplicities = map(1:d.n_atoms) do i
        prim_idx = d.mapping_to_primitive[i]
        count(==(prim_idx), d.std_mapping_to_primitive)
    end
    sitesymmetries = SiteSymmetry.(d.site_symmetry_symbols, multiplicities, d.wyckoffs)

    ret = Crystal(lat_vecs, positions, classes, types, symops, spacegroup, sitesymmetries, symprec)
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

    if hall_cell == FastDipole.monoclinic
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

    for hall_number in 1:n_hall_numbers
        sgt = Spglib.get_spacegroup_type(hall_number)

        if (replace(symbol, " "=>"") == sgt.international_short || 
            symbol in [sgt.hall_symbol, sgt.international, sgt.international_full])

            # Some Hall numbers may be incompatible with unit cell of provided
            # lattice vectors; skip them.
            is_compatible = true

            cell = cell_type(lat_vecs)
            hall_cell = cell_type(hall_number)
            allowed_cells = all_compatible_cells(hall_cell)

            # Special handling of trigonal space groups
            if FastDipole.is_trigonal_symmetry(hall_number)
                # Trigonal symmetry must have either hexagonal or rhombohedral
                # cell, according to the Hall number.
                is_latvecs_valid = cell in [FastDipole.rhombohedral, FastDipole.hexagonal]
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

            if hall_cell == FastDipole.monoclinic
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

        println("The symbol '$symbol' is ambiguous! Returning all crystals:")
        for (i, (hall_number, c)) in enumerate(zip(hall_numbers, crysts))
            sgt = Spglib.get_spacegroup_type(hall_number)
            hm_symbol = sgt.international
            choice = sgt.choice
            n_atoms = length(c.positions)
            i_str = @sprintf "%2d" i
            natoms_str = @sprintf "%2d" n_atoms
            println("   $i_str. \"$hm_symbol\", setting=\"$choice\", generates $natoms_str atoms")
        end
        println()
        println("Note: To disambiguate, you may pass a named parameter, setting=\"...\".")
        println()
        return crysts
    end
end

"Build Crystal from explicit set of symmetry operations and a minimal set of positions "
function crystal_from_symops(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}, symops::Vector{SymOp}, spacegroup::String; symprec=1e-5)
    all_positions = Vec3[]
    all_types = String[]
    classes = Int[]
    
    for i = eachindex(positions)
        for s = symops
            x = wrap_to_unit_cell(transform(s, positions[i]); symprec)

            idx = findfirst(y -> is_same_position(x, y; symprec), all_positions)
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

    # Check that symops are present in Spglib-inferred space group
    cryst′ = crystal_from_inferred_symmetry(lat_vecs, all_positions, all_types; symprec)
    for s in symops
        is_found = any(cryst′.symops) do s′
            isapprox(s, s′; atol=symprec)
        end
        if !is_found
            println("WARNING: User provided symmetry operation could not be inferred by Spglib,")
            println("which likely indicates a non-conventional unit cell.")
            break
        end
    end

    # TODO: Extract pointgroups from spglib

    ret = Crystal(lat_vecs, all_positions, classes, all_types, symops, spacegroup, nothing, symprec)
    sort_sites!(ret)
    validate(ret)
    return ret
end

"""
    subcrystal(cryst, types) :: Crystal

Filter sublattices of a `Crystal` by atom `types`, keeping the space group
unchanged.
"""
function subcrystal(cryst::Crystal, types::Vararg{String, N}) where {N}
    for s in types
        if !(s in cryst.types)
            error("types string '$s' is not present in crystal.")
        end
    end
    idxs = findall(in(types), cryst.types)
    classes = unique(cryst.classes[idxs])
    return subcrystal(cryst, classes...)
end

"""
    subcrystal(cryst, classes) :: Crystal

Filter sublattices of `Crystal` by equivalence `classes`, keeping the space
group unchanged.
"""
function subcrystal(cryst::Crystal, classes::Vararg{Int, N}) where {N}
    for c in classes
        if !(c in cryst.classes)
            error("Class '$c' is not present in crystal.")
        end
    end

    idxs = findall(in(classes), cryst.classes)
    new_positions = cryst.positions[idxs]
    new_types = cryst.types[idxs]
    new_classes = cryst.classes[idxs]

    if idxs != 1:maximum(idxs)
        println("WARNING, atoms are being renumbered.")
    end

    ret = Crystal(cryst.lat_vecs, new_positions, new_classes, new_types,
                  cryst.symops, cryst.spacegroup, cryst.sitesymmetries, cryst.symprec)
    return ret
end

function Base.display(cryst::Crystal)
    printstyled("Crystal info\n"; bold=true, color=:underline)
    println(cryst.spacegroup)

    if is_standard_form(cryst.lat_vecs)
        (a, b, c, α, β, γ) = lattice_params(cryst.lat_vecs)
        @printf "Lattice params a=%.4g, b=%.4g, c=%.4g, α=%.4g°, β=%.4g°, γ=%.4g°\n" a b c α β γ
    else
        println("Lattice vectors:")
        for a in eachcol(cryst.lat_vecs)
            @printf "   [%.4g %.4g %.4g]\n" a[1] a[2] a[3]
        end
    end

    @printf "Volume %.4g\n" cell_volume(cryst)

    for c in unique(cryst.classes)
        i = findfirst(==(c), cryst.classes)
        print("Class $c")
        if cryst.types[i] != ""
            print(", type '$(cryst.types[i])'")
        end
        if !isnothing(cryst.sitesymmetries)
            # TODO: simplify in Julia 1.7
            symbol = cryst.sitesymmetries[c].symbol
            multiplicity = cryst.sitesymmetries[c].multiplicity
            wyckoff = cryst.sitesymmetries[c].wyckoff
            print(", Site sym. '$symbol', Wyckoff $multiplicity$wyckoff")
        end
        println(":")

        for i in findall(==(c), cryst.classes)
            r = cryst.positions[i]
            @printf "   %d. [%.4g, %.4g, %.4g]\n" i r[1] r[2] r[3]
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
