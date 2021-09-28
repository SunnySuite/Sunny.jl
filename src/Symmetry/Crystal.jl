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

"""
    Crystal

A type holding all geometry and symmetry information needed to represent
 a three-dimensional crystal.
"""
struct Crystal
    lat_vecs             :: Mat3                # Lattice vectors as columns
    positions            :: Vector{Vec3}        # Full set of atoms, fractional coords
    equiv_atoms          :: Vector{Int}         # Index to equivalent atom type
    types                :: Vector{String}      # Types for each atom
    symops               :: Vector{SymOp}       # Symmetry operations
    spacegroup           :: String              # Description of space group
    symprec              :: Float64             # Tolerance to imperfections in symmetry
end

nbasis(cryst::Crystal) = length(cryst.positions)
cell_volume(cryst::Crystal) = abs(det(lat.lat_vecs))
lattice_params(cryst::Crystal) = lattice_params(cryst.lat_vecs)
lattice_vectors(cryst::Crystal) = cryst.lat_vecs


function Crystal(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec=1e-5)
    return crystal_from_inferred_symmetry(lat_vecs, positions, types; symprec)
end

function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, symbol::String; symprec=1e-5)
    crystal_from_symbol(lat_vecs, base_positions, base_types, symbol; symprec)
end

function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, spacegroup_number::Int; symprec=1e-5)
    symbol = international_short_names[spacegroup_number]
    crystal_from_symbol(lat_vecs, base_positions, base_types, symbol; symprec)
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
    symops = map(SymOp, Rs, Ts)
end


"""
    crystal_from_inferred_symmetry(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec=1e-5)

Construct a `Crystal` using explicit geometry information, with all symmetry information
automatically inferred. `positions` should be a list of site positions (in fractional
coordinates) within the unit cell defined by lattice vectors which are the columns of `lat_vecs`.
"""
function crystal_from_inferred_symmetry(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec=1e-5)
    positions = wrap_to_unit_cell.(positions; symprec)

    cell = Spglib.Cell(lat_vecs, hcat(positions...), types)
    d = Spglib.get_dataset(cell, symprec)
    equiv_atoms = d.equivalent_atoms
    symops = symops_from_spglib(d.rotations, d.translations)
    spacegroup = spacegroup_name(d.hall_number)
    
    # Sort atoms so that they are contiguous in equivalence classes
    p = sortperm(equiv_atoms)
    positions .= positions[p]
    types .= types[p]

    # Base constructor
    ret = Crystal(lat_vecs, positions, equiv_atoms, types, symops, spacegroup, symprec)
    validate(ret)
    return ret
end


# Build Crystal using the space group denoted by a unique Hall number. The complete
# list is given at http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
function crystal_from_hall_number(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, hall_number::Int; symprec=1e-5)
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

    return crystal_from_symops(lat_vecs, base_positions, base_types, symops, spacegroup; symprec)
end

# Make best effort to build Crystal from symbolic representation of spacegroup
function crystal_from_symbol(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, symbol::String; symprec=1e-5)
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
                cryst = crystal_from_hall_number(lat_vecs, base_positions, base_types, hall_number; symprec)
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
        println("The symbol '$symbol' is ambiguous! Returning all crystals:")
        for (i, (hall_number, c)) in enumerate(zip(hall_numbers, crysts))
            sgt = Spglib.get_spacegroup_type(hall_number)
            hm_symbol = sgt.international
            choice = sgt.choice
            n_atoms = length(c.positions)
            i_str = @sprintf "%2d" i
            hall_str = @sprintf "%3d" hall_number
            natoms_str = @sprintf "%2d" n_atoms
            println("   $i_str. Hall number $hall_str generates $natoms_str atoms ('$hm_symbol' setting '$choice')")
        end
        println()
        println("Note: The constructor `crystal_from_hall_number()` would be unambiguous.")
        println()
        return crysts
    end
end

"Build Crystal from explicit set of symmetry operations and a minimal set of positions "
function crystal_from_symops(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, symops::Vector{SymOp}, spacegroup::String; symprec=1e-5)
    positions = Vec3[]
    types = String[]
    equiv_atoms = Int[]
    
    for i = eachindex(base_positions)
        for s = symops
            x = wrap_to_unit_cell(transform(s, base_positions[i]); symprec)

            idx = findfirst(y -> is_same_position(x, y; symprec), positions)
            if isnothing(idx)
                push!(positions, x)
                push!(types, base_types[i])
                push!(equiv_atoms, i)
            else
                j = equiv_atoms[idx]
                if i != j
                    error("Base positions $(base_positions[i]) and $(base_positions[j]) are symmetry equivalent.")
                end
            end
        end
    end

    # Check that symops are present in Spglib-inferred space group
    cryst′ = crystal_from_inferred_symmetry(lat_vecs, positions, types; symprec)
    for s in symops
        is_found = any(cryst′.symops) do s′
            isapprox(s, s′; atol=symprec)
        end
        if !is_found
            println("Warning: User provided symmetry operation could not be inferred by Spglib,")
            println("which likely indicates a non-conventional unit cell.")
            break
        end
    end

    ret = Crystal(lat_vecs, positions, equiv_atoms, types, symops, spacegroup, symprec)
    validate(ret)
    return ret
end

"""
    subcrystal(cryst, types) :: Crystal

Filter sublattices of a `Crystal` by types, keeping the symmetry group of the original `Crystal`.
"""
function subcrystal(cryst::Crystal, types::Vararg{String, N}) where {N}
    for s in types
        if !(s in cryst.types)
            error("types string '$s' is not present in crystal.")
        end
    end
    indxs = findall(in(types), cryst.types)
    equiv_idxs = unique(cryst.equiv_atoms[indxs])
    return subcrystal(cryst, equiv_idxs...)
end

"""
    subcrystal(cryst, equiv_idxs) :: Crystal

Filter sublattices of a `Crystal` by a list of indexes into `cryst.equiv_atoms`,
 keeping the symmetry group of the original `Crystal`.
"""
function subcrystal(cryst::Crystal, equiv_idxs::Vararg{Int, N}) where {N}
    for equiv_idx in equiv_idxs
        if !(equiv_idx in cryst.equiv_atoms)
            error("Equivalent index '$equiv_idx' is not present in crystal.")
        end
    end

    # Keep only atoms matching given "equivalent indices"
    new_positions = empty(cryst.positions)
    new_types = empty(cryst.types)
    new_equiv_atoms = empty(cryst.equiv_atoms)

    for (i, equiv_idx) in enumerate(equiv_idxs)
        idxs = findall(isequal(equiv_idx), cryst.equiv_atoms)
        append!(new_positions, cryst.positions[idxs])
        append!(new_types, cryst.types[idxs])
        append!(new_equiv_atoms, fill(i, length(idxs)))
    end

    return Crystal(cryst.lat_vecs, new_positions, new_equiv_atoms,
                   new_types, cryst.symops, cryst.spacegroup, cryst.symprec)
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

    println("Atoms:")
    for i in eachindex(cryst.positions)
        print("   $i. ")
        if length(unique(cryst.types)) > 1
            print("Types '$(cryst.types[i])', ")
        end
        if length(unique(cryst.equiv_atoms)) > length(unique(cryst.types))
            print("Class $(cryst.equiv_atoms[i]), ")
        end
        r = cryst.positions[i]
        @printf "Coords [%.4g, %.4g, %.4g]\n" r[1] r[2] r[3]
    end
    println()
end

function validate(cryst::Crystal)
    # Atoms must be sorted by equivalence class
    sortperm(cryst.equiv_atoms) == eachindex(cryst.equiv_atoms)

    # Equivalent atoms must have the same types
    for i in eachindex(cryst.positions)
        for j in eachindex(cryst.positions)
            if cryst.equiv_atoms[i] == cryst.equiv_atoms[j]
                @assert cryst.types[i] == cryst.types[j]
            end
        end
    end

    # Check symmetry rotations are orthogonal, up to periodicity
    for s in cryst.symops
        R = cryst.lat_vecs * s.R * inv(cryst.lat_vecs)
        # Due to possible imperfections in the lattice vectors, only require
        # that R is approximately orthogonal
        @assert norm(R*R' - I) < cryst.symprec "Lattice vectors and symmetry operations are incompatible."
    end

    # TODO: Check that space group is closed and that symops have inverse?
end
