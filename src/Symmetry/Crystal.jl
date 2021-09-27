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
    isapprox(s1.R, s2.R; atol) && isapprox(s1.T, s2.T; atol)
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
    types                :: Vector{String}      # Species for each atom
    symops               :: Vector{SymOp}       # Symmetry operations
    hall_number          :: Union{Int, Nothing} # Hall number
    symprec              :: Float64             # Tolerance to imperfections in symmetry
end

nbasis(cryst::Crystal) = length(cryst.positions)
cell_volume(cryst::Crystal) = abs(det(lat.lat_vecs))
lattice_params(cryst::Crystal) = lattice_params(cryst.lat_vecs)
lattice_vectors(cryst::Crystal) = cryst.lat_vecs

"""
    Crystal(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec=1e-5)

Construct a `Crystal` using explicit geometry information, with all symmetry information
automatically inferred. `positions` should be a list of site positions (in fractional
coordinates) within the unit cell defined by lattice vectors which are the columns of `lat_vecs`.
"""
function Crystal(lat_vecs::Mat3, positions::Vector{Vec3}, types::Vector{String}; symprec=1e-5)
    positions = wrap_to_unit_cell.(positions; symprec)

    cell = Spglib.Cell(lat_vecs, hcat(positions...), types)
    d = Spglib.get_dataset(cell, symprec)
    equiv_atoms = d.equivalent_atoms
    Rs = Mat3.(transpose.(eachslice(d.rotations, dims=3)))
    Ts = Vec3.(eachcol(d.translations))
    symops = map(SymOp, Rs, Ts)
    
    # Sort atoms so that they are contiguous in equivalence classes
    p = sortperm(equiv_atoms)
    positions .= positions[p]
    types .= types[p]

    # Base constructor
    ret = Crystal(lat_vecs, positions, equiv_atoms, types, symops, d.hall_number, symprec)
    validate(ret)
    return ret
end

# Build Crystal using the space group denoted by a unique Hall number. The complete
# list is given at http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
# TODO: Make hall_number a named parameter. But then we have to avoid a conflict
# with the constructor above. Maybe use unique names for every constructor?
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, hall_number::Int; symprec=1e-5)
    cell = cell_type(lat_vecs)
    hall_cell = cell_type(hall_number)
    allowed_cells = all_compatible_cells(hall_cell)
    @assert cell in allowed_cells "Hall number $hall_number requires a $hall_cell cell, but found $cell."

    if hall_cell == monoclinic
        is_compatible = is_compatible_monoclinic_cell(lat_vecs, hall_number)
        @assert is_compatible "Lattice vectors define a monoclinic cell that is incompatible with Hall number $hall_number."
    end

    rotations, translations = Spglib.get_symmetry_from_database(hall_number)
    Rs = Mat3.(transpose.(eachslice(rotations, dims=3)))
    Ts = Vec3.(eachcol(translations))
    symops = map(SymOp, Rs, Ts)

    # Use symops to fill in symmetry-related atom positions
    return Crystal(lat_vecs, base_positions, base_types, symops; hall_number, symprec)
end

# Make best effort to build Crystal from symbolic representation of spacegroup
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, symbol::String; symprec=1e-5)
    # See "Complete list of space groups" at Seto's home page:
    # http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
    n_space_groups = 530

    crysts = Crystal[]
    for hall_number in 1:n_space_groups
        sgt = Spglib.get_spacegroup_type(hall_number)

        if (replace(symbol, " "=>"") == sgt.international_short || 
            symbol in [sgt.hall_symbol, sgt.international, sgt.international_full, "I:"*string(sgt.number)])

            # Some Hall numbers may be incompatible with unit cell of provided
            # lattice vectors; skip them.
            is_compatible = true

            cell = cell_type(lat_vecs)
            hall_cell = cell_type(hall_number)
            allowed_cells = all_compatible_cells(hall_cell)

            # Special handling of trigonal space groups
            if is_trigonal_symmetry(hall_number)
                # Trigonal symmetry must have either hexagonal or rhombohedral
                # cell, according to the Hall number.
                is_latvecs_valid = cell in [rhombohedral, hexagonal]
                @assert is_latvecs_valid "Symbol $symbol requires a rhomobohedral or hexagonal cell, but found $cell."
                is_compatible = cell in allowed_cells
            else
                # For all other symmetry types, there is a unique cell for each Hall number
                @assert cell in allowed_cells "Symbol $symbol requires a $hall_cell cell, but found $cell."
            end

            if hall_cell == monoclinic
                is_compatible = is_compatible_monoclinic_cell(lat_vecs, hall_number)
            end

            if is_compatible
                # Build crystal using symops for hall_number from database
                c = Crystal(lat_vecs, base_positions, base_types, hall_number; symprec)
                push!(crysts, c)
            end
        end
    end

    if length(crysts) == 0
        error("Could not find symbol '$symbol' in database.")
    elseif length(crysts) == 1
        return first(crysts)
    else
        sort!(crysts, by=c->length(c.positions))

        println("Warning, the symbol '$symbol' is ambiguous. It could refer to:")
        for c in crysts
            hall_number = c.hall_number
            hm_symbol = Spglib.get_spacegroup_type(hall_number).international
            n_atoms = length(c.positions)
            println("   HM symbol '$hm_symbol' (Hall number $hall_number), which generates $n_atoms atoms")
        end
        println()
        println("Selecting Hall number $(first(crysts).hall_number). You may wish to specify")
        println("an alternative Hall number in place of the symbol '$symbol'.")
        return first(crysts)
    end
end

"Build Crystal from explicit set of symmetry operations and a minimal set of positions "
function Crystal(lat_vecs::Mat3, base_positions::Vector{Vec3}, base_types::Vector{String}, symops::Vector{SymOp}; hall_number=nothing, symprec=1e-5)
    # Spglib can map a Hall number to symops (via `get_symmetry_from_database`)
    # and can map symops to a Hall number (via `get_hall_number_from_symmetry`).
    # Unfortunately the round trip is not the identity function:
    #   https://github.com/spglib/spglib/issues/132
    # As a workaround, allow the caller to pass an explicit `hall_number`. If
    # `nothing`, then Spglib can infer the Hall number.
    if isnothing(hall_number)
        rotation = zeros(3, 3, length(symops))
        translation = zeros(3, length(symops))
        for (i, s) = enumerate(symops)
            rotation[:, :, i] = s.R'
            translation[:, i] = s.T
        end
        # If rotation matrices are not all integral-valued, we are working with
        #  a nonstandard unit cell and can't infer the Hall number.
        if all(trunc.(rotation) .== rotation)
            hall_number = Int(Spglib.get_hall_number_from_symmetry(
                rotation, translation, length(symops)
            ))
        else
            @warn "Provided symmetry ops are for a nonstandard unit cell."
        end
    end
    
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
    if !isnothing(hall_number)
        cryst′ = Crystal(lat_vecs, positions, types; symprec)
        for s in symops
            @assert any(cryst′.symops) do s′
                isapprox(s, s′; atol=symprec)
            end
        end
    end

    # Call base constructor
    ret = Crystal(lat_vecs, positions, equiv_atoms, types, symops, hall_number, symprec)
    validate(ret)
    return ret
end

"""
    subcrystal(cryst, types) :: Crystal

Filter sublattices of a `Crystal` by types, keeping the symmetry group of the original `Crystal`.
"""
function subcrystal(cryst::Crystal, types::String) :: Crystal
    if !in(types, cryst.types)
        error("Species string '$types' is not present in crystal.")
    end
    subindexes = findall(isequal(types), cryst.types)
    new_positions = cryst.positions[subindexes]
    new_equiv_atoms = cryst.equiv_atoms[subindexes]
    # Reduce all of the equivalent atom indexes to count again from 1, 2,...
    unique_equiv = unique(new_equiv_atoms)
    for (i, unique) in enumerate(unique_equiv)
        equiv_sites = findall(isequal(unique), new_equiv_atoms)
        new_equiv_atoms[equiv_sites] .= i
    end
    new_types = fill(types, length(subindexes))

    return Crystal(cryst.lat_vecs, new_positions, new_equiv_atoms,
                   new_types, cryst.symops, cryst.hall_number, cryst.symprec)
end

"""
    subcrystal(cryst, equiv_idxs) :: Crystal

Filter sublattices of a `Crystal` by a list of indexes into `cryst.equiv_atoms`,
 keeping the symmetry group of the original `Crystal`.
"""
function subcrystal(cryst::Crystal, equiv_idxs::Vector{Int}) :: Crystal
    for equiv_idx in equiv_idxs
        if !in(equiv_idx, cryst.equiv_atoms)
            error("Equivalent index '$equiv_idx' is not present in crystal.")
        end
    end
    new_positions = empty(cryst.positions)
    new_types = empty(cryst.types)
    new_equiv_atoms = empty(cryst.equiv_atoms)

    for (i, equiv_idx) in enumerate(equiv_idxs)
        subindexes = findall(isequal(equiv_idx), cryst.equiv_atoms)
        append!(new_positions, cryst.positions[subindexes])
        append!(new_types, cryst.types[subindexes])
        append!(new_equiv_atoms, fill(i, length(subindexes)))
    end

    return Crystal(cryst.lat_vecs, new_positions, new_equiv_atoms,
                   new_types, cryst.symops, cryst.hall_number, cryst.symprec)
end

subcrystal(cryst::Crystal, equiv_idxs::Vararg{Int}) = subcrystal(cryst, equiv_idxs)

function Base.display(cryst::Crystal)
    printstyled("Crystal info\n"; bold=true, color=:underline)
    if isnothing(cryst.hall_number)
        println("Nonstandard unit cell, unknown space group")
    else
        sgt = Spglib.get_spacegroup_type(cryst.hall_number)
        # println("Hall group '$(sgt.hall_symbol)' (Hall number $(cryst.hall_number))")
        println("H-M space group '$(sgt.international)' ($(sgt.number))")
    end

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
            print("Species '$(cryst.types[i])', ")
        end
        if length(unique(cryst.equiv_atoms)) > length(unique(cryst.types))
            print("Class $(cryst.equiv_atoms[i]), ")
        end
        r = cryst.positions[i]
        @printf "Coords [%.4g, %.4g, %.4g]\n" r[1] r[2] r[3]
    end
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
