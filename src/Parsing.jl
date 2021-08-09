using CrystalInfoFramework
using FilePaths

struct ValidateError <: Exception end
Base.showerror(io::IO, e::ValidateError) = print(io, e)

function _parse_lattice(config::Dict{String, Any}) :: Lattice
    try
        dim = config["dimension"]
        lat_vecs = config["lattice_vectors"]
        basis_vecs = config["basis_vectors"]
        lattice_size = config["lattice_size"]

        @assert length(lat_vecs) == dim
        @assert all(v -> length(v) == dim, lat_vecs)
        @assert all(v -> length(v) == dim, basis_vecs)
        @assert length(lattice_size) == dim

        lat_vecs = SMatrix{dim, dim, Float64, dim*dim}(hcat(lat_vecs...))
        basis_vecs = [SVector{dim}(bv) for bv in basis_vecs]
        # basis_vecs = hcat(basis_vecs...
        lattice_size = SVector{dim}(lattice_size)

        return Lattice{dim, dim*dim, dim+1}(lat_vecs, basis_vecs, lattice_size) 
    catch err
        if isa(err, KeyError)
            throw(ValidateError("lattice config missing mandatory key: $(err.key)"))
        else
            # Re-raise all other kinds of errors
            throw(err)
        end
    end
end

function _parse_interactions(config::Dict{String, Any}, lattice::Lattice) :: Vector{Interaction}
    interactions = Vector{Interaction}()
    crystal = Crystal(lattice)
    if haskey(config, "field")
        for field_int in config["field"]
            push!(interactions, ExternalField(field_int["strength"]))
        end
    end

    if haskey(config, "pair")
        for pair_int in config["pair"]
            push!(interactions, PairInteraction(
                pair_int["strength"],
                crystal,
                pair_int["dist"],
                pair_int["class"],
            ))
        end
    end

    if haskey(config, "onsite")
        for easy_ax in config["onsite"]
            push!(interactions, OnSite(easy_ax["J"]))
        end
    end
    return interactions
end

function parse_config(filename::String) :: SpinSystem
    try
        config = TOML.tryparsefile(filename)
        lattice = _parse_lattice(config["lattice"])
        interactions = _parse_interactions(config["model"], lattice)

        return SpinSystem(lattice, interactions)
    catch err
        if isa(err, TOML.ParserError)
            println("Parse error on line $(err.line), column $(err.column) of $(filename).")
        end
    end
end

"Strips trailing uncertainty values from a String, then parses as a Float"
function _parse_cif_float(str::String) :: Float64
    i = findfirst('(', str)
    str = isnothing(i) ? str : str[1:i-1]
    return parse(Float64, str)
end

"Parses a symmetry transformation from a String"
function _parse_op(str::String) :: Symmetry.SymOp
    res = zeros(Float64, 3, 4)
    R = zeros(3, 3)
    T = zeros(3)
    for (i, str_row) = enumerate(map(strip, split(str, ",")))
        (sign, coord, rest) = match(r"(\-?)([xyz])(.*)", str_row).captures
        (fsign, numer, denom) = isempty(rest) ? ("+", "0", "1") : match(r"([\+\-]+)(\N)/(\N)", rest).captures
        sign = (sign == "-" ? -1 : +1)
        fsign = (fsign == "-" ? -1 : +1)
        coord = Dict("x" => 1, "y" => 2, "z" => 3)[coord]
        (numer, denom) = parse.(Float64, [numer, denom])
        R[i, coord] = sign
        T[i] = fsign * numer/denom
    end
    return Symmetry.SymOp(R, T)
end

function Crystal(filename::AbstractPath)
    cif = Cif(filename)
    # For now, assumes there is only one data collection per .cif
    cif = cif[first(keys(cif))]

    a = _parse_cif_float(cif["_cell_length_a"][1])
    b = _parse_cif_float(cif["_cell_length_b"][1])
    c = _parse_cif_float(cif["_cell_length_c"][1])
    α = _parse_cif_float(cif["_cell_angle_alpha"][1])
    β = _parse_cif_float(cif["_cell_angle_beta"][1])
    γ = _parse_cif_float(cif["_cell_angle_gamma"][1])
    bravais_lat = Lattice(a, b, c, α, β, γ, [1, 1, 1])

    geo_table = get_loop(cif, "_atom_site_fract_x")
    xs = map(_parse_cif_float, geo_table[!, "_atom_site_fract_x"])
    ys = map(_parse_cif_float, geo_table[!, "_atom_site_fract_y"])
    zs = map(_parse_cif_float, geo_table[!, "_atom_site_fract_z"])
    sitetypes = geo_table[!, "_atom_site_type_symbol"]
    sitetypes = Vector{String}(sitetypes)
    sitelabels = geo_table[!, "_atom_site_label"]
    unique_atoms = Vec3.(zip(xs, ys, zs))

    # By default, assume P1 symmetry
    symmetries = [ Symmetry.SymOp(I(3), @SVector zeros(3)) ]
    # There's two possible headers symmetry operations could be listed under
    for sym_header in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
        if sym_header in keys(cif)
            sym_table = get_loop(cif, sym_header)
            symmetries = map(_parse_op, sym_table[!, sym_header])
        end
    end

    # By default, assume P1 symmetry
    groupnum = 1
    # Two possible headers for symmetry group number
    for group_header in ("_space_group_it_number", "_symmetry_int_tables_number")
        if group_header in keys(cif)
            groupnum = parse(Int, cif[group_header][1])
        end
    end

    hall_symbol = nothing
    hall_header = "_space_group_name_Hall"
    if hall_header in keys(cif)
        hall_symbol = cif[hall_header]
    end

    # Symmetry preferences: Explicit List > Hall Symbol > Infer
    if length(symmetries) > 1
        # Assume all symmetries have been provided
        return Crystal(bravais_lat.lat_vecs, unique_atoms, sitetypes, symmetries)
    elseif !isnothing(hall_symbol)
        return Crystal(bravais_lat.lat_vecs, unique_atoms, sitetypes, hall_symbol)
    else
        # Infer the symmetries automatically
        return Crystal(bravais_lat.lat_vecs, unique_atoms, sitetypes)
    end
end
