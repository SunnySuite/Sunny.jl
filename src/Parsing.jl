using CrystalInfoFramework
using FilePaths
using Printf

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

function _parse_op(s::AbstractString) :: Symmetry.SymOp
    R, T = xyzt2components(s, Val(3))
    # Wrap the translations to be within range [0, 1)
    T = mod.(T, 1)
    Symmetry.SymOp(R, T)
end

# These two functions are from Crystalline.jl
# TODO: Add licensing info before public release
function xyzt2components(s::AbstractString, ::Val{D}) where D
    xyzts = split(s, ',')
    length(xyzts) == D || throw(DimensionMismatch("incompatible matrix size and string format"))

    # initialize zero'd MArrays for rotation/translation parts (allocation will be elided)
    W = zero(MMatrix{D, D, Float64}) # rotation
    w = zero(MVector{D, Float64})    # translation
    
    # "fill in" `W` and `w` according to content of `xyzts`
    xyzt2components!(W, w, xyzts)

    # convert to SArrays (elides allocation since `xyzt2components!` is inlined)
    return SMatrix(W), SVector(w)
end

@inline function xyzt2components!(W::MMatrix{D,D,T}, w::MVector{D,T},
                                  xyzts::AbstractVector{<:AbstractString}) where {D,T<:Real}

    chars = D == 3 ? ('x','y','z') : D == 2 ? ('x','y') : ('x',)
    @inbounds for (i,s) in enumerate(xyzts)
        # rotation/inversion/reflection part
        firstidx = nextidx = firstindex(s)
        while (idx = findnext(c -> c ∈ chars, s, nextidx)) !== nothing
            c = s[idx]
            j = c=='x' ? 1 : (c=='y' ? 2 : 3)
            
            if idx == firstidx
                W[i,j] = one(T)
            else
                previdx = prevind(s, idx)
                while (c′=s[previdx]; isspace(s[previdx]))
                    previdx = prevind(s, previdx)
                    previdx ≤ firstidx && break
                end
                if c′ == '+' || isspace(c′)
                    W[i,j] = one(T)
                elseif c′ == '-'
                    W[i,j] = -one(T)
                else
                    throw(ArgumentError("failed to parse provided string representation"))
                end
            end
            nextidx = nextind(s, idx)
        end

        # nonsymmorphic part/fractional translation part
        lastidx = lastindex(s)
        if nextidx ≤ lastidx # ⇒ stuff "remaining" in `s`; a nonsymmorphic part
            slashidx = findnext(==('/'), s, nextidx)
            if slashidx !== nothing # interpret as integer fraction
                num = SubString(s, nextidx, prevind(s, slashidx))
                den = SubString(s, nextind(s, slashidx), lastidx)
                w[i] = convert(T, parse(Int, num)/parse(Int, den))
            else                    # interpret at number of type `T`
                w[i] = parse(T, SubString(s, nextidx, lastidx))
            end
        end
    end

    return W, w
end

"""
    Crystal(filename::AbstractString; symprec=1e-5)

Parse a `Crystal` from a `.cif` file located at the path of `filename`.
"""
function Crystal(filename::AbstractString; symprec=nothing)
    cif = Cif(Path(filename))
    # For now, assumes there is only one data collection per .cif
    cif = cif[first(keys(cif))]

    a = _parse_cif_float(cif["_cell_length_a"][1])
    b = _parse_cif_float(cif["_cell_length_b"][1])
    c = _parse_cif_float(cif["_cell_length_c"][1])
    α = _parse_cif_float(cif["_cell_angle_alpha"][1])
    β = _parse_cif_float(cif["_cell_angle_beta"][1])
    γ = _parse_cif_float(cif["_cell_angle_gamma"][1])
    lat_vecs = lattice_vectors(a, b, c, α, β, γ)

    geo_table = get_loop(cif, "_atom_site_fract_x")
    xs = map(_parse_cif_float, geo_table[!, "_atom_site_fract_x"])
    ys = map(_parse_cif_float, geo_table[!, "_atom_site_fract_y"])
    zs = map(_parse_cif_float, geo_table[!, "_atom_site_fract_z"])
    sitetypes = geo_table[!, "_atom_site_type_symbol"]
    sitetypes = Vector{String}(sitetypes)
    sitelabels = geo_table[!, "_atom_site_label"]
    unique_atoms = Vec3.(zip(xs, ys, zs))

    # Try to infer symprec from coordinate strings
    # TODO: Use uncertainty information if available from .cif
    if isnothing(symprec)
        strs = vcat(geo_table[!, "_atom_site_fract_x"], geo_table[!, "_atom_site_fract_y"], geo_table[!, "_atom_site_fract_z"])
        elems = vcat(xs, ys, zs)
        # guess fractional errors by assuming each elem is a fraction with simple denominator (2, 3, or 4)
        errs = map(elems) do x
            c = 12
            return abs(rem(x*c, 1, RoundNearest)) / c
        end
        (err, i) = findmax(errs)
        if err < 1e-12
            println("Warning: Precision parameter is unspecified.")
            println("All coordinate strings seem to be simple fractions. Setting symprec=1e-12.")
            symprec = 1e-12
        elseif 1e-12 < err < 1e-4
            println("Warning: Precision parameter is unspecified.")
            symprec = 15err
            @printf "Inferring that coordinate string '%s' has error %.1e. Setting symprec=%.1e.\n" strs[i] err symprec
        else
            println("Error: Please specify an explicit `symprec` parameter to load this file, '$filename'")
            return Nothing
        end
    end

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
        # Use explicitly provided symmetries
        return Crystal(lat_vecs, unique_atoms, sitetypes, symmetries; symprec)
    elseif !isnothing(hall_symbol)
        # Read symmetries from database for Hall symbol
        return Crystal(lat_vecs, unique_atoms, sitetypes, hall_symbol; symprec)
    else
        # Infer the symmetries automatically
        return Crystal(lat_vecs, unique_atoms, sitetypes; symprec)
    end
end
