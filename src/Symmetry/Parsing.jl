# Functions for parsing Crystals / SymOps from text / .cif files

# Strips trailing uncertainty values from a String, then parses as a Float
function parse_cif_float(str::String) :: Float64
    i = findfirst('(', str)
    str = isnothing(i) ? str : str[1:i-1]
    return parse(Float64, str)
end


function parse_number_or_fraction(s)
    # Parse a number or fraction
    cnt = count(==('/'), s)
    if cnt == 0
        return parse(Float64, s)
    elseif cnt == 1
        n, d = parse.(Float64, split(s, '/'))
        return n/d
    else
        error("Cannot parse '$s' as a number.")
    end
end


function parse_op(str::AbstractString) :: SymOp
    D = 3
    R = zeros(D, D)
    T = zeros(D)

    @assert length(split(str, ",")) == D

    # The substring s describes the i'th row of the output
    for (i, s) in enumerate(split(str, ","))
        # Remove all spaces
        s = replace(s, " " => "")
        # In weird circumstances, double negatives may be generated. Get rid of them.
        s = replace(s, "--" => "+")
        # Trick to split on either + or - symbols
        s = replace(s, "-" => "+-")

        # Each term t describes the numerical value for one element of R or T
        for t in split(s, '+', keepempty=false)
            j = findfirst(last(t), "xyz")
            if isnothing(j)
                T[i] += parse_number_or_fraction(t)
            else
                coeff = t[begin:end-1]
                R[i, j] += begin
                    if coeff == ""
                        1
                    elseif coeff == "-"
                        -1
                    else
                        parse_number_or_fraction(coeff)
                    end
                end
            end
        end
    end

    # Wrap translations
    T = mod.(T, 1)

    return SymOp(R, T)
end


# Reads the crystal from a `.cif` file located at the path `filename`.
function Crystal(filename::AbstractString; symprec=nothing, override_symmetry=false)
    cif = CIF.Cif(filename)
    # For now, assumes there is only one data collection per .cif
    cif = cif[first(keys(cif))]
    oneof(fields...) = findfirstval(in(keys(cif)), fields)

    a = parse_cif_float(cif["_cell_length_a"][1])
    b = parse_cif_float(cif["_cell_length_b"][1])
    c = parse_cif_float(cif["_cell_length_c"][1])
    α = parse_cif_float(cif["_cell_angle_alpha"][1])
    β = parse_cif_float(cif["_cell_angle_beta"][1])
    γ = parse_cif_float(cif["_cell_angle_gamma"][1])
    latvecs = lattice_vectors(a, b, c, α, β, γ)

    geo_table = CIF.get_loop(cif, "_atom_site_fract_x")
    xs = parse_cif_float.(geo_table[:, "_atom_site_fract_x"])
    ys = parse_cif_float.(geo_table[:, "_atom_site_fract_y"])
    zs = parse_cif_float.(geo_table[:, "_atom_site_fract_z"])
    positions = Vec3.(zip(xs, ys, zs))

    types = nothing
    if "_atom_site_type_symbol" in keys(cif)
        types = String.(geo_table[:, "_atom_site_type_symbol"])
    end

    if "_atom_site_label" in keys(cif)
        labels = String.(geo_table[:, "_atom_site_label"])
    else
        labels = types
    end

    # Try to infer symprec from coordinate strings. TODO: Use uncertainty
    # information if available from .cif
    if isnothing(symprec)
        elems = vcat(xs, ys, zs)
        # guess fractional errors by assuming each elem is a fraction with simple denominator (2, 3, or 4)
        errs = map(elems) do x
            c = 12
            return abs(rem(x*c, 1, RoundNearest)) / c
        end
        (err, _) = findmax(errs)
        if err < 1e-12
            @info """Precision parameter is unspecified, but all coordinates seem to be simple fractions.
                     Setting symprec=1e-12."""
            symprec = 1e-12
        elseif 1e-12 < err < 1e-4
            symprec = 15err
            err_str = @sprintf "%.1e" err
            symprec_str = @sprintf "%.1e" symprec
            @info """Precision parameter is unspecified, but coordinate string '$s' seems to have error $err_str.
                     Setting symprec=$symprec_str."""
        else
            error("Cannot infer precision. Please provide an explicit `symprec` parameter to load '$filename'")
        end
    end

    symops = nothing
    sym_header = oneof("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
    if !isnothing(sym_header)
        sym_table = CIF.get_loop(cif, sym_header)
        symops = parse_op.(sym_table[:, sym_header])
    end

    sg_number = nothing
    sg_number_header = oneof("_space_group_it_number", "_symmetry_int_tables_number")
    if !isnothing(sg_number_header)
        sg_number = parse(Int, cif[sg_number_header][1])
    end

    hall_symbol = nothing
    hall_header = "_space_group_name_hall"
    if hall_header in keys(cif)
        hall_symbol = cif[hall_header][1]
    end

    hm_symbol = nothing
    hm_header = "_space_group_name_h-m_alt"
    if hm_header in keys(cif)
        hm_symbol = cif[hm_header][1]
    end

    # If we're overriding symmetry, then distinct labels may need to be merged
    # into equivalent sites. E.g., site-symmetry may be broken due to magnetic
    # order in an .mcif -- undo that.
    classes = override_symmetry ? types : labels

    # If chemical cell symops are missing, attempt to build a crystal from
    # magnetic symops
    if isnothing(symops)
        # IUCR standard or legacy field name, respectively
        mcif_fields = (;
            magn_operation_xyz=oneof("_space_group_symop_magn_operation.xyz", "_space_group_symop.magn_operation_xyz"),
            magn_centering_xyz=oneof("_space_group_symop_magn_centering.xyz", "_space_group_symop.magn_centering_xyz"),
        )
        if !isnothing(mcif_fields.magn_operation_xyz)
            sym_table = CIF.get_loop(cif, mcif_fields.magn_operation_xyz)
            operations = MSymOp.(sym_table[:, mcif_fields.magn_operation_xyz])
            sym_table = CIF.get_loop(cif, mcif_fields.magn_centering_xyz)
            centerings = MSymOp.(sym_table[:, mcif_fields.magn_centering_xyz])
            symops = vec([m1*m2 for m1 in operations, m2 in centerings])

            # Convert each MSymOp to SymOp, dropping parity field. The purpose
            # is to reconstruct all atom positions in crystal_from_symops.
            symops = [SymOp(s.R, s.T) for s in symops]

            # Fill atom positions by symmetry and infer symmetry operations
            orbits = crystallographic_orbits_distinct(symops, positions; symprec)
            all_positions = reduce(vcat, orbits)
            all_types = repeat_multiple(classes, length.(orbits))
            supercell = crystal_from_inferred_symmetry(latvecs, all_positions, all_types; symprec, check_cell=false)

            if override_symmetry
                # Don't idealize because later `set_dipoles_from_mcif!` will
                # need to search for atom positions using this specific setting.
                return standardize(supercell; idealize=false)
            else
                @warn """Loading the magnetic cell as chemical cell for TESTING PURPOSES only.
                         Set the option `override_symmetry=true` to infer the standard chemical
                         cell and its spacegroup symmetries."""
                return supercell
            end
        end
    end

    if !isnothing(symops)
        isnothing(sg_number) && error("Spacegroup number not specified.")
        hall_number = hall_number_from_symops(sg_number, symops)
        isnothing(hall_number) && error("Symmetry operations do not match an ITA setting for spacegroup $sg_number.")
        sg_label = spacegroup_label(hall_number)
        sg_setting = mapping_to_standard_setting(hall_number)
        sg = Spacegroup(symops, sg_label, sg_number, sg_setting)
        ret = crystal_from_spacegroup(latvecs, positions, classes, sg; symprec)
    elseif !isnothing(hall_symbol)
        # Use symmetries for Hall symbol
        ret = Crystal(latvecs, positions, hall_symbol; types=classes, symprec)
    else
        error("Spacegroup symops are unknown.")
    end

    return override_symmetry ? standardize(ret) : ret
end
