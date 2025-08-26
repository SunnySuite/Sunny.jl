# Reading CIF and mCIF files

function parse_cif_float_with_err(str::String; maybe_fractional)
    parts = split(str, '(')
    number_str = parts[1]
    number = parse(Float64, number_str)
    sigfigs = in('.', number_str) ? length(split(number_str, '.')[2]) : 0

    if length(parts) == 1
        # E.g., 0.2 has an assumed error bar of ±0.05
        err = 10.0^(-sigfigs) / 2

        # Components of each Wyckoff position in an ITA standard cell may be an
        # exact multiple of {1/2, 1/3, 1/4, 1/6, 1/8}. For example, the position
        # 0.125 = 1/8 is likely to be exact. Such strings do not contribute to
        # the estimated error when `maybe_fractional=true`.
        if maybe_fractional
            smallest_remainder = minimum(abs(rem(number * c, 1, RoundNearest)) / c for c in (2, 3, 4, 6, 8))
            if smallest_remainder < 1e-12
                err = 0.0
            end
        end
    else
        # E.g., 0.2(3) has an error bar of ±0.3
        @assert parts[2][end] == ')'
        err_str = parts[2][begin:end-1]
        err = 10.0^(-sigfigs) * parse(Int, err_str)
    end

    return (number, err)
end

function parse_cif_float(str::String)
    return parse_cif_float_with_err(str; maybe_fractional=false)[1]
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


# Cannot use crystallographic_orbit(::WyckoffExpr) because Wyckoffs and
# spacegroup setting data may not be available during CIF loading.
function crystallographic_orbit(position::Vec3; symops::Vector{SymOp}, symprec)
    orbit = Vec3[]
    for s = symops
        x = wrap_to_unit_cell(transform(s, position); symprec)
        if !any(y -> is_periodic_copy(x, y; symprec), orbit)
            push!(orbit, x)
        end
    end
    return orbit
end


# Reads the crystal from a CIF or mCIF located at the path `filename`. See
# extended doc string in Crystal.jl.
function Crystal(filename::AbstractString; keep_supercell=false, symprec=nothing, override_symmetry=nothing)
    if !isnothing(override_symmetry)
        @warn "Ignoring deprecated `override_symmetry` option"
    end

    cif = CIF.Cif(filename)
    # Accessor `oneof` assumes collections in CIF are unique
    cif = cif[first(keys(cif))]
    oneof(fields...) = findfirstval(in(keys(cif)), fields)

    # The user-provided symprec may be nothing, so we also infer a symprec from
    # the coordinate positions represented as strings.
    symprec0 = symprec
    symprec = 0.0

    (a, err_a) = parse_cif_float_with_err(cif["_cell_length_a"][1]; maybe_fractional=false)
    (b, err_b) = parse_cif_float_with_err(cif["_cell_length_b"][1]; maybe_fractional=false)
    (c, err_c) = parse_cif_float_with_err(cif["_cell_length_c"][1]; maybe_fractional=false)
    α = parse_cif_float(cif["_cell_angle_alpha"][1])
    β = parse_cif_float(cif["_cell_angle_beta"][1])
    γ = parse_cif_float(cif["_cell_angle_gamma"][1])
    latvecs = lattice_vectors(a, b, c, α, β, γ)

    geo_table = CIF.get_loop(cif, "_atom_site_fract_x")
    positions = Vec3[]
    for i in axes(geo_table, 1)
        (x, err_x) = parse_cif_float_with_err(geo_table[i, "_atom_site_fract_x"]; maybe_fractional=true)
        (y, err_y) = parse_cif_float_with_err(geo_table[i, "_atom_site_fract_y"]; maybe_fractional=true)
        (z, err_z) = parse_cif_float_with_err(geo_table[i, "_atom_site_fract_z"]; maybe_fractional=true)
        push!(positions, Vec3(x, y, z))
        symprec = max(symprec, err_x, err_y, err_z)
    end

    types = nothing
    if "_atom_site_type_symbol" in keys(cif)
        types = String.(geo_table[:, "_atom_site_type_symbol"])
    end

    if "_atom_site_label" in keys(cif)
        labels = String.(geo_table[:, "_atom_site_label"])
    else
        labels = types
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

    #=
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
    =#

    classes = labels

    # If chemical cell symops are missing, attempt to build a crystal from
    # magnetic symops
    from_mcif = false
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

            # Some magnetically inequivalent sites may become equivalent in the
            # chemical cell. This analysis requires type labels for the atoms.
            isnothing(types) && error("Missing _atom_site_type_symbol data")
            classes = types

            # If the supercell dimensions are not all equal, then it is possible
            # that the spacegroup symmetry requires them to have a precise
            # ratio. Any uncertainty in (a, b, c) should therefore be
            # incorporated into symprec. See test/cifs/ZnFe2O4.cif for a
            # real-world example.
            if !allequal((a, b, c))
                symprec = max(symprec, err_a/a, err_b/b, err_c/c)
            end

            # Remember that symops were projected from msymops
            from_mcif = true
        end
    end

    if isnothing(symops)
        error("Missing spacegroup symmetry operations")
    end

    # If user-provided symprec0 is unavailable, use inferred symprec. The 10x
    # fudge factor and 1e-8 lower bound give some safety margin.
    symprec = @something symprec0 max(10*symprec, 1e-8)

    # Fill atom positions by symmetry and infer symmetry operations
    orbits = crystallographic_orbit.(positions; symops, symprec)
    validate_orbits(positions, orbits; symprec)
    all_positions = reduce(vcat, orbits)
    all_types = repeat_multiple(classes, length.(orbits))

    # TODO: Check orbit lengths against _atom_site_symmetry_multiplicity data if
    # available.

    # Infer symmetry from full list of atoms. Disable cell check because we will
    # customize the warning message below.
    ret = crystal_from_inferred_symmetry(latvecs, all_positions, all_types; symprec, check_cell=false)
    sort_sites!(ret)

    if from_mcif
        if !keep_supercell
            # Disable idealization when loading an mCIF because a subsequent
            # call to `set_dipoles_from_mcif!` must be able to search for atom
            # positions using the mCIF setting.
            return standardize(ret; idealize=false)
        else
            @warn "Use `keep_supercell=true` for testing purposes only! Inferred symmetries \
                   are unreliable."
            return ret
        end
    else
        symprec_str = number_to_simple_string(symprec; digits=2)
        suggest_symprec = isnothing(symprec0) ? " Try overriding `symprec` (inferred $symprec_str)." : ""

        if !isnothing(sg_number) && sg_number != ret.sg.number
            error("Inferred spacegroup $(ret.sg.number) differs from $sg_number in CIF.$suggest_symprec")
        end

        symops_consistent = isapprox(symops, ret.sg.symops; atol=symprec)
        hall_number_inferred = hall_number_from_symops(ret.sg.number, ret.sg.symops; atol=symprec)

        if symops_consistent
            if !isnothing(hall_number_inferred)
                # The inferred ret.sg data (setting and symops) may be slightly
                # perturbed from the ITA standard tables. This happens for
                # test/cifs/UPt3.cif, with spacegroup 194. Wyckoff 6h has a site at
                # position (x,2x,1/4). The CIF stores x and 2x using strings with
                # truncated decimal expansions, 0.333 and 0.667. The ratio of these
                # numbers slightly deviates from 2, causing Spglib to infer a
                # slightly perturbed setting. Get clean symmetry data using tables
                # for the inferred Hall number.
                sg = Spacegroup(hall_number_inferred)
                ret = crystal_from_spacegroup(latvecs, positions, classes, sg; symprec)
            else
                @info "Cell appears non-standard. Consider `standardize(cryst)` and then \
                       `reshape_supercell(sys, shape)` for calculations on an arbitrarily \
                       shaped system."
            end
        else
            @warn "Inconsistent symmetry operations! This may occur with an incomplete CIF, \
                   a non-standard setting, or failed inference.$suggest_symprec"
        end
        return ret
    end
end
