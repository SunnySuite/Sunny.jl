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
    cif = CIF.Cif(Path(filename))
    # For now, assumes there is only one data collection per .cif
    cif = cif[first(keys(cif))]

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

    multiplicities = nothing
    if "_atom_site_symmetry_multiplicity" in names(geo_table)
        multiplicities = parse.(Int, geo_table[:, "_atom_site_symmetry_multiplicity"])
    end

    wyckoffs = nothing
    if "_atom_site_wyckoff_symbol" in names(geo_table)
        wyckoffs = String.(geo_table[:, "_atom_site_wyckoff_symbol"])
    end

    # Try to infer symprec from coordinate strings
    # TODO: Use uncertainty information if available from .cif
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

    symmetries = nothing
    for sym_header in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
        if sym_header in keys(cif)
            sym_table = CIF.get_loop(cif, sym_header)
            symmetries = parse_op.(sym_table[:, sym_header])
        end
    end

    groupnum = nothing
    for group_header in ("_space_group_it_number", "_symmetry_int_tables_number")
        if group_header in keys(cif)
            groupnum = parse(Int, cif[group_header][1])
        end
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

    spacegroup = if !isnothing(hm_symbol)
        if !isnothing(groupnum)
            "'$hm_symbol' ($groupnum)"
        else
            "'$hm_symbol'"
        end
    else
        if !isnothing(groupnum)
            "($groupnum)"
        else
            ""
        end
    end

    if override_symmetry
        # Ignore all symmetry data in the CIF. The use of atom `types` (not site
        # `labels`) is necessary to allow merging of distinct labels into a
        # single symmetry-equivalent site.
        return Crystal(latvecs, positions; types, symprec)
    elseif !isnothing(symmetries)
        # Use explicitly provided symmetries
        return crystal_from_symops(latvecs, positions, labels, symmetries, spacegroup; symprec)
    elseif !isnothing(hall_symbol)
        # Use symmetries for Hall symbol
        return Crystal(latvecs, positions, hall_symbol; types=labels, symprec)
    elseif !isnothing(groupnum)
        # Use symmetries for international group number
        return Crystal(latvecs, positions, groupnum; types=labels, symprec)
    else
        # Infer the symmetries automatically, trusting that distinct CIF labels
        # correspond to symmetry-inequivalent sites.
        return Crystal(latvecs, positions; types=labels, symprec)
    end
end


"""
    set_dipoles_from_mcif!(sys::System, filename::AbstractString)

Load the magnetic supercell according to an MCIF file. System `sys` must already
be resized to the correct supercell dimensions.
"""
function set_dipoles_from_mcif!(sys::System, filename::AbstractString)
    cif = CIF.Cif(Path(filename))
    # For now, assumes there is only one data collection per .cif
    cif = cif[first(keys(cif))]

    a = parse_cif_float(cif["_cell_length_a"][1])
    b = parse_cif_float(cif["_cell_length_b"][1])
    c = parse_cif_float(cif["_cell_length_c"][1])
    α = parse_cif_float(cif["_cell_angle_alpha"][1])
    β = parse_cif_float(cif["_cell_angle_beta"][1])
    γ = parse_cif_float(cif["_cell_angle_gamma"][1])
    supervecs = sys.crystal.latvecs .* sys.latsize
    supervecs2 = lattice_vectors(a, b, c, α, β, γ)
    if !isapprox(supervecs, supervecs2; rtol=sys.crystal.symprec)
        error("""Invalid supercell dimensions,
                   System: $supervecs
                   MCIF:   $supervecs2
                 """)
    end

    sym_header = "_space_group_symop.magn_operation_xyz"
    sym_table = CIF.get_loop(cif, sym_header)
    magn_operations = MSymOp.(sym_table[:, sym_header])

    sym_header = "_space_group_symop.magn_centering_xyz"
    sym_table = CIF.get_loop(cif, sym_header)
    magn_centerings = MSymOp.(sym_table[:, sym_header])

    dip_table = CIF.get_loop(cif, "_atom_site_moment_label")
    labels = String.(dip_table[:, "_atom_site_moment_label"])
    Mxs = parse_cif_float.(dip_table[:, "_atom_site_moment_crystalaxis_x"])
    Mys = parse_cif_float.(dip_table[:, "_atom_site_moment_crystalaxis_y"])
    Mzs = parse_cif_float.(dip_table[:, "_atom_site_moment_crystalaxis_z"])
    moments = Vec3.(zip(Mxs, Mys, Mzs))

    geom_table = CIF.get_loop(cif, "_atom_site_label")
    all_labels = String.(geom_table[:, "_atom_site_label"])
    xs = parse_cif_float.(geom_table[:, "_atom_site_fract_x"])
    ys = parse_cif_float.(geom_table[:, "_atom_site_fract_y"])
    zs = parse_cif_float.(geom_table[:, "_atom_site_fract_z"])
    idxs = indexin(labels, all_labels)
    positions = Vec3.(zip(xs[idxs], ys[idxs], zs[idxs]))

    set_dipoles_from_mcif_aux!(sys; positions, moments, magn_operations, magn_centerings)
end

# The keyword parameters follow the MCIF format. As a consequence, all (x,y,z)
# components are provided in the coordinate system defined by the lattice
# vectors of the supercell. Similarly, the symmetry operations act in this
# coordinate system.
function set_dipoles_from_mcif_aux!(sys; positions, moments, magn_operations, magn_centerings)
    supervecs = sys.crystal.latvecs .* sys.latsize
    
    # Use the zero vector as a marker for unvisited sites
    fill!(sys.dipoles, zero(Vec3))

    for (r, μ) in zip(positions, moments)
        for s1 in magn_operations, s2 in magn_centerings
            # Symmetry operation that composes centering an then a
            # crystallographic operation
            s = s1 * s2

            # Apply symmetry operation to position, and then map to original
            # fractional coordinates
            r_new = s.R * r + s.T
            r_new = orig_crystal(sys).latvecs \ supervecs * r_new

            # Apply symmetry operation to magnetic moment, and then map to
            # Cartesian coordinates. Because the moment is a pseudo-vector, it
            # is invariant to an inversion in space. For this reason, remove the
            # determinant of s.R. If the parity s.p = ±1 is negative, this
            # implies time reversal, which reverses the magnetic dipole.
            μ_new = (s.R * det(s.R) * s.p) * μ
            μ_new = supervecs * μ_new

            site = try 
                position_to_site(sys, r_new)
            catch _
                rethrow(ErrorException("MCIF position $r_new is missing in the chemical cell"))
            end

            # Get spin dipole by inverting the `magnetic_moment` transformation
            dipole = inv(- sys.units.μB * sys.gs[site]) * μ_new
            s_prev = sys.dipoles[site]
            set_dipole!(sys, dipole, site)
            s_prev == zero(Vec3) || s_prev ≈ sys.dipoles[site] || error("Conflicting dipoles at site $site")
        end
    end

    unvisited = [site for site in eachsite(sys) if iszero(sys.dipoles[site])]
    if !isempty(unvisited)
        error("Missing dipoles for sites $(collect(Tuple.(unvisited)))")
    end
end
