"""
    set_dipoles_from_mcif!(sys::System, filename::AbstractString)

Load the magnetic supercell according to an mCIF file. System `sys` must already
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

    # TODO: Tolerance to permutations (with sign flips) of lattice vectors
    if !isapprox(supervecs, supervecs2; rtol=sys.crystal.symprec)
        tol = 0.1 * sys.crystal.symprec # Tolerance might need tuning
        orig_cryst = orig_crystal(sys)

        primvecs = @something orig_cryst.prim_latvecs orig_cryst.latvecs

        suggestion = if all(isinteger.(rationalize.(primvecs \ supervecs2; tol)))
            suggested_shape = rationalize.(orig_cryst.latvecs \ supervecs2; tol)
            suggestion = if isdiag(suggested_shape)
                sz = fractional_vec3_to_string(diag(suggested_shape))
                error("Use `resize_supercell(sys, $sz)` to get compatible system")
            else
                shp = fractional_mat3_to_string(suggested_shape)
                error("Use `reshape_supercell(sys, $shp)` to get compatible system")
            end
        else
            error("""System dimensions are incompatible with mCIF cell,
                         System: $supervecs
                         mCIF:   $supervecs2
                     """)
        end
    end

    oneof(fields...) = findfirstval(in(keys(cif)), fields)
    # The first entry is the IUCR standard field name. If missing, search for
    # alternate field names that appear in legacy files.
    mcif_fields = (;
        magn_operation_xyz=oneof("_space_group_symop_magn_operation.xyz", "_space_group_symop.magn_operation_xyz"),
        magn_centering_xyz=oneof("_space_group_symop_magn_centering.xyz", "_space_group_symop.magn_centering_xyz"),
        moment_label=oneof("_atom_site_moment.label", "_atom_site_moment_label"),
        moment_crystalaxis_x=oneof("_atom_site_moment.crystalaxis_x", "_atom_site_moment_crystalaxis_x"),
        moment_crystalaxis_y=oneof("_atom_site_moment.crystalaxis_y", "_atom_site_moment_crystalaxis_y"),
        moment_crystalaxis_z=oneof("_atom_site_moment.crystalaxis_z", "_atom_site_moment_crystalaxis_z"),
    )

    sym_table = CIF.get_loop(cif, mcif_fields.magn_operation_xyz)
    magn_operations = MSymOp.(sym_table[:, mcif_fields.magn_operation_xyz])
    sym_table = CIF.get_loop(cif, mcif_fields.magn_centering_xyz)
    magn_centerings = MSymOp.(sym_table[:, mcif_fields.magn_centering_xyz])

    dip_table = CIF.get_loop(cif, mcif_fields.moment_label)
    labels = String.(dip_table[:, mcif_fields.moment_label])

    Mxs = parse_cif_float.(dip_table[:, mcif_fields.moment_crystalaxis_x])
    Mys = parse_cif_float.(dip_table[:, mcif_fields.moment_crystalaxis_y])
    Mzs = parse_cif_float.(dip_table[:, mcif_fields.moment_crystalaxis_z])
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

# The keyword parameters follow the mCIF format. As a consequence, all (x,y,z)
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
                rethrow(ErrorException("mCIF position $r_new is missing in the chemical cell"))
            end

            # Get spin dipole by inverting the `magnetic_moment` transformation
            dipole = - sys.gs[site] \ μ_new
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
