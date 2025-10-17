"""
    set_dipoles_from_mcif!(sys::System, filename::AbstractString; symprec=nothing)

Load the magnetic moments from data in the mCIF at `filename`.

The system shape must exactly match the magnetic cell described by the mCIF. If
mismatched, an error message will provide instructions for calling
[`reshape_supercell`](@ref).
"""
function set_dipoles_from_mcif!(sys::System, filename::AbstractString; symprec=nothing)
    cryst = orig_crystal(sys)
    mcif_cryst, symprec = Crystal(filename; symprec, keep_supercell=true)

    if mcif_cryst.sg.number != cryst.sg.number
        error("Inferred mCIF spacegroup $(mcif_cryst.sg.number) differs from $(cryst.sg.number). Try changing `symprec`?")
    end

    # Standard lattice vectors of the two crystals (possibly different Cartesian
    # coordinate systems)
    stdvecs = cryst.latvecs * inv(cryst.sg.setting.R)
    mcif_stdvecs = mcif_cryst.latvecs * inv(mcif_cryst.sg.setting.R)

    # Orthogonal transformation that maps between Cartesian coordinate systems
    R = stdvecs / mcif_stdvecs
    @assert R' * R ≈ I

    # System lattice vectors
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.dims))

    # mCIF lattice vectors rotated into consistent Cartesian frame
    mcif_supervecs_rotated = R * mcif_cryst.latvecs

    # Verify mcif_supervecs are integer multiple of cryst primitive cells.
    primvecs = cryst.latvecs * primitive_cell(cryst)
    @assert all_integer(primvecs \ mcif_supervecs_rotated; tol=1e-12)

    # Check consistency of supercells at appropriate tolerance.
    if !isapprox(supervecs, mcif_supervecs_rotated)
        suggested_shape = rationalize.(cryst.latvecs \ mcif_supervecs_rotated; tol=1e-12)
        if isdiag(suggested_shape)
            diag_strs = number_to_math_string.(diag(suggested_shape))
            sz = "("*join(diag_strs, ", ")*")"
            error("Use `resize_supercell(sys, $sz)` to get compatible system")
        else
            shp = mat3_to_string(suggested_shape)
            error("Use `reshape_supercell(sys, $shp)` to get compatible system")
        end
    end

    (; magn_operations, magn_centerings, moments, positions) = parse_mcif_data(filename)

    # Map from mcif_cryst positions to cryst positions. This works by
    # round-tripping through the "ITA standard" chemical cell.
    map_mcif_coord = MSymOp(inv(cryst.sg.setting) * mcif_cryst.sg.setting)

    # Use the zero vector as a marker for unvisited sites
    fill!(sys.dipoles, zero(Vec3))

    for (r, μ) in zip(positions, moments)
        for s1 in magn_operations, s2 in magn_centerings
            # Symmetry operation that composes centering and then a
            # crystallographic operation
            s = s1 * s2

            # Apply symmetry operation to position, and then map to original
            # fractional coordinates
            r_new = transform(map_mcif_coord * s, r)

            # Apply symmetry operation to magnetic moment, and then map to
            # Cartesian coordinates. Because the moment is a pseudo-vector, it
            # is invariant to an inversion in space. For this reason, remove the
            # determinant of s.R. If the parity s.p = ±1 is negative, this
            # implies time reversal, which reverses the magnetic dipole.
            μ_new = transform_dipole(s, μ)
            μ_new = supervecs * μ_new

            site = try
                # Since we are using mCIF positions and magn symops directly,
                # allow for finite symprec in the position lookup.
                position_to_site(sys, r_new; tol=symprec)
            catch _
                rethrow(ErrorException("Magnetic symops inconsistent with spacegroup symmetry"))
            end

            # Get spin dipole by inverting the `magnetic_moment` transformation
            dipole = - sys.gs[site] \ μ_new
            s_prev = sys.dipoles[site]
            set_dipole!(sys, dipole, site)
            s_prev == zero(Vec3) || s_prev ≈ sys.dipoles[site] || error("Magnetic symops internally inconsistent")
        end
    end

    any(iszero, sys.dipoles) && error("Magnetic orbits do not cover all sites")
    return
end


function parse_mcif_data(filename)
    cif = CIF.Cif(filename)
    # For now, assumes there is only one data collection per .cif
    cif = cif[first(keys(cif))]

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

    # All positions in given in fractional coordinates of the mCIF cell, i.e.,
    # supervecs
    idxs = indexin(labels, all_labels)
    positions = Vec3.(zip(xs[idxs], ys[idxs], zs[idxs]))

    return (; magn_operations, magn_centerings, moments, positions)
end
