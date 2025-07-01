function characteristic_length_between_atoms(cryst::Crystal)
    # Detect if atom displacements are on a submanifold (aligned line or plane)
    ps = cryst.positions[1:end-1] .- Ref(cryst.positions[end])
    any_nonzero = map(1:3) do i
        any(p -> !iszero(p[i]), ps)
    end
    vecs = eachcol(cryst.latvecs)[findall(any_nonzero)]

    # Take nth root of appropriate hypervolume per atom
    if length(vecs) == 0
        ℓ = Inf                            # For a single atom, use ℓ0 below
    elseif length(vecs) == 1
        ℓ = norm(vecs[1]) / natoms(cryst)  # Atoms aligned with single lattice vector
    elseif length(vecs) == 2
        ℓ = sqrt(norm(vecs[1] × vecs[2]) / natoms(cryst))
    elseif length(vecs) == 3
        ℓ = cbrt(abs(det(cryst.latvecs)) / natoms(cryst))
    else
        error("Internal error")
    end

    # An upper bound is the norm of the smallest lattice vector.
    ℓ0 = minimum(norm.(eachcol(cryst.latvecs)))

    return min(ℓ0, ℓ)
end

# Like `reference_bonds` but supply a number of bonds
function reference_bonds_upto(cryst, nbonds, ndims)
    # Calculate heuristic maximum distance
    min_a = minimum(norm.(eachcol(cryst.latvecs)))
    nclasses = length(unique(cryst.classes))
    max_dist = 2 * min_a * (nbonds / (nclasses*natoms(cryst)))^(1/ndims)

    # Find bonds up to distance, without self-bonds
    refbonds = filter(reference_bonds(cryst, max_dist)) do b
        return !(b.i == b.j && iszero(b.n))
    end

    # Verify max_dist heuristic
    if length(refbonds) > 10nbonds
        @warn "Found $(length(refbonds)) bonds using max_dist of $max_dist"
    end

    return first(refbonds, nbonds)
end

function propagate_reference_bond_for_cell(cryst, b_ref)
    symops = Sunny.canonical_group_order(cryst.sg.symops)

    found = map(_ -> Bond[], cryst.positions)
    for s in symops
        b = Sunny.transform(cryst, s, b_ref)
        # If this bond hasn't been found, add it to the list
        if !(b in found[b.i]) && !(reverse(b) in found[b.j])
            push!(found[b.i], b)
        end
    end

    return reduce(vcat, found)
end

# Get largest exchange interaction scale. For symmetric part, this is the
# largest eigenvalue. For antisymmetric part, this is an empirical rescaling of
# the norm of the DM vector. (Note that a typical system has small DM vector
# magnitude relative to the symmetric exchange, and the heuristics for visual
# size are taking this into account.)
function exchange_magnitude(interactions)
    ret = -Inf
    for int in interactions, pc in int.pair
        J = pc.bilin * Mat3(I)
        sym = maximum(abs.(eigvals(Hermitian(J+J')/2)))
        dm = norm(Sunny.extract_dmvec(J))
        ret = max(ret, sym + 2dm)
    end
    return ret
end

# Return an axis scaling and quaternion rotation corresponding to J
function exchange_decomposition(J)
    # Absolute value of eigenvalues control scaling of ellipsoidal axis, with
    # ellipsoid volume depicting interaction strength.
    vals, vecs = eigen(Hermitian(J+J')/2)

    # If vecs includes a reflection, then permute columns
    if det(vecs) < 0
        vals = [vals[2], vals[1], vals[3]]
        vecs = hcat(vecs[:,2], vecs[:,1], vecs[:,3])
    end

    # Now vecs is a pure rotation
    @assert vecs'*vecs ≈ I && det(vecs) ≈ 1

    # Quaternion that rotates Cartesian coordinates into principle axes of J.
    axis, angle = Sunny.matrix_to_axis_angle(Mat3(vecs))
    q = iszero(axis) ? Makie.Quaternionf(0,0,0,1) : Makie.qrotation(axis, angle)

    return (vals, q)
end

function draw_exchange_geometries(; ax, visible, ionradius, pts, scaled_exchanges)

    ### Ellipsoids for symmetric exchanges

    # Dimensionless scalings and rotations associated with principle axes
    decomps = exchange_decomposition.(scaled_exchanges)
    scalings = map(x -> x[1], decomps)
    rotation = map(x -> x[2], decomps)

    # Enlarge scalings so that the maximum scaling _cubed_ denotes magnitude
    scalings = map(scalings) do scal
        szmax = maximum(abs.(scal))
        cbrt(szmax) * (scal/szmax)
    end

    markersize = map(scalings) do scal
        # Make sure ellipsoids don't get flattened to zero
        szmax = maximum(abs.(scal))
        ionradius * Makie.Vec3f([max(abs(x), szmax/4) for x in scal])
    end

    # Draw ellipsoidal bounding box
    color = map(scalings) do x
        y = sum(x) / sum(abs.(x)) # -1 ≤ y ≤ 1
        c = 0.8
        d = c+(1-c)*abs(y) # c ≤ d ≤ 1
        y > 0 ? Makie.RGBf(c, c, d) : Makie.RGBf(d, c, c)
    end
    Makie.meshscatter!(ax, pts; color, markersize, rotation, specular=0, diffuse=1.5, visible, inspectable=false)

    # Draw dots using cylinders
    cylinders = map(eachcol(Sunny.Mat3(I))) do x
        p = Makie.GeometryBasics.Point(x...)
        Makie.GeometryBasics.Cylinder(-p, p, 0.3)
    end
    for dim in 1:3
        color = map(scalings) do x
            x[dim] < 0 ? :red : :blue
        end

        # Apply some additional scaling so that all the dots on a given
        # ellipsoid have a roughly constant linear size
        rescalings = map(scalings) do x
            c = sqrt(abs(x[dim]) / maximum(abs.(x)))
            [dim == 1 ? 1 : c,
             dim == 2 ? 1 : c,
             dim == 3 ? 1 : c]
        end
        markersize2 = [ms .* rs for (ms, rs) in zip(markersize, rescalings)]

        Makie.meshscatter!(ax, pts; color, markersize=markersize2, rotation, marker=cylinders[dim], visible, inspectable=false)
    end

    ### Cones for DM vectors

    # DM vectors normalized by the overall exchange scale
    dmvecs = Sunny.extract_dmvec.(scaled_exchanges)

    # Filter nonzero DM vectors and retain associated exchange scalings
    nonzeros_indices = findall(norm(dmvec) > 0 for dmvec in dmvecs)
    dm_pts = pts[nonzeros_indices]
    dmvecs = dmvecs[nonzeros_indices]
    scalings_nonzero = scalings[nonzeros_indices]

    # The largest possible ellipsoid occurs in the case of `scalings ==
    # [1,1,1]`, yielding a sphere with size `ionradius`.
    ellipsoid_radii = @. ionradius * norm(scalings_nonzero) / √3

    # Size of cone scales like cube root of DM vector magnitude
    markersize = @. 2ionradius * cbrt(norm(dmvecs))

    # Offset cone according to size of exchange ellipsoid
    dm_pts = @. Makie.Point3f(dm_pts + 1.1ellipsoid_radii * normalize(dmvecs))

    # Draw cones
    marker = Makie.GeometryBasics.Cone(Makie.Point3f(0, 0, 0), Makie.Point3f(0, 0, 1), 0.5)
    Makie.meshscatter!(ax, dm_pts; marker, markersize, rotation=dmvecs, color=:magenta,
                       specular=0.0, diffuse=1.15, visible, inspectable=false)
end

function draw_bonds(; ax, visible, ionradius, exchange_mag, cryst, interactions, bonds, refbond, color)
    # The bond is directed if it allows for a DM exchange. This will be
    # visualized with an arrow.
    isdirected = begin
        basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, refbond)
        any(J -> J ≈ -J', basis)
    end

    if isempty(bonds)
        bonds = [refbond]
        show_ghosts = false
    else
        show_ghosts = true
    end

    # Map each bond to line segments in global coordinates
    segments = map(bonds) do b
        (; ri, rj) = Sunny.BondPos(cryst, b)
        Makie.Point3f.(Ref(cryst.latvecs) .* (ri, rj))
    end

    # Add periodic "ghost" bonds if we're not showing only reference bonds
    if show_ghosts
        # Indices for the bonds which most be repeated
        ghosts = findall(b -> !iszero(b.n), bonds)

        # Concatenate ghosts to the end of arrays
        append!(bonds, bonds[ghosts])
        append!(color, color[ghosts])

        # Ghost bonds are offset by -n multiples of lattice vectors
        ghost_segments = map(ghosts) do i
            offset = - cryst.latvecs * bonds[i].n
            Makie.Point3f.(segments[i] .+ Ref(offset))
        end
        append!(segments, ghost_segments)
    end

    # String for each bond b′. Like print_bond(b′), but shorter.
    bond_labels = map(bonds) do b
        dist = Sunny.global_distance(cryst, b)
        dist_str = Sunny.number_to_simple_string(dist; digits=4, atol=1e-12)

        if isnothing(interactions)
            basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, b; b_ref=refbond)
            basis_strs = Sunny.coupling_basis_strings(zip('A':'Z', basis); digits=4, atol=1e-12)
            J_matrix_str = Sunny.formatted_matrix(basis_strs; prefix="J:  ")
            antisym_basis_idxs = findall(J -> J ≈ -J', basis)
            if !isempty(antisym_basis_idxs)
                antisym_basis_strs = Sunny.coupling_basis_strings(collect(zip('A':'Z', basis))[antisym_basis_idxs]; digits=4, atol=1e-12)
                dmvecstr = join([antisym_basis_strs[2,3], antisym_basis_strs[3,1], antisym_basis_strs[1,2]], ", ")
                J_matrix_str *= "\nDM: [$dmvecstr]"
            end
        else
            c = Sunny.search_pair_couplings_for_bond(interactions[b.i].pair, b)
            J = isnothing(c) ? zero(Mat3) : c.bilin * Mat3(I)
            basis_strs = Sunny.number_to_simple_string.(J; digits=3)
            J_matrix_str = Sunny.formatted_matrix(basis_strs; prefix="J:  ")
            if J ≉ J'
                dmvec = Sunny.extract_dmvec(J)
                dmvecstr = join(Sunny.number_to_simple_string.(dmvec; digits=3), ", ")
                J_matrix_str *= "\nDM: [$dmvecstr]"
            end
            if !isnothing(c) && (!iszero(c.biquad) || !isempty(c.general.data))
                J_matrix_str *= "\n  + higher order terms"
            end
        end

        return """
            $b
            Distance $dist_str
            $J_matrix_str
            """
    end
    inspector_label(_plot, index, _position) = bond_labels[index]

    # Draw cylinders or arrows for each bond
    shaftradius = 0.125ionradius # * hasarrowhead
    tipradius = 1.8shaftradius
    tiplength = isdirected ? 4.4tipradius : 0.0
    dirs = [normalize(rj-ri) for (ri, rj) in segments]
    pts0 = @. getindex.(segments, 1) + ionradius*dirs
    pts1 = @. getindex.(segments, 2) - ionradius*dirs
    Makie.arrows3d!(ax, pts0, pts1-pts0; markerscale=1, minshaftlength=0, tiplength, tipradius, shaftradius,
                    color, diffuse=3, transparency=true, visible, inspectable=true, inspector_label)

    # Draw exchange interactions if data is available
    if exchange_mag > 0
        pts = [(ri+rj)/2 for (ri, rj) in segments]
        exchanges = map(bonds) do b
            Sunny.get_exchange_from_interactions(interactions[b.i], b)
        end
        draw_exchange_geometries(; ax, visible, ionradius, pts, scaled_exchanges=exchanges/exchange_mag)
    end

    return
end


# Return true if `type` doesn't uniquely identify the site equivalence class
function is_type_degenerate(cryst, i)
    typ = cryst.types[i]
    same_typ = findall(==(typ), cryst.types)
    return !allequal(cryst.classes[same_typ])
end

# Construct atom labels for use in DataInspector
function label_atoms(cryst; ismagnetic, sys)
    return map(1:natoms(cryst)) do i
        typ = cryst.types[i]
        rstr = Sunny.fractional_vec3_to_string(cryst.positions[i])
        ret = []

        (; multiplicity, letter) = Sunny.get_wyckoff(cryst, i)
        wyckstr = "$multiplicity$letter"
        typstr = isempty(typ) ? "" : "'$typ', "
        push!(ret, typstr * "Wyckoff $wyckstr, $rstr")

        if ismagnetic
            if isnothing(sys)
                # See similar logic in print_site()
                refatoms = [b.i for b in Sunny.reference_bonds(cryst, 0.0)]
                i_ref = Sunny.findfirstval(i_ref -> Sunny.is_related_by_symmetry(cryst, i, i_ref), refatoms)
                R_site = Sunny.rotation_between_sites(cryst, i, i_ref)
                push!(ret, Sunny.allowed_g_tensor_string(cryst, i_ref; R_site, prefix="Aniso: ", digits=8, atol=1e-12))
            else
                site = Sunny.map_atom_to_other_system(cryst, i, sys)
                stvexp = Sunny.get_stevens_expansion_at(sys, site)
                aniso = Sunny.unrenormalize_quadratic_anisotropy(stvexp, sys, site)
                basis_strs = Sunny.number_to_simple_string.(aniso; digits=3)
                push!(ret, Sunny.formatted_matrix(basis_strs; prefix="Aniso: "))
                if !iszero(stvexp.c4) || !iszero(stvexp.c6)
                    push!(ret, "  + higher order terms")
                end
            end
        end
        join(ret, "\n")
    end
end

function scaled_dipole_to_arrow_length(dipole, lengthscale, tiplength)
    s = norm(dipole)
    dir = dipole / s
    baselen = s * lengthscale
    # If dipole becomes very short, its magnitude will be depicted as the
    # _volume_ of the arrow tip. This is achieved by scaling tiplength (space
    # reserved for arrow tip) by a reduction factor c ~ cbrt(s) ≤ 1
    r = 4baselen/tiplength
    c = cbrt(min(r, 1))
    return (baselen + c * tiplength) * dir
end

function draw_atoms_or_dipoles(; ax, full_crystal_toggle, dipole_menu, cryst, sys, class_colors, ionradius, ndims, ghost_radius)
    selection = isnothing(dipole_menu) ? Makie.Observable("No dipoles") : dipole_menu.selection
    show_spin_dipoles = Makie.lift(==("Spin dipoles"), selection)
    show_magn_dipoles = Makie.lift(==("Magnetic dipoles"), selection)

    # Draw magnetic and non-magnetic ions
    for ismagnetic in (false, true)
        if ismagnetic
            xtal = cryst
            relevant_classes = cryst.classes
        else
            isnothing(cryst.root) && continue
            # Relevant classes are those present in cryst.root but not cryst
            xtal = cryst.root
            relevant_classes = setdiff(xtal.classes, cryst.classes)
        end

        for isghost in (false, true)
            if isghost
                (idxs, offsets) = Sunny.all_offsets_within_distance(xtal.latvecs, xtal.positions, cell_center(ndims); max_dist=ghost_radius, nonzeropart=true)
                alpha = 0.08
            else
                idxs = eachindex(xtal.positions)
                offsets = [zero(Vec3) for _ in idxs]
                alpha = 1.0
            end

            # Reduce idxs and offsets to include only atom indices according to
            # `relevant_classes`, as set by `ismagnetic`
            downselect = findall(in(relevant_classes), xtal.classes[idxs])
            isempty(downselect) && continue
            idxs = idxs[downselect]
            offsets = offsets[downselect]

            # Information for drawing atoms in xtal labeled by idxs
            color = [(class_colors[c], alpha) for c in xtal.classes[idxs]]
            rs = xtal.positions[idxs] .+ offsets
            pts = [xtal.latvecs * r for r in rs]

            # Labels for non-ghost atoms
            inspector_label = nothing
            if !isghost
                labels = label_atoms(xtal; ismagnetic, sys)[idxs]
                inspector_label = (_plot, index, _position) -> labels[index]
            end

            # Show dipoles. Mostly consistent with code in plot_spins.
            if !isnothing(sys) && ismagnetic
                sites = Sunny.position_to_site.(Ref(sys), rs)
                g0 = norm(sys.gs) / sqrt(length(sys.gs) * 3)
                N0 = norm(sys.Ns) / sqrt(length(sys.Ns))
                s0 = (N0 - 1) / 2
                spin_dipoles = sys.dipoles[sites] / s0
                magn_dipoles = magnetic_moment.(Ref(sys), sites) / (s0*g0)
                for (dipoles, visible) in [(spin_dipoles, show_spin_dipoles), (magn_dipoles, show_magn_dipoles)]
                    a0 = 5ionradius
                    shaftradius = 0.06a0
                    tipradius = 0.2a0
                    tiplength = 0.4a0
                    lengthscale = 0.6a0
                    vecs = scaled_dipole_to_arrow_length.(dipoles, lengthscale, tiplength)

                    # Draw arrows
                    shaftcolor = (:white, alpha)
                    tipcolor = (:gray, alpha)
                    Makie.arrows3d!(ax, pts, vecs; align=0.37, markerscale=1, tipradius, shaftradius,
                                    tiplength, minshaftlength=0, tipcolor, shaftcolor, diffuse=1.15,
                                    transparency=isghost, visible, inspectable=false)
                end
            end

            # Show atoms as spheres
            markersize = ionradius * (ismagnetic ? 1 : 1/2)

            visible = ismagnetic ? true : full_crystal_toggle.active
            Makie.meshscatter!(ax, pts; markersize, color, diffuse=1.15, transparency=isghost, visible, inspectable=!isghost, inspector_label)

            # White numbers for real, magnetic ions
            if !isghost && ismagnetic
                text = repr.(eachindex(pts))
                Makie.text!(ax, pts; text, color=:white, fontsize=16, align=(:center, :center), depth_shift=-1f0)
            end
        end
    end
end


function Sunny.view_crystal(cryst::Crystal, max_dist::Number)
    @warn "view_crystal(cryst, max_dist) is deprecated! Use `view_crystal(cryst)` instead. See also optional `ghost_radius` argument."
    Sunny.view_crystal(cryst; ghost_radius=max_dist)
end

"""
    view_crystal(crystal::Crystal; refbonds=10, orthographic=false, ghost_radius=nothing, ndims=3, compass=true)
    view_crystal(sys::System; ...)

Launches a graphical user interface to visualize the [`Crystal`](@ref) unit
cell. If a [`System`](@ref) is provided, then the 3×3 exchange matrices for each
bond will be depicted graphically.

 - `refbonds`: By default, calculate up to 10 reference bonds using the
   `reference_bonds` function. An explicit list of reference bonds may also be
   provided.
 - `orthographic`: Use orthographic camera perspective.
 - `ghost_radius`: Show periodic images up to a given distance. Defaults to the
   cell size.
 - `ndims`: Spatial dimensions of system (1, 2, or 3).
 - `compass`: If true, draw Cartesian axes in bottom left.
"""
function Sunny.view_crystal(cryst::Crystal; refbonds=10, orthographic=false, ghost_radius=nothing, ndims=3, compass=true, size=(768, 512), dims=nothing)
    isnothing(dims) || error("Use notation `ndims=$dims` instead of `dims=$dims`")
    view_crystal_aux(cryst, nothing; refbonds, orthographic, ghost_radius, ndims, compass, size)
end

function Sunny.view_crystal(sys::System; refbonds=10, orthographic=false, ghost_radius=nothing, ndims=3, compass=true, size=(768, 512), dims=nothing)
    isnothing(dims) || error("Use notation `ndims=$dims` instead of `dims=$dims`")
    Sunny.is_homogeneous(sys) || error("Cannot plot interactions for inhomogeneous system.")
    view_crystal_aux(orig_crystal(sys), sys;
                     refbonds, orthographic, ghost_radius, ndims, compass, size)
end

function view_crystal_aux(cryst, sys; refbonds, orthographic, ghost_radius, ndims, compass, size)
    interactions = isnothing(sys) ? nothing : Sunny.interactions_homog(something(sys.origin, sys))

    # Dict that maps atom class to color
    class_colors = build_class_colors(cryst)

    # Distance to show periodic images
    if isnothing(ghost_radius)
        ghost_radius = cell_diameter(cryst.latvecs, ndims)/2
    end

    # Use provided reference bonds or find from symmetry analysis
    if refbonds isa Number
        @assert isinteger(refbonds)
        activate_refbonds = false
        refbonds = reference_bonds_upto(cryst, Int(refbonds), ndims)
    elseif refbonds isa AbstractArray{Bond}
        activate_refbonds = true
    else
        error("Parameter `refbonds` must be an integer or a `Bond` list.")
    end

    refbonds_dists = [Sunny.global_distance(cryst, b) for b in refbonds]

    # Radius of the magnetic ions. Sets a length scale for other objects too.
    ionradius = let
        # The root crystal may contain non-magnetic ions. If present, these
        # should reduce the characteristic length scale.
        ℓ0 = characteristic_length_between_atoms(something(cryst.root, cryst))
        # If there exists a very short bond distance, then appropriately reduce the
        # length scale
        ℓ0 = min(ℓ0, 0.8minimum(refbonds_dists))
        # Small enough to fit everything
        0.2ℓ0
    end

    fig = Makie.Figure(; size)

    # Main scene
    ax = Makie.LScene(fig[1, 1], show_axis=false)

    # Show Cartesian axes, with link to main camera
    if compass
        axcompass = add_cartesian_compass(fig, ax)
    end

    # Set of widgets
    widget_list = Makie.GridLayout(; tellheight=false, valign=:top)
    fig[1, 2] = widget_list
    widget_cnt = 0
    fontsize = 16

    # Controls for camera perspective
    menu = Makie.Menu(fig; options=["Perspective", "Orthographic"], default=(orthographic ? "Orthographic" : "Perspective"), fontsize)
    button = Makie.Button(fig; label="Reset", fontsize)
    Makie.onany(button.clicks, menu.selection; update=true) do _, mselect
        orthographic = mselect == "Orthographic"
        # Zoom out a little bit extra according to nn bond distance
        ℓ0=minimum(refbonds_dists)
        orient_camera!(ax, cryst.latvecs; ghost_radius, orthographic, ndims, ℓ0)
        compass && register_compass_callbacks(axcompass, ax)
    end
    widget_list[widget_cnt+=1, 1] = Makie.hgrid!(menu, button)

    # Controls for dipoles
    if !isnothing(sys)
        dipole_menu = Makie.Menu(fig; options=["No dipoles", "Spin dipoles", "Magnetic dipoles"], fontsize)
        widget_list[widget_cnt+=1, 1] = dipole_menu
    else
        dipole_menu = nothing
    end

    # Set up grid of toggles
    toggle_grid = widget_list[widget_cnt+=1,1] = Makie.GridLayout()
    toggle_cnt = 0
    buttoncolor = Makie.RGB(0.2, 0.2, 0.2)
    framecolor_active = Makie.RGB(0.7, 0.7, 0.7)
    framecolor_inactive = Makie.RGB(0.9, 0.9, 0.9)

    # Toggle for non-magnetic ions
    if !isnothing(cryst.root)
        full_crystal_toggle = Makie.Toggle(fig; active=true, buttoncolor, framecolor_inactive, framecolor_active)
        toggle_grid[toggle_cnt+=1, 1:2] = [full_crystal_toggle, Makie.Label(fig, "Full crystal"; fontsize, halign=:left)]
    else
        full_crystal_toggle = nothing
    end
    draw_atoms_or_dipoles(; ax, full_crystal_toggle, dipole_menu, cryst, sys, class_colors, ionradius, ndims, ghost_radius)

    exchange_mag = isnothing(interactions) ? 0.0 : exchange_magnitude(interactions)

    # Toggle on/off atom reference bonds
    bond_colors = [getindex_cyclic(seaborn_bright, i) for i in eachindex(refbonds)]
    toggle = Makie.Toggle(fig; active=activate_refbonds, buttoncolor, framecolor_inactive, framecolor_active)
    for (refbond, bond_color) in zip(refbonds, bond_colors)
        color = set_alpha(bond_color, 0.25)
        draw_bonds(; ax, visible=toggle.active, ionradius, exchange_mag, cryst, interactions, bonds=[], refbond, color)
    end
    toggle_grid[toggle_cnt+=1, 1:2] = [toggle, Makie.Label(fig, "Reference bonds"; fontsize, halign=:left)]

    # Toggle on/off bonds within each class
    for (i, (refbond, bond_color)) in enumerate(zip(refbonds, bond_colors))
        active = (i == 1)
        framecolor_active = set_alpha(bond_color, 0.7)
        framecolor_inactive = set_alpha(bond_color, 0.15)
        toggle = Makie.Toggle(fig; active, buttoncolor, framecolor_inactive, framecolor_active)
        bonds = propagate_reference_bond_for_cell(cryst, refbond)
        color = fill(set_alpha(bond_color, 0.25), length(bonds))
        draw_bonds(; ax, visible=toggle.active, ionradius, exchange_mag, cryst, interactions, bonds, refbond, color)
        bondstr = "Bond($(refbond.i), $(refbond.j), $(refbond.n))"
        toggle_grid[toggle_cnt+=1, 1:2] = [toggle, Makie.Label(fig, bondstr; fontsize, halign=:left)]
    end

    # Show cell volume
    Makie.linesegments!(ax, cell_wireframe(cryst.latvecs, ndims); color=:teal, linewidth=1.5, inspectable=false)
    pos = [(3/4)*Makie.Point3f(p) for p in eachcol(cryst.latvecs)[1:ndims]]
    text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:ndims]
    Makie.text!(ax, pos; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
                glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-1f0)

    # Add inspector for pop-up information. Use a monospaced font provided
    # available in Makie.jl/assets/fonts/. The depth needs to be almost exactly
    # 1e4, but not quite, otherwise only a white background will be shown.
    Makie.DataInspector(ax; indicator_color=:gray, fontsize, font="Deja Vu Sans Mono", depth=(1e4 - 1))

    return fig
end
