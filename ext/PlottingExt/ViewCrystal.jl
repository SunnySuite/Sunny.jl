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

function anisotropy_on_site(sys, i)
    interactions = isnothing(sys) ? nothing : Sunny.interactions_homog(something(sys.origin, sys))
    onsite = interactions[i].onsite
    if onsite isa Sunny.HermitianC64
        onsite = Sunny.StevensExpansion(Sunny.matrix_to_stevens_coefficients(onsite))
    end
    return onsite :: Sunny.StevensExpansion
end

# Get the quadratic anisotropy as a 3×3 exchange matrix for atom `i` in the
# chemical cell.
function quadratic_anisotropy(sys, i)
    # Get certain Stevens expansion coefficients
    (; c0, c2) = anisotropy_on_site(sys, i)

    # Undo RCS renormalization for quadrupolar anisotropy for spin-s
    if sys.mode == :dipole
        s = (sys.Ns[i] - 1) / 2
        c2 = c2 / Sunny.rcs_factors(s)[2] # Don't mutate c2 in-place!
    end

    # Stevens quadrupole operators expressed as 3×3 spin bilinears
    quadrupole_basis = [
        [1 0 0; 0 -1 0; 0 0 0],    # 𝒪₂₂  = SˣSˣ - SʸSʸ
        [0 0 1; 0 0 0; 1 0 0] / 2, # 𝒪₂₁  = (SˣSᶻ + SᶻSˣ)/2
        [-1 0 0; 0 -1 0; 0 0 2],   # 𝒪₂₀  = 2SᶻSᶻ - SˣSˣ - SʸSʸ
        [0 0 0; 0 0 1; 0 1 0] / 2, # 𝒪₂₋₁ = (SʸSᶻ + SᶻSʸ)/2
        [0 1 0; 1 0 0; 0 0 0],     # 𝒪₂₋₂ = SˣSʸ + SʸSˣ
    ]

    # The c0 coefficient incorporates a factor of S². For quantum spin
    # operators, S² = s(s+1) I. For the large-s classical limit, S² = s² is a
    # scalar.
    S² = if sys.mode == :dipole_uncorrected
        # Undoes extraction in `operator_to_stevens_coefficients`. Note that
        # spin magnitude s² is set to κ², as originates from `onsite_coupling`
        # for p::AbstractPolynomialLike.
        sys.κs[i]^2
    else
        # Undoes extraction in `matrix_to_stevens_coefficients` where 𝒪₀₀ = I.
        s = (sys.Ns[i]-1) / 2
        s * (s+1)
    end

    return c2' * quadrupole_basis + only(c0) * I / S²
end

function coupling_on_bond(interactions, b)
    isnothing(interactions) && return zero(Mat3)
    pairs = interactions[b.i].pair
    indices = findall(pc -> pc.bond == b, pairs)
    return isempty(indices) ? nothing : pairs[only(indices)]
end

# Get the 3×3 exchange matrix for bond `b`
function exchange_on_bond(interactions, b)
    coupling = coupling_on_bond(interactions, b)
    return isnothing(coupling) ? zero(Mat3) : coupling.bilin * Mat3(I)
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

function draw_exchange_geometries(; ax, obs, ionradius, pts, scaled_exchanges)

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
    o = Makie.meshscatter!(ax, pts; color, markersize, rotation, specular=0, diffuse=1.5, inspectable=false)
    Makie.connect!(o.visible, obs)

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
    
        o = Makie.meshscatter!(ax, pts; color, markersize=markersize2, rotation, marker=cylinders[dim], inspectable=false)
        Makie.connect!(o.visible, obs)            
    end

    ### Cones for DM vectors. Because they tend to be weaker in magnitude,
    ### we apply some heuristic amplification to the arrow size.

    dmvecs = Sunny.extract_dmvec.(scaled_exchanges)
    dirs = @. Makie.Vec3f(normalize(dmvecs))
    # The largest possible ellipsoid occurs in the case of `scalings ==
    # [1,1,1]`, yielding a sphere with size `ionradius`.
    ellipsoid_radii = @. ionradius * norm(scalings) / √3
    arrowsize = @. 2ionradius * cbrt(norm(dmvecs)) # size of arrow head
    dm_pts = @. pts + 1.1ellipsoid_radii * dirs
    o = Makie.arrows!(ax, dm_pts, dirs; lengthscale=0, arrowsize, diffuse=1.15, color=:magenta, specular=0.0, inspectable=false) 
    Makie.connect!(o.visible, obs)
end

function draw_bonds(; ax, obs, ionradius, exchange_mag, cryst, interactions, bonds, refbonds, color)

    # Map each bond to line segments in global coordinates
    segments = map(bonds) do b
        (; ri, rj) = Sunny.BondPos(cryst, b)
        Makie.Point3f.(Ref(cryst.latvecs) .* (ri, rj))
    end

    # If the bonds are distinct from the refbonds, then add periodic "ghost" images
    if bonds !== refbonds
        # Indices for the bonds which most be repeated
        ghosts = findall(b -> !iszero(b.n), bonds)

        # Concatenate ghosts to the end of arrays
        bonds = vcat(bonds, bonds[ghosts])
        refbonds = vcat(refbonds, refbonds[ghosts])
        color = vcat(color, color[ghosts])

        # Ghost bonds are offset by -n multiples of lattice vectors
        segments = vcat(segments, map(ghosts) do i
            offset = - cryst.latvecs * bonds[i].n
            segments[i] .+ Ref(offset)
        end)
    end

    # String for each bond b′. Like print_bond(b′), but shorter.
    bond_labels = map(zip(bonds, refbonds)) do (b, b_ref)
        dist = Sunny.global_distance(cryst, b)
        dist_str = Sunny.number_to_simple_string(dist; digits=4, atol=1e-12)

        if isnothing(interactions)
            basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, b; b_ref)
            basis_strs = Sunny.coupling_basis_strings(zip('A':'Z', basis); digits=4, atol=1e-12)
            J_matrix_str = Sunny.formatted_matrix(basis_strs; prefix="J:  ")
            antisym_basis_idxs = findall(J -> J ≈ -J', basis)
            if !isempty(antisym_basis_idxs)
                antisym_basis_strs = Sunny.coupling_basis_strings(collect(zip('A':'Z', basis))[antisym_basis_idxs]; digits=4, atol=1e-12)
                dmvecstr = join([antisym_basis_strs[2,3], antisym_basis_strs[3,1], antisym_basis_strs[1,2]], ", ")
                J_matrix_str *= "\nDM: [$dmvecstr]"
            end
        else
            J = exchange_on_bond(interactions, b)
            basis_strs = Sunny.number_to_simple_string.(J; digits=3)
            J_matrix_str = Sunny.formatted_matrix(basis_strs; prefix="J:  ")
            if J ≉ J'
                dmvec = Sunny.extract_dmvec(J)
                dmvecstr = join(Sunny.number_to_simple_string.(dmvec; digits=3), ", ")
                J_matrix_str *= "\nDM: [$dmvecstr]"
            end
            c = coupling_on_bond(interactions, b)
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

    # A bond has an arrowhead if it allows DM interactions
    hasarrowhead = map(bonds) do b
        basis = Sunny.basis_for_symmetry_allowed_couplings(cryst, b)
        any(J -> J ≈ -J', basis)
    end

    # Draw cylinders or arrows for each bond
    linewidth = 0.25ionradius
    arrowwidth = 1.8linewidth
    arrowlength = 2.2arrowwidth
    disps = [rj-ri for (ri, rj) in segments]
    dirs = normalize.(disps)
    pts = @. getindex.(segments, 1) + ionradius*dirs
    arrowsize = hasarrowhead .* Ref(Makie.Vec3f(arrowwidth, arrowwidth, arrowlength))
    lengthscale = @. norm(disps) - 2ionradius - hasarrowhead*arrowlength
    o = Makie.arrows!(ax, pts, dirs; arrowsize, lengthscale, linewidth, color, diffuse=3,
                      transparency=true, inspectable=true, inspector_label)
    Makie.connect!(o.visible, obs)

    # Draw exchange interactions if data is available
    if exchange_mag > 0
        pts = [(ri+rj)/2 for (ri, rj) in segments]
        exchanges = exchange_on_bond.(Ref(interactions), bonds)
        draw_exchange_geometries(; ax, obs, ionradius, pts, scaled_exchanges=exchanges/exchange_mag)
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
                aniso = quadratic_anisotropy(sys, i)
                basis_strs = Sunny.number_to_simple_string.(aniso; digits=3)
                push!(ret, Sunny.formatted_matrix(basis_strs; prefix="Aniso: "))
                (; c4, c6) = anisotropy_on_site(sys, i)
                if !iszero(c4) || !iszero(c6)
                    push!(ret, "  + higher order terms")
                end
            end
        end
        join(ret, "\n")
    end
end

function draw_atoms_or_dipoles(; ax, full_crystal_toggle, dipole_menu, cryst, sys, class_colors, ionradius, ndims, ghost_radius)
    selection = isnothing(dipole_menu) ? Makie.Observable("No dipoles") : dipole_menu.selection
    show_spin_dipoles = Makie.lift(==("Spin dipoles"), selection)
    show_magn_dipoles = Makie.lift(==("Magnetic dipoles"), selection)
    show_atom_spheres = Makie.lift(==("No dipoles"), selection)

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
                for (dipoles, obs) in [(spin_dipoles, show_spin_dipoles), (magn_dipoles, show_magn_dipoles)]
                    a0 = 5ionradius
                    arrowsize = 0.4a0
                    linewidth = 0.12a0
                    lengthscale = 0.6a0
                    markersize = 0.9ionradius
                    arrow_fractional_shift = 0.6

                    vecs = lengthscale * dipoles
                    pts_shifted = pts - arrow_fractional_shift * vecs

                    # Draw arrows
                    linecolor = (:white, alpha)
                    arrowcolor = (:gray, alpha)
                    o = Makie.arrows!(ax, Makie.Point3f.(pts_shifted), Makie.Vec3f.(vecs); arrowsize, linewidth, linecolor, arrowcolor, diffuse=1.15, transparency=isghost, inspectable=false)
                    Makie.connect!(o.visible, obs)

                    # Small sphere inside arrow to mark atom position
                    o = Makie.meshscatter!(ax, pts; markersize, color, diffuse=1.15, transparency=isghost, inspectable=!isghost, inspector_label)
                    Makie.connect!(o.visible, obs)
                end
            end

            # Show atoms as spheres
            markersize = ionradius * (ismagnetic ? 1 : 1/2)
            o = Makie.meshscatter!(ax, pts; markersize, color, diffuse=1.15, transparency=isghost, inspectable=!isghost, inspector_label)
            Makie.connect!(o.visible, ismagnetic ? show_atom_spheres : full_crystal_toggle.active)

            # White numbers for real, magnetic ions
            if !isghost && ismagnetic
                text = repr.(eachindex(pts))
                o = Makie.text!(ax, pts; text, color=:white, fontsize=16, align=(:center, :center), depth_shift=-1f0)
                !ismagnetic && Makie.connect!(o.visible, full_crystal_toggle.active)
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
        custombonds = false
        refbonds = reference_bonds_upto(cryst, Int(refbonds), ndims)
    elseif refbonds isa AbstractArray{Bond}
        custombonds = true
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
    active = custombonds
    toggle = Makie.Toggle(fig; active, buttoncolor, framecolor_inactive, framecolor_active)
    color = set_alpha.(bond_colors, 0.25)
    draw_bonds(; ax, obs=toggle.active, ionradius, exchange_mag, cryst, interactions, bonds=refbonds, refbonds, color)
    toggle_grid[toggle_cnt+=1, 1:2] = [toggle, Makie.Label(fig, "Reference bonds"; fontsize, halign=:left)]
    
    # Toggle on/off bonds within each class
    for (i, (b, bond_color)) in enumerate(zip(refbonds, bond_colors))
        active = (i == 1)
        framecolor_active = set_alpha(bond_color, 0.7)
        framecolor_inactive = set_alpha(bond_color, 0.15)
        toggle = Makie.Toggle(fig; active, buttoncolor, framecolor_inactive, framecolor_active)
        bonds = propagate_reference_bond_for_cell(cryst, b)
        refbonds = fill(b, length(bonds))
        color = fill(set_alpha(bond_color, 0.25), length(bonds))
        draw_bonds(; ax, obs=toggle.active, ionradius, exchange_mag, cryst, interactions, bonds, refbonds, color)
        bondstr = "Bond($(b.i), $(b.j), $(b.n))"
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
