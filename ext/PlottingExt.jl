module PlottingExt

using Sunny
import Sunny: Vec3, orig_crystal, natoms # Private functions

using LinearAlgebra
import Makie
import Colors: RGB

const seaborn_bright = [
    RGB{Float64}(0.00784313725490196,0.24313725490196078,1.0),
    RGB{Float64}(1.0,0.48627450980392156,0.0),
    RGB{Float64}(0.10196078431372549,0.788235294117647,0.2196078431372549),
    RGB{Float64}(0.9098039215686274,0.0,0.043137254901960784),
    RGB{Float64}(0.5450980392156862,0.16862745098039217,0.8862745098039215),
    RGB{Float64}(0.6235294117647059,0.2823529411764706,0.0),
    RGB{Float64}(0.9450980392156862,0.2980392156862745,0.7568627450980392),
    RGB{Float64}(0.6392156862745098,0.6392156862745098,0.6392156862745098),
    RGB{Float64}(1.0,0.7686274509803922,0.0),
    RGB{Float64}(0.0,0.8431372549019608,1.0),
]

const seaborn_muted = [
    RGB{Float64}(0.2823529411764706,0.47058823529411764,0.8156862745098039),
    RGB{Float64}(0.9333333333333333,0.5215686274509804,0.2901960784313726),
    RGB{Float64}(0.41568627450980394,0.8,0.39215686274509803),
    RGB{Float64}(0.8392156862745098,0.37254901960784315,0.37254901960784315),
    RGB{Float64}(0.5843137254901961,0.4235294117647059,0.7058823529411765),
    RGB{Float64}(0.5490196078431373,0.3803921568627451,0.23529411764705882),
    RGB{Float64}(0.8627450980392157,0.49411764705882355,0.7529411764705882),
    RGB{Float64}(0.4745098039215686,0.4745098039215686,0.4745098039215686),
    RGB{Float64}(0.8352941176470589,0.7333333333333333,0.403921568627451),
    RGB{Float64}(0.5098039215686274,0.7764705882352941,0.8862745098039215),
]

getindex_cyclic(a, i) = a[mod1(i, length(a))] 

function fill_colors(c::AbstractArray, sz)
    size(c) == sz || error("Colors array must have size $sz.")
    if eltype(c) <: Number
        c = numbers_to_colors(c)
    end
    return c
end
fill_colors(c, sz) = fill(c, sz)


# TODO: See if we can reuse Makie internal functions, e.g. `numbers_to_colors`,
# https://github.com/MakieOrg/Makie.jl/blob/ac02141c4c87dbf71d06b301f6dc18f5719e6d05/src/colorsampler.jl#L154-L177
"""
    numbers_to_colors(x; colorrange=nothing, colormap=:viridis)

Converts each number in `x` to a color according to a given `colormap`. The data
in `x` will be scaled according to `colorrange`, which defaults to the min and
max values of `x`. This function mirrors the color conventions of Makie [1].

[1] https://docs.makie.org/stable/documentation/colors/.
"""
function numbers_to_colors(xs; colorrange=nothing, colormap=:viridis)
    (cmin, cmax) = @something colorrange extrema(xs)
    colors = Makie.cgrad(colormap)
    return [colors[(x - cmin) / (cmax - cmin)] for x in xs]
end


function orient_camera!(ax, latvecs, dist; orthographic=false)
    # Set up camera
    lookat = latvecs * Vec3(0.5, 0.5, 0.5)

    # This hack seems needed to get good feeling Makie camera rotations
    eyeposition = lookat - dist * (diagm([1, 1, 0]) * latvecs * Vec3(1, 1, 0))
    upvector = Vec3(0, 0, 1)
    @assert norm((lookat-eyeposition) ⋅ upvector) < 1e-12

    projectiontype = orthographic ? Makie.Orthographic : Makie.Perspective
    # TODO: investigate cam3d_cad! later
    Makie.cam3d!(ax.scene; lookat, eyeposition, upvector, projectiontype)

    # TODO: Make this work. It gives weird results when applied to the reshaped
    # system.
    #=
    lookat = latvecs * Vec3(0.5, 0.5, 0.5)
    eyeposition = lookat - dist * normalize(latvecs * Vec3(-1, -1, -1))

    z = Vec3(0, 0, 1)
    delta = normalize(lookat - eyeposition)
    upvector = normalize(z - delta * (delta ⋅ z))
    @assert norm(delta ⋅ upvector) < 1e-12

    Makie.cam3d_cad!(ax.scene; eyeposition, lookat, upvector)
    =#
end


function cell_wireframe(latvecs)
    vecs = Makie.Point3f0.(eachcol(latvecs))
    ret = Tuple{Makie.Point3f0, Makie.Point3f0}[]

    origin = zero(Makie.Point3f0)

    for j in 0:1, k in 0:1
        shift = j*vecs[2]+k*vecs[3]
        push!(ret, (origin+shift, vecs[1]+shift))
    end
    for i in 0:1, k in 0:1
        shift = i*vecs[1]+k*vecs[3]
        push!(ret, (origin+shift, vecs[2]+shift))
    end
    for i in 0:1, j in 0:1
        shift = i*vecs[1]+j*vecs[2]
        push!(ret, (origin+shift, vecs[3]+shift))
    end
    return ret
end


# TODO: We could rewrite `all_bonds_for_atom` to use this function.
function all_images_within_distance(latvecs, rs, r0s; min_dist=0, max_dist, include_zeros=false)
    # box_lengths[i] represents the perpendicular distance between two parallel
    # boundary planes spanned by lattice vectors a_j and a_k (where indices j
    # and k differ from i)
    box_lengths = [a⋅b/norm(b) for (a,b) = zip(eachcol(latvecs), eachrow(inv(latvecs)))]
    n_max = round.(Int, max_dist ./ box_lengths, RoundUp)

    # optionally initialize to include all (0,0,0) images
    images = include_zeros ? [[zero(Vec3)] for _ in rs] : [Vec3[] for _ in rs]

    # loop over all center points
    for r0 in r0s
        # loop over each atom in primary cell or system
        for (ns, r) in zip(images, rs)
            # loop over image cells or systems
            for n1 in -n_max[1]:n_max[1]
                for n2 in -n_max[2]:n_max[2]
                    for n3 in -n_max[3]:n_max[3]
                        # track list of periodic offsets where the atom image is
                        # within distance bounds
                        n = Vec3(n1, n2, n3)
                        dist = norm(latvecs * (r + n - r0))
                        if min_dist <= dist <= max_dist && !(n in ns)
                            push!(ns, n)
                        end
                    end
                end
            end
        end
    end

    return images
end


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


"""
    plot_spins(sys::System; arrowscale=1.0, color=:red, show_cell=true,
               orthographic=false, ghost_radius=0)

Plot the spin configuration defined by `sys`. Optional parameters include:

  - `arrowscale`: Scale all arrows by dimensionless factor.
  - `color`: Arrow color. May be a numeric value per site in system.
  - `show_cell`: Show original crystallographic unit cell.
  - `orthographic`: Use camera with orthographic projection.
  - `ghost_radius`: Show translucent periodic images up to a given distance
    (length units).
"""
function Sunny.plot_spins(sys::System; arrowscale=1.0, stemcolor=:lightgray, color=:red, show_cell=true,
                          orthographic=false, ghost_radius=0, resolution=(768, 512), show_axis=false, rescale=1.0)
    fig = Makie.Figure(; resolution)
    ax = Makie.LScene(fig[1, 1]; show_axis)
    plot_spins!(ax, sys; arrowscale, stemcolor, color, show_cell, orthographic, ghost_radius, rescale)
    return fig
end

"""
    plot_spins!(ax, sys::System; arrowscale=1.0, color=:red,
                show_cell=true, orthographic=false, ghost_radius=0)

Like [`plot_spins`](@ref) but will draw into the given Makie Axis, `ax`.
"""
function plot_spins!(ax, sys::System; arrowscale=1.0, stemcolor=:lightgray, color=:red, show_cell=true,
                     orthographic=false, ghost_radius=0, rescale=1.0)
    # TODO: Why can't this move to the bottom?
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.latsize))
    orient_camera!(ax, supervecs, 1.0; orthographic)

    ### Plot spins ###

    # Show bounding box of magnetic supercell in gray (this needs to come first
    # to set a scale for the scene in case there is only one atom).
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.latsize))
    Makie.linesegments!(ax, cell_wireframe(supervecs); color=:gray, linewidth=rescale*1.5)

    # Bounding box of original crystal unit cell in teal
    if show_cell
        Makie.linesegments!(ax, cell_wireframe(orig_crystal(sys).latvecs); color=:teal, linewidth=rescale*1.5)
    end

    # Infer characteristic length scale between sites
    ℓ0 = characteristic_length_between_atoms(orig_crystal(sys))

    # Quantum spin-S, averaged over all sites. Will be used to normalize
    # dipoles.
    S0 = (sum(sys.Ns)/length(sys.Ns) - 1) / 2

    # Parameters defining arrow shape
    a0 = arrowscale * ℓ0
    arrowsize = 0.4a0
    linewidth = 0.12a0
    lengthscale = 0.6a0
    markersize = 0.8linewidth
    arrow_fractional_shift = 0.6
   
    # Make sure colors are indexable by site
    color0 = fill_colors(color, size(sys.dipoles))

    # Find all sites within max_dist of the system center
    rs = [supervecs \ global_position(sys, site) for site in eachsite(sys)]
    r0 = [0.5, 0.5, 0.5]
    images = all_images_within_distance(supervecs, rs, [r0]; max_dist=ghost_radius, include_zeros=true)

    # Require separate drawing calls with `transparency=true` for ghost sites
    for (isghost, alpha) in ((false, 1.0), (true, 0.08))
        pts = Makie.Point3f0[]
        vecs = Makie.Vec3f0[]
        arrowcolor = Tuple{eltype(color0), Float64}[]
        for site in eachindex(images)
            vec = (lengthscale / S0) * sys.dipoles[site]
            # Loop over all periodic images of site within radius
            for n in images[site]
                # If drawing ghosts, require !iszero(n), and vice versa
                iszero(n) == isghost && continue
                pt = supervecs * (rs[site] + n)
                push!(pts, Makie.Point3f0(pt))
                push!(vecs, Makie.Vec3f0(vec))
                push!(arrowcolor, (color0[site], alpha))
            end
        end
        shifted_pts = pts - arrow_fractional_shift * vecs
        linecolor = @something (stemcolor, alpha) arrowcolor
        Makie.arrows!(ax, shifted_pts, vecs; arrowsize, linewidth, linecolor, arrowcolor, transparency=isghost)

        # Small sphere inside arrow to mark atom position
        Makie.meshscatter!(ax, pts; markersize, color=linecolor, transparency=isghost)
    end

    if show_cell
        # Labels for lattice vectors. This needs to come last for
        # `overdraw=true` to work.
        pos = [(3/4)*Makie.Point3f0(p) for p in collect(eachcol(orig_crystal(sys).latvecs))]
        text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:3]
        Makie.text!(ax, pos; text, color=:black, fontsize=rescale*20, font=:bold, glowwidth=4.0,
                    glowcolor=(:white, 0.6), align=(:center, :center), overdraw=true)
    end

    return ax
end

"""
    view_crystal(crystal::Crystal, max_dist::Real; show_axis=true, orthographic=false)

An interactive crystal viewer, with bonds up to `max_dist`.
"""
function Sunny.view_crystal(cryst::Crystal, max_dist; show_axis=true, orthographic=false,
                            spherescale=0.2, resolution=(768, 512), rescale=1.0)
    fig = Makie.Figure(; resolution)
    ax = Makie.LScene(fig[1, 1], show_axis=false)

    # TODO: Why can't this move to the bottom?
    orient_camera!(ax, cryst.latvecs, 1; orthographic)

    # Show cell volume and label lattice vectors (this needs to come first to
    # set a scale for the scene in case there is only one atom).
    Makie.linesegments!(ax, cell_wireframe(cryst.latvecs); color=:teal, linewidth=rescale*1.5, inspectable=false)

    # Draw Cartesian axes
    axissize = (1/3)*minimum(norm.(eachcol(cryst.latvecs)))
    axes = Makie.arrows!(ax, Makie.Point3f0.(fill([0,0,0], 3)), axissize*Makie.Point3f0.([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                         arrowsize=0.5axissize, linewidth=0.15axissize, color=[:red, :orange, :yellow], inspectable=false, visible=show_axis)

    # Map atom classes to indices that run from 1..nclasses
    unique_classes = unique(cryst.classes)
    class_indices = [findfirst(==(c), unique_classes) for c in cryst.classes]

    # Show atoms
    ℓ0 = characteristic_length_between_atoms(cryst)
    markersize = spherescale * ℓ0
    max_dist = max(max_dist, ℓ0 + 1e-6)
    images = all_images_within_distance(cryst.latvecs, cryst.positions, cryst.positions; max_dist, include_zeros=true)
    atom_labels = nothing
    for (isghost, alpha) in ((true, 0.15), (false, 1.0))
        pts = Makie.Point3f0[]
        color = RGB{Float64}[]
        for i in eachindex(images), n in images[i]
            # If drawing ghosts, require !iszero(n), and vice versa
            iszero(n) == isghost && continue
            push!(pts, cryst.latvecs * (cryst.positions[i] + n))
            push!(color, getindex_cyclic(seaborn_muted, class_indices[i]))
        end
        Makie.meshscatter!(ax, pts; markersize, color, alpha, inspectable=false, transparency=isghost)

        # Atom indices
        if !isghost
            text = repr.(eachindex(pts))
            atom_labels = Makie.text!(ax, pts; text, color=:white, fontsize=rescale*14, align=(:center, :center),
                                      overdraw=true, visible=true)
        end
    end

    # Get up to 10 reference bonds, without self bonds
    refbonds = filter(reference_bonds(cryst, max_dist)) do b
        return !(b.i == b.j && iszero(b.n))
    end
    refbonds = first(refbonds, 10)

    function all_segments_for_bond(b, color, visible)
        # Prune bonds
        bonds = filter(Sunny.all_symmetry_related_bonds(cryst, b)) do b
            if iszero(collect(b.n))
                # Bonds within the unit cell must not be self bonds, and must not be
                # duplicated.
                return b.i != b.j && Sunny.bond_parity(b)
            else
                # Bonds between two unit cells can always be include
                return true
            end            
        end

        # Capture results of `print_bond` for each b′
        bond_labels = map(bonds) do b′
            basis = Sunny.basis_for_exchange_on_bond(cryst, b′; b_ref=b)
            basis_strs = Sunny.coupling_basis_strings(zip('A':'Z', basis); digits=12, atol=1e-12)
            rows = map(row -> join(row, "   "), eachrow(basis_strs))
            J_matrix_str = "[" * join(rows, ";\n") * "]"
            return "$b′\n$J_matrix_str"
        end

        # Map each bond to line segments in global coordinates
        segments = map(bonds) do b
            (; ri, rj) = Sunny.BondPos(cryst, b)
            Makie.Point3f0.(Ref(cryst.latvecs) .* (ri, rj))
        end
        
        # TODO: Report bug of ÷2 indexing
        inspector_label(plot, index, position) = bond_labels[index ÷ 2]
        s = Makie.linesegments!(ax, segments; color, linewidth=rescale*3,
                                inspectable=true, inspector_label, visible)
        return [s]
    end

    layout = Makie.GridLayout(; tellheight=false, valign=:top)

    # Toggle on/off Cartesian axes
    fontsize = rescale*16
    toggle_cnt = 0
    axes_toggle = Makie.Toggle(fig; active=axes.visible[], buttoncolor=:gray)
    Makie.connect!(axes.visible, axes_toggle.active)
    axes_labels = Makie.GridLayout()
    axes_labels[1, 1] = Makie.Label(fig, "Show"; fontsize)
    axes_labels[1, 2] = Makie.Label(fig, "x"; color=RGB(0.90, 0.0, 0.0), font=:bold, fontsize)
    axes_labels[1, 3] = Makie.Label(fig, "y"; color=RGB(0.90, 0.5, 0.0), font=:bold, fontsize)
    axes_labels[1, 4] = Makie.Label(fig, "z"; color=RGB(0.90, 0.85, 0.0), font=:bold, fontsize)
    layout[toggle_cnt+=1, 1:2] = [axes_toggle, axes_labels]

    # Toggle on/off atom indices
    atom_labels_toggle = Makie.Toggle(fig; active=true, buttoncolor=:gray)
    Makie.connect!(atom_labels.visible, atom_labels_toggle.active)
    layout[toggle_cnt+=1, 1:2] = [atom_labels_toggle, Makie.Label(fig, "Show atom indices"; fontsize)]
    
    # Toggle on/off bonds
    for (i, b) in enumerate(refbonds)
        color = getindex_cyclic(seaborn_bright, i)
        active = (i == 1)
        toggle = Makie.Toggle(fig; active, framecolor_active=color, buttoncolor=:gray)
        observables = all_segments_for_bond(b, color, active)
        for o in observables
            Makie.connect!(o.visible, toggle.active)
        end
        # Equivalent:
        # Makie.on(x -> segments.visible[] = x, toggle.active; update=true)

        layout[toggle_cnt+=1, 1:2] = [toggle, Makie.Label(fig, repr(b); fontsize)]
    end

    fig[1, 2] = layout


    # Label lattice vectors. Putting this last helps with visibility (Makie
    # v0.19)
    pos = [(3/4)*Makie.Point3f0(p) for p in collect(eachcol(cryst.latvecs))]
    text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:3]
    Makie.text!(ax, pos; text, color=:black, fontsize=rescale*20, font=:bold, glowwidth=4.0,
                glowcolor=(:white, 0.6), align=(:center, :center), overdraw=true)

    # Add inspector for pop-up information. Putting this last helps with
    # visibility (Makie v0.19)
    Makie.DataInspector(ax; fontsize)

    return fig
end


function draw_level!(ax,n_level,level,center,radius,dir,z; arrows = true, linewidth, lengthscale, arrowsize)
    if level == n_level || level == 1
        top_level = level == n_level
        col = map(x -> Makie.Colors.HSVA(rad2deg(angle(x[level])),1,1,abs2(x[level])),z)
        if arrows
          Makie.arrows!(ax,center,(top_level ? radius : -radius) .* dir,color = col; linewidth, arrowsize)
        else
          Makie.scatter!(ax,center .+ (top_level ? radius : -radius) .* dir,color = col)
        end
    else
        theta = range(0,2π,length=16)
        for i in eachindex(center)
            normal_dir = norm(dir[i] × [0,0,1]) < 1e-4 ? [1,0,0] : [0,0,1]

            codir1 = normalize(dir[i] × normal_dir)
            codir2 = normalize(codir1 × dir[i])
            l = (n_level - 1)/2
            m = (level - 1) - l
            phi = acos(m/l)
            pts = Vector{Makie.Point3f}(undef,length(theta))
            for j = eachindex(theta)
                pts[j] = center[i] .+ sin(phi) .* radius .* (cos(theta[j]) .* codir1 .+ sin(theta[j]) .* codir2) .+ radius .* (m/l) .* dir[i]
            end
            Makie.lines!(pts,color = Makie.Colors.HSVA(rad2deg(angle(z[i][level])),1,1,abs2(z[i][level])); linewidth)
        end
    end
end

function plot_coherents(sys::System{N};scale = 1., quantization_axis = nothing, use_arrows = true, resolution=(768, 512)) where N

    ℓ0 = characteristic_length_between_atoms(orig_crystal(sys))

    # Parameters defining arrow shape
    a0 = scale * ℓ0
    radius = 0.4a0
    arrowsize = 0.4a0
    linewidth = 0.12a0
    lengthscale = 0.6a0
    markersize = 0.52a0
    #arrow_fractional_shift = 0.6


    n_level = length(sys.coherents[1])

    fig = Makie.Figure(; resolution)
    ax = Makie.LScene(fig[1, 1])

    # Hack to get the `one_dim_chain.jl` examples working. TODO: make
    # `orient_camera!` work better.
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.latsize))
    lookat = sum(eachcol(supervecs)/2)
    eyeposition = lookat - [0, 1, 0]
    Makie.cam3d_cad!(ax.scene; lookat, eyeposition, projectiontype=Makie.Orthographic)
    
    centers = [Makie.Point3f(Sunny.global_position(sys,site)) for site in eachsite(sys)][:]
    Makie.scatter!(ax,centers,color = :black,marker='x';markersize)

    dir = zeros(Makie.Point3f,length(sys.coherents))
    opacity = sys.coherents[:]
    for (i,site) in enumerate(eachsite(sys))
      z = sys.coherents[site]
      v = if isnothing(quantization_axis)
        normalize(Sunny.expected_spin(z))
      else
        quantization_axis
      end
      S = spin_operators(sys,site[4])
      spin_operator = S[1] .* v[1] .+ S[2] .* v[2] .+ S[3] .* v[3]
      basis_rotation = eigvecs(spin_operator;sortby = λ -> -real(λ))
      dir[i] = Makie.Point3f(v...)
      opacity[i] = basis_rotation' * z
    end

    for level = 1:n_level
        draw_level!(ax,n_level,level,centers,radius,dir,opacity;linewidth,lengthscale,arrowsize, arrows = use_arrows)
    end

    fig
end


function scatter_bin_centers(params;axes)
    labels = ["Qx [r.l.u]","Qy [r.l.u.]","Qz [r.l.u.]","E [meV]"]
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1],xlabel = labels[axes[1]], ylabel = labels[axes[2]])
    scatter_bin_centers!(ax,params;axes)
    fig
end

function scatter_bin_centers!(ax,params;axes)
    bcs = axes_bincenters(params)
    xs = Vector{Float64}(undef,0)
    ys = Vector{Float64}(undef,0)
    for xx = bcs[axes[1]], yy = bcs[axes[2]]
        push!(xs,xx)
        push!(ys,yy)
    end
    Makie.scatter!(ax,xs,ys,marker='x',markersize=10,color = :black)
end


# The purpose of __init__() below is to make all the internal functions of
# PlottingExt accessible to developers of Sunny.
#
# The standard and recommended use of Julia package extensions is to add methods
# to existing functions.
# https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions).
# For public exports, we create a function stub in Sunny.jl using the syntax
# `function f end`. Then the implementation is provided in this extension module
# as `function Sunny.f() ... end`.
#
# For non-public functions, however, it is undesirable fill Sunny.jl with stubs
# that will be irrelevant to most users. Access to such internal functions will
# instead be provided through the global variable `Sunny.Plotting`, which is set
# below. Note that `@__MODULE__` references the current extension module, here
# `PlottingExt`.
#
# Without the global variable `Sunny.Plotting`, one would need to use something
# like `Base.get_extension(Sunny, :PlottingExt)` to find the extension module.
function __init__()
    Sunny.Plotting = @__MODULE__
end

end
