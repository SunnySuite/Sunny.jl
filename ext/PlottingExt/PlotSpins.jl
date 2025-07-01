# Wrapper over `FigureLike` to support both `show` and `notify`.
struct NotifiableFigure
    notifier :: Makie.Observable{Nothing}
    figure :: Makie.FigureLike
end
Base.showable(mime::MIME, fig::NotifiableFigure) = showable(mime, fig.figure)
Base.show(io::IO, ::MIME"text/plain", fig::NotifiableFigure) = print(io, "(Notifiable) " * repr(fig.figure))
Base.show(io::IO, m::MIME, fig::NotifiableFigure) = show(io, m, fig.figure)
Base.display(fig::NotifiableFigure; kwargs...) = display(fig.figure; kwargs...)
Base.notify(fig::NotifiableFigure) = notify(fig.notifier)
Makie.record(func, nf::NotifiableFigure, path, iter; kwargs...) = Makie.record(func, nf.figure, path, iter; kwargs...)

"""
    plot_spins(sys::System; arrowscale=1.0, color=:red, colorfn=nothing,
               colormap=:viridis, colorrange=nothing, show_cell=true, orthographic=false,
               ghost_radius=0, ndims=3, compass=true)

Plot the spin configuration defined by `sys`. Optional parameters are:

  - `arrowscale`: Scale all arrows by dimensionless factor.
  - `color`: Arrow colors. May be symbolic or numeric. If scalar, will be shared
    among all sites.
  - `colorfn`: Function that dynamically maps from a site index to a numeric
    color value. Useful for animations.
  - `colormap`, `colorrange`: Used to populate colors from numbers following
    Makie conventions.
  - `show_cell`: Show original crystallographic unit cell.
  - `orthographic`: Use orthographic camera perspective.
  - `ghost_radius`: Show periodic images up to a given distance (length units).
  - `ndims`: Spatial dimensions of system (1, 2, or 3).
  - `compass`: If true, draw Cartesian axes in bottom left.

Calling `notify` on the return value will animate the figure.
"""
function Sunny.plot_spins(sys::System; size=(768, 512), compass=true, kwargs...)
    fig = Makie.Figure(; size)
    ax = Makie.LScene(fig[1, 1]; show_axis=false)
    notifier = Makie.Observable(nothing)
    Sunny.plot_spins!(ax, sys; notifier, kwargs...)
    compass && add_cartesian_compass(fig, ax)
    return NotifiableFigure(notifier, fig)
end

"""
    plot_spins!(ax, sys::System; opts...)

Mutating variant of [`plot_spins`](@ref) that allows drawing into a single panel
of a Makie figure.

# Example

```julia
fig = Figure()
plot_spins!(fig[1, 1], sys1)
plot_spins!(fig[2, 1], sys2)
display(fig)
```
"""
function Sunny.plot_spins!(ax, sys::System; notifier=Makie.Observable(nothing), arrowscale=1.0, stemcolor=:lightgray, color=:red,
                     colorfn=nothing, colormap=:viridis, colorrange=nothing, show_cell=true, orthographic=false,
                     ghost_radius=0, ndims=3, dims=nothing)
    isnothing(dims) || error("Use notation `ndims=$dims` instead of `dims=$dims`")

    if ndims == 2
        sys.dims[3] == 1 || error("System not two-dimensional in (a₁, a₂)")
    elseif ndims == 1
        sys.dims[[2,3]] == [1,1] || error("System not one-dimensional in (a₁)")
    end

    # Show bounding box of magnetic supercell in gray (this needs to come first
    # to set a scale for the scene in case there is only one atom).
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.dims))
    Makie.linesegments!(ax, cell_wireframe(supervecs, ndims); color=:gray, linewidth=1.5)

    # Infer characteristic length scale between sites
    ℓ0 = characteristic_length_between_atoms(orig_crystal(sys))

    # Quantum spin-s averaged over sites. Will be used to normalize dipoles.
    N0 = norm(sys.Ns) / sqrt(length(sys.Ns))
    s0 = (N0 - 1) / 2

    # Parameters defining arrow shape
    a0 = arrowscale * ℓ0
    markersize = 0.1*a0
    shaftradius = 0.06a0
    tipradius = 0.2a0
    tiplength = 0.4a0
    lengthscale = 0.6a0

    # Positions in fractional coordinates of supercell vectors
    rs = [supervecs \ global_position(sys, site) for site in eachsite(sys)]

    for isghost in (false, true)
        if isghost
            alpha = 0.08
            (idxs, offsets) = Sunny.all_offsets_within_distance(supervecs, rs, cell_center(ndims); max_dist=ghost_radius, nonzeropart=true)
        else
            alpha = 1.0
            idxs = eachindex(rs)
            offsets = [zero(Vec3) for _ in idxs]
        end

        # Every call to RGBf constructor allocates, so pre-calculate color
        # arrays to speed animations
        cmap_with_alpha = set_alpha.(Makie.to_colormap(colormap), Ref(alpha))
        numeric_colors = zeros(size(sys.dipoles))
        rgba_colors = zeros(Makie.RGBAf, size(sys.dipoles))

        if isnothing(colorfn)
            # In this case, we can precompute the fixed `rgba_colors` array
            # according to `color`
            if color isa AbstractArray
                @assert length(color) == length(sys.dipoles)
                if eltype(color) <: Number
                    dyncolorrange = @something colorrange extrema(color)
                    numbers_to_colors!(rgba_colors, color, cmap_with_alpha, dyncolorrange)
                else
                    map!(rgba_colors, color) do c
                        set_alpha(Makie.to_color(c), alpha)
                    end
                end
            else
                c = set_alpha(Makie.to_color(color), alpha)
                fill!(rgba_colors, c)
            end
        end

        # These observables will be reanimated upon calling `notify(notifier)`.
        pts = Makie.Observable(Makie.Point3f[])
        offset_pts = Makie.Observable(Makie.Point3f[])
        vecs = Makie.Observable(Makie.Vec3f[])
        tipcolor = Makie.Observable(Makie.RGBAf[])

        Makie.on(notifier, update=true) do _
            empty!.((vecs[], offset_pts[], tipcolor[]))

            # Dynamically adapt `rgba_colors` according to `colorfn`
            if !isnothing(colorfn)
                numeric_colors .= colorfn.(CartesianIndices(sys.dipoles))
                dyncolorrange = @something colorrange extrema(numeric_colors)
                numbers_to_colors!(rgba_colors, numeric_colors, cmap_with_alpha, dyncolorrange)
            end

            for (site, n) in zip(idxs, offsets)
                (; offset, vec) = scaled_dipole_to_arrow_geometry(sys.dipoles[site]/s0, lengthscale, tiplength)
                pt = supervecs * (rs[site] + n)
                push!(pts[], Makie.Point3f(pt))
                push!(offset_pts[], pt + Makie.Point3f(offset))
                push!(vecs[], Makie.Vec3f(vec))
                push!(tipcolor[], rgba_colors[site])
            end
            # Trigger Makie redraw
            notify.((offset_pts, vecs, tipcolor))
            # isnothing(color) || notify(arrowcolor)
        end

        # Draw arrows
        shaftcolor = (stemcolor, alpha)
        if !isempty(offset_pts[])
            Makie.arrows3d!(ax, offset_pts, vecs; align=0, markerscale=1, minshaftlength=0, tipradius,
                            shaftradius, tiplength, tipcolor, shaftcolor, diffuse=1.15, transparency=isghost)
        end

        # Small sphere inside arrow to mark atom position
        Makie.meshscatter!(ax, pts; markersize, color=shaftcolor, diffuse=1.15, transparency=isghost)
    end

    # Bounding box of original crystal unit cell in teal
    if show_cell
        Makie.linesegments!(ax, cell_wireframe(orig_crystal(sys).latvecs, ndims); color=:teal, linewidth=1.5)
        pos = [(3/4)*Makie.Point3f(p) for p in eachcol(orig_crystal(sys).latvecs)[1:ndims]]
        text = [Makie.rich("a", Makie.subscript(repr(i))) for i in 1:ndims]
        Makie.text!(ax, pos; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
                    glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-1f0)
    end

    orient_camera!(ax, supervecs; ghost_radius, ℓ0, orthographic, ndims)

    return ax
end
