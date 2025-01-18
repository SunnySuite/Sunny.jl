function show_qspace_obj(ax, recipvecs, ℓ0, params::Sunny.BinningParameters)
    bin_colors = [:red,:blue,:green]
    line_alpha = 0.3
    bin_line_width = 0.5

    @assert iszero(params.covectors[1:3, 4]) && iszero(params.covectors[4, 1:3])
    bes = Sunny.axes_binedges(params)
    M = recipvecs * inv(params.covectors[1:3, 1:3])

    for dir in 1:3
        ix = [2, 3, 1][dir]
        iy = [3, 1, 2][dir]
        iz = dir

        # The grid of q points making up the lowest side of the histogram
        # along the iz direction
        grid = Vector{Float64}[]
        grid_sparse = Vector{Float64}[]
        for i in 1:length(bes[ix]), j in 1:length(bes[iy])
            is_i_edge = (i == 1 || i == length(bes[ix]))
            is_j_edge = (j == 1 || j == length(bes[iy]))
            grid_point = [bes[ix][i], bes[iy][j], bes[iz][1]][invperm([ix, iy, iz])]
            if is_i_edge && is_j_edge # Corner case; render special for outline
                push!(grid_sparse, grid_point)
                continue
            end
            push!(grid, grid_point)
        end
        offset = [0, 0, bes[iz][end] - bes[iz][1]][invperm([ix, iy, iz])]

        if !isempty(grid)
            segs = map(x -> (Makie.Point3f(M * x), Makie.Point3f(M * (x .+ offset))), grid[:])
            Makie.linesegments!(ax, segs, color=bin_colors[dir], linewidth=bin_line_width, alpha=line_alpha, inspectable=false)
        end

        segs = map(x -> (Makie.Point3f(M * x), Makie.Point3f(M * (x .+ offset))), grid_sparse[:])
        Makie.linesegments!(ax, segs; color=:black, linewidth=2.5, inspectable=false)
    end
end

function show_qspace_obj(ax, recipvecs, ℓ0, qpts::Sunny.AbstractQPoints)
    pts = Makie.Point3f.(Ref(recipvecs) .* vec(qpts.qs))
    Makie.meshscatter!(ax, pts; markersize=ℓ0*0.01, color=:red, inspectable=false)
end


"""
    view_bz(crystal::Crystal, objs...; orthographic=false, compass=true)

**Experimental**

Launches a graphical user interface to visualize reciprocal space with respect
to the [`Crystal`](@ref) unit cell. The first Brilliouin zone for the primitive
lattice is shown as a convex polyhedron. High-symmetry points and paths between
them are inspectable, consistent with the data in
[`print_irreducible_bz_paths`](@ref). Additional `objs` may be passed, e.g., to
visualize custom paths or grids in reciprocal space.

 - `orthographic`: Use orthographic camera perspective.
 - `compass`: If true, draw Cartesian axes in bottom left.
"""
function Sunny.view_bz(cryst::Crystal, objs...; orthographic=false, compass=true, size=(768, 512))
    fig = Makie.Figure(; size)

    # Main scene
    ax = Makie.LScene(fig[1, 1], show_axis=false)

    # Show Cartesian axes, with link to main camera
    if compass
        axcompass = add_cartesian_compass(fig, ax)
    end

    ℓ0 = maximum(norm.(eachcol(cryst.recipvecs)))

    # Set of widgets
    widget_list = Makie.GridLayout(; tellheight=false, valign=:top)
    fig[1, 2] = widget_list
    widget_cnt = 0
    fontsize = 16

    # Controls for camera perspective
    menu = Makie.Menu(fig; options=["Perspective", "Orthographic"], default=(orthographic ? "Orthographic" : "Perspective"), fontsize)
    button = Makie.Button(fig; label="Reset", fontsize)
    Makie.onany(button.clicks, menu.selection; update=true) do _, mselect
        b1, b2, _ = eachcol(cryst.recipvecs)
        lookat = zero(Makie.Point3f)
        camshiftdir = normalize(b1 + b2)
        upvector = normalize(b1 × b2)
        camdist = 1.5 * maximum(norm.(eachcol(cryst.recipvecs)))
        orthographic = (mselect == "Orthographic")
        orient_camera!(ax; lookat, camshiftdir, upvector, camdist, orthographic)
        compass && register_compass_callbacks(axcompass, ax)
    end
    widget_list[widget_cnt+=1, 1] = Makie.hgrid!(menu, button)

    # Show reciprocal vectors
    bs = collect(Makie.Point3f.(eachcol(cryst.recipvecs)))
    b_segments = [(zero(Makie.Point3f), b) for b in bs]
    Makie.linesegments!(ax, b_segments; color=:teal, linewidth=1.5, inspectable=false)

    text = [Makie.rich("b", Makie.subscript(repr(i))) for i in 1:3]
    Makie.text!(ax, bs; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
                glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-1f0)

    # Calculate and show Brillouin zone
    bzcell = Brillouin.cartesianize!(Brillouin.wignerseitz(eachcol(Sunny.prim_recipvecs(cryst))))
    segments = Makie.Point3f[]
    for face in bzcell
        append!(segments, face)
        push!(segments, face[1]) # cycle to first point
        push!(segments, Makie.Point3f(NaN)) # separator between faces
    end
    Makie.lines!(ax, segments; inspectable=false)

    try
        (; paths, points) = Sunny.irreducible_bz_paths(cryst)

        # Show high symmetry paths
        for path in paths
            Makie.lines!(ax, [cryst.recipvecs * points[pt] for pt in path]; inspectable=false, color=:pink, linewidth=4)
        end

        # High symmetry q-points
        names = keys(points)
        pts = [Makie.Point3f(cryst.recipvecs * points[nm]) for nm in names]
        Makie.text!(ax, pts; text=String.(names), fontsize=16, font=:bold, color=:purple, glowwidth=4.0,
                    glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-0.9f0)
    catch e
        if startswith(e.msg, "Triclinic")
            @error """Triclinic lattice angles must currently be all-acute or all-obtuse.
                          See https://github.com/thchr/Brillouin.jl/issues/34
                   """
        else
            rethrow()
        end
    end

    # Show other objects
    for obj in objs
        show_qspace_obj(ax, cryst.recipvecs, ℓ0, obj)
    end

    # See similar command in `view_crystal`
    Makie.DataInspector(ax; indicator_color=:gray, fontsize, font="Deja Vu Sans Mono", depth=(1e4 - 1))

    return fig
end
