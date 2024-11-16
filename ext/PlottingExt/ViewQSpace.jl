function Sunny.viz_qqq_path(params::Sunny.BinningParameters; kwargs...)
    f = Makie.Figure()
    ax = Makie.LScene(f[1,1]; show_axis=false)
    Makie.cam3d!(ax.scene; projectiontype=Makie.Orthographic)
    viz_qqq_path!(ax, params; kwargs...)
  
    aabb_lo, aabb_hi = Sunny.binning_parameters_aabb(params)
    lo = min.(0, floor.(Int64, aabb_lo))
    hi = max.(0, ceil.(Int64, aabb_hi))
    Makie.scatter!(ax, map(x -> Makie.Point3f(lo .+ x.I .- 1), CartesianIndices(ntuple(i -> 1 + hi[i] - lo[i], 3)))[:], color=:black)
    global_axes = [(Makie.Point3f(-1,0,0), Makie.Point3f(1,0,0)),
                   (Makie.Point3f(0,-1,0), Makie.Point3f(0,1,0)),
                   (Makie.Point3f(0,0,-1), Makie.Point3f(0,0,1))]
    Makie.linesegments!(ax, global_axes, color=:black)
    Makie.text!(1.1, 0, 0; text="𝐛₁ [R.L.U.]")
    Makie.text!(0, 1.1, 0; text="𝐛₂ [R.L.U.]")
    Makie.text!(0, 0, 1.1; text="𝐛₃ [R.L.U.]")
    display(f)
    ax
end

function viz_qqq_path!(ax, params::Sunny.BinningParameters; line_alpha=0.3, color=nothing, colorrange=nothing, bin_colors=[:red,:blue,:green], bin_line_width=0.5)
    @assert iszero(params.covectors[1:3, 4]) && iszero(params.covectors[4, 1:3])
    bes = Sunny.axes_binedges(params)
    M = inv(params.covectors[1:3, 1:3])
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
            push!(grid,grid_point)
        end
        offset = [0, 0, bes[iz][end] - bes[iz][1]][invperm([ix, iy, iz])]

        if !isempty(grid)
            segs = map(x -> (Makie.Point3f(M * x), Makie.Point3f(M * (x .+ offset))), grid[:])
            Makie.linesegments!(ax, segs, color=bin_colors[dir], linewidth=bin_line_width, alpha=line_alpha)
        end

        segs = map(x -> (Makie.Point3f(M * x), Makie.Point3f(M * (x .+ offset))), grid_sparse[:])
        Makie.linesegments!(ax, segs; color=isnothing(color) ? :black : color, linewidth=2.5, colorrange)
    end
end


function Sunny.view_qspace(cryst::Crystal, orthographic=false, compass=true, size=(768, 512))
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
        b1, b2, _ = eachcol(cryst.recipvecs)
        lookat = zero(Makie.Point3f0)
        camshiftdir = normalize(b1 + b2)
        upvector = normalize(b1 × b2)
        camdist = 1.5 * maximum(norm.(eachcol(cryst.recipvecs)))
        orthographic = (mselect == "Orthographic")
        orient_camera!(ax; lookat, camshiftdir, upvector, camdist, orthographic)
        compass && register_compass_callbacks(axcompass, ax)
    end
    widget_list[widget_cnt+=1, 1] = Makie.hgrid!(menu, button)

    # Show reciprocal vectors
    bs = collect(Makie.Point3f0.(eachcol(cryst.recipvecs)))
    b_segments = [(zero(Makie.Point3f0), b) for b in bs]
    Makie.linesegments!(ax, b_segments; color=:teal, linewidth=1.5, inspectable=false)

    text = [Makie.rich("b", Makie.subscript(repr(i))) for i in 1:3]
    Makie.text!(ax, bs; text, color=:black, fontsize=20, font=:bold, glowwidth=4.0,
                glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-1f0)

    # Map direct basis to ITA standard setting via [aₛ bₛ cₛ] = [a, b, c] P⁻¹.
    std_latvecs = cryst.latvecs / cryst.sg.setting.R
    # Map path and 1st BZ cell back to original Cartesian coordinates.
    bzpath = Brillouin.cartesianize!(Brillouin.irrfbz_path(cryst.sg.number, eachcol(std_latvecs)))
    bzcell = Brillouin.cartesianize!(Brillouin.wignerseitz(Brillouin.basis(bzpath)))

    segments = Makie.Point3f0[]
    for face in bzcell
        append!(segments, face)
        push!(segments, face[1]) # cycle to first point
        push!(segments, Makie.Point3f0(NaN)) # separator between faces
    end
    Makie.lines!(ax, segments; inspectable=false)

    for path in bzpath.paths
        Makie.lines!(ax, [bzpath.points[pt] for pt in path]; inspectable=false, color=:pink, linewidth=4)
    end

    names = keys(bzpath.points)
    pts = Makie.Point3f.(values(bzpath.points))
    Makie.text!(ax, pts; text=String.(names), fontsize=16, font=:bold, color=:purple, glowwidth=4.0,
                glowcolor=(:white, 0.6), align=(:center, :center), depth_shift=-0.9f0)

    # Following `view_crystal` settings
    Makie.DataInspector(ax; indicator_color=:gray, fontsize, font="Deja Vu Sans Mono", depth=(1e4 - 1))

    return fig
end
