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
