"""Plotting functions for lattices and spins on lattices.
"""

using GLMakie

function plot_lattice(lattice::Lattice{2}; color=:blue, markersize=20, linecolor=:grey, linewidth=1.0, kwargs...)
    f = Figure()
    ax = Axis(f[1, 1])
    ax.autolimitaspect = 1
    hidespines!(ax)
    hidedecorations!(ax)

    # Plot the unit cell mesh
    plot_cells!(lattice; color=linecolor, linewidth=linewidth)

    # Plot markers at each site
    sites = reinterpret(reshape, Float64, collect(lattice))
    xs, ys = vec(sites[1, 1:end, 1:end, 1:end]), vec(sites[2, 1, 1:end, 1:end])
    scatter!(xs, ys; color=color, markersize=markersize, kwargs...)
    f
end

# 3D is a bit wonky at the moment - Axis3 doesn't seem to work with scatter!
# For now, have to plot sites below the unit cell grid
function plot_lattice(lattice::Lattice{3}; color=:blue, markersize=100, linecolor=:grey, linewidth=1.0, kwargs...)
    # f = Figure()
    # ax = Axis3(f[1, 1], viewmode=:fit)
    # hidespines!(ax)
    # hidedecorations!(ax)

    # Plot markers at each site
    sites = reinterpret(reshape, Float64, collect(lattice))
    xs = vec(sites[1, 1:end, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end, 1:end])
    zs = vec(sites[3, 1:end, 1:end, 1:end, 1:end])
    f = scatter(xs, ys, zs; color=color, markersize=markersize, show_axis=false, kwargs...)

    # For some odd reason, the sites will not appear unless this happens afterwards
    # Plot the unit cell mesh
    plot_cells!(lattice; color=linecolor, linewidth=linewidth)
    f
end

# TODO: Base.Cartesian could combine these functions
"Plot the outlines of the unit cells of a lattice"
function plot_cells!(lattice::Lattice{2}; color=:grey, linewidth=1.0, kwargs...)
    lattice = brav_lattice(lattice)

    xs, ys = Vector{Float64}(), Vector{Float64}()
    nx, ny = lattice.size
    for j in 1:ny
        left_pt, right_pt = lattice[1, 1, j], lattice[1, nx, j]
        push!(xs, left_pt[1])
        push!(xs, right_pt[1])
        push!(ys, left_pt[2])
        push!(ys, right_pt[2])
    end
    for i in 1:nx
        bot_pt, top_pt = lattice[1, i, 1], lattice[1, i, ny]
        push!(xs, bot_pt[1])
        push!(xs, top_pt[1])
        push!(ys, bot_pt[2])
        push!(ys, top_pt[2])
    end

    linesegments!(xs, ys; color=color, linewidth=linewidth)
end

"Plot the outlines of the unit cells of a lattice"
function plot_cells!(lattice::Lattice{3}; color=:grey, linewidth=1.0, kwargs...)
    lattice = brav_lattice(lattice)

    xs, ys, zs = Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    nx, ny, nz = lattice.size
    for j in 1:ny
        for k in 1:nz
            bot_pt, top_pt = lattice[1, 1, j, k], lattice[1, nx, j, k]
            push!(xs, bot_pt[1])
            push!(xs, top_pt[1])
            push!(ys, bot_pt[2])
            push!(ys, top_pt[2])
            push!(zs, bot_pt[3])
            push!(zs, top_pt[3])
        end
        for i in 1:nx
            left_pt, right_pt = lattice[1, i, j, 1], lattice[1, i, j, nz]
            push!(xs, left_pt[1])
            push!(xs, right_pt[1])
            push!(ys, left_pt[2])
            push!(ys, right_pt[2])
            push!(zs, left_pt[3])
            push!(zs, right_pt[3])
        end
    end
    for k in 1:nz
        for i in 1:nx
            left_pt, right_pt = lattice[1, i, 1, k], lattice[1, i, ny, k]
            push!(xs, left_pt[1])
            push!(xs, right_pt[1])
            push!(ys, left_pt[2])
            push!(ys, right_pt[2])
            push!(zs, left_pt[3])
            push!(zs, right_pt[3])
        end
    end

    linesegments!(xs, ys, zs; color=color, linewidth=linewidth)
end


# No support for higher than 3D visualization, sorry!