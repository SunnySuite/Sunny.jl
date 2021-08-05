"""Plotting functions for lattices and spins on lattices.
"""

import GLMakie

function plot_lattice!(ax, lattice::Lattice{2}; color=:blue, markersize=20, linecolor=:grey, linewidth=1.0, kwargs...)
    # Plot the unit cell mesh
    plot_cells!(ax, lattice; color=linecolor, linewidth=linewidth)

    # Plot markers at each site
    sites = reinterpret(reshape, Float64, collect(lattice))
    xs, ys = vec(sites[1, 1:end, 1:end, 1:end]), vec(sites[2, 1:end, 1:end, 1:end])
    GLMakie.scatter!(ax, xs, ys; color=color, markersize=markersize, show_axis=false, kwargs...)
end

function plot_lattice!(ax, lattice::Lattice{3}; color=:blue, markersize=200, linecolor=:grey, linewidth=1.0, kwargs...)
    # Plot markers at each site
    sites = reinterpret(reshape, Float64, collect(lattice))
    xs = vec(sites[1, 1:end, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end, 1:end])
    zs = vec(sites[3, 1:end, 1:end, 1:end, 1:end])
    GLMakie.scatter!(ax, xs, ys, zs; color=color, markersize=markersize, show_axis=false, kwargs...)

    # For some odd reason, the sites will not appear unless this happens afterwards
    # Plot the unit cell mesh
    plot_cells!(ax, lattice; color=linecolor, linewidth=linewidth)
end

function _setup_2d()
    f = GLMakie.Figure()
    ax = GLMakie.Axis(f[1, 1])
    ax.autolimitaspect = 1
    GLMakie.hidespines!(ax)
    GLMakie.hidedecorations!(ax)
    return f, ax
end

function _setup_3d()
    f = GLMakie.Figure()
    lscene = GLMakie.LScene(f[1, 1], scenekw=(camera=GLMakie.cam3d!, raw=false))
    return f, lscene
end


function plot_lattice(lattice::Lattice{2}; kwargs...)
    f, ax = _setup_2d()
    plot_lattice!(ax, lattice; kwargs...)
    f
end

# 3D is a bit wonky at the moment - Axis3 doesn't seem to work with scatter!,
#   so you need to use LScene
function plot_lattice(lattice::Lattice{3}; kwargs...)
    f, lscene = _setup_3d()
    plot_lattice!(lscene, lattice; kwargs...)
end

# Warning: These functions assume that bondlists have "duplicates". I.e. both
#            Bond(i, j, n) and Bond(j, i, -n) appear.

function plot_bonds(lattice::Lattice{2}, ints::Vector{<:PairInt{2}}; bondwidth=4, kwargs...)
    f, ax = _setup_2d()

    colors = GLMakie.to_colormap(:Dark2_8, 8)
    # Sort interactions so that longer bonds are plotted first
    sort!(ints, by=int->-int.dist)

    toggles = Vector{GLMakie.Toggle}()
    labels = Vector{GLMakie.Label}()

    # Plot the lattice
    plot_lattice!(ax, lattice; kwargs...)
    # Plot the bonds, all relative to a central atom on the first sublattice
    basis_idx = 1
    cent_cell = CartesianIndex(div.(lattice.size .+ 1, 2)...)
    cent_pt = lattice[basis_idx, cent_cell]
    for (n, int) in enumerate(ints)
        xs = Vector{Float64}()
        ys = Vector{Float64}()
        for bond in int.bonds
            if bond.i == basis_idx
                new_cell = modc1(cent_cell + bond.celloffset, lattice.size)
                bond_pt = lattice[bond.j, new_cell]
                push!(xs, cent_pt[1])
                push!(ys, cent_pt[2])
                push!(xs, bond_pt[1])
                push!(ys, bond_pt[2])
            end
        end
        bondlabel = "J" * string(int.dist)
        if !isnothing(int.class)
            bondlabel *= "_" * string(int.class)
        end
        color = colors[mod1(n, 8)]
        seg = GLMakie.linesegments!(xs, ys; linewidth=bondwidth, label=bondlabel, color=color)

        tog = GLMakie.Toggle(f, active=true)
        GLMakie.connect!(seg.visible, tog.active)
        push!(toggles, tog)
        push!(labels, GLMakie.Label(f, bondlabel))
    end
    GLMakie.axislegend()
    f[1, 2] = GLMakie.grid!(hcat(toggles, labels), tellheight=false)
    f
end

function plot_bonds(lattice::Lattice{3}, ints::Vector{<:PairInt{3}}; bondwidth=4, kwargs...)
    f, ax = _setup_3d()

    colors = GLMakie.to_colormap(:Dark2_8, 8)
    # Sort interactions so that longer bonds are plotted first
    sort!(ints, by=int->-int.dist)

    toggles = Vector{GLMakie.Toggle}()
    labels = Vector{GLMakie.Label}()

    # Plot the lattice
    plot_lattice!(ax, lattice; kwargs...)
    # Plot the bonds, all relative to a central atom on the first sublattice
    basis_idx = 1
    cent_cell = CartesianIndex(div.(lattice.size .+ 1, 2)...)
    cent_pt = lattice[basis_idx, cent_cell]
    for (n, int) in enumerate(ints)
        xs = Vector{Float64}()
        ys = Vector{Float64}()
        zs = Vector{Float64}()
        for bond in int.bonds
            if bond.i == basis_idx
                new_cell = modc1(cent_cell + bond.celloffset, lattice.size)
                bond_pt = lattice[bond.j, new_cell]
                push!(xs, cent_pt[1])
                push!(ys, cent_pt[2])
                push!(zs, cent_pt[3])
                push!(xs, bond_pt[1])
                push!(ys, bond_pt[2])
                push!(zs, bond_pt[3])
            end
        end
        bondlabel = "J" * string(int.dist)
        if !isnothing(int.class)
            bondlabel *= "_" * string(int.class)
        end
        color = colors[mod1(n, 8)]
        seg = GLMakie.linesegments!(xs, ys, zs; linewidth=bondwidth, label=bondlabel, color=color)

        tog = GLMakie.Toggle(f, active=true)
        GLMakie.connect!(seg.visible, tog.active)
        push!(toggles, tog)
        push!(labels, GLMakie.Label(f, bondlabel))
    end
    GLMakie.axislegend()
    f[1, 2] = GLMakie.grid!(hcat(toggles, labels), tellheight=false)
    f
end

# TODO: Base.Cartesian could combine these functions
"Plot the outlines of the unit cells of a lattice"
function plot_cells!(ax, lattice::Lattice{2}; color=:grey, linewidth=1.0, kwargs...)
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

    GLMakie.linesegments!(ax, xs, ys; color=color, linewidth=linewidth)
end

"Plot the outlines of the unit cells of a lattice"
function plot_cells!(ax, lattice::Lattice{3}; color=:grey, linewidth=1.0, kwargs...)
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

    GLMakie.linesegments!(ax, xs, ys, zs; color=color, linewidth=linewidth)
end

function plot_spins(sys::SpinSystem{2}; linecolor=:grey, arrowcolor=:red, linewidth=0.1, arrowsize=0.3, kwargs...)
    sites = reinterpret(reshape, Float64, collect(sys.lattice))
    spins = 0.1 .* reinterpret(reshape, Float64, collect(sys.sites))

    xs = vec(sites[1, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end])
    zs = zero(xs)
    us = vec(spins[1, 1:end, 1:end, 1:end])
    vs = vec(spins[2, 1:end, 1:end, 1:end])
    ws = vec(spins[3, 1:end, 1:end, 1:end])

    GLMakie.arrows(
        xs, ys, zs, us, vs, ws;
        linecolor=linecolor, arrowcolor=arrowcolor, linewidth=linewidth, arrowsize=arrowsize,
        show_axis=false, kwargs...    
    )
end

function plot_spins(sys::SpinSystem{3}; linecolor=:grey, arrowcolor=:red, linewidth=0.1, arrowsize=0.3, kwargs...)
    sites = reinterpret(reshape, Float64, collect(sys.lattice))
    spins = 0.2 .* reinterpret(reshape, Float64, collect(sys.sites))

    xs = vec(sites[1, 1:end, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end, 1:end])
    zs = vec(sites[3, 1:end, 1:end, 1:end, 1:end])
    us = vec(spins[1, 1:end, 1:end, 1:end, 1:end])
    vs = vec(spins[2, 1:end, 1:end, 1:end, 1:end])
    ws = vec(spins[3, 1:end, 1:end, 1:end, 1:end])

    GLMakie.arrows(
        xs, ys, zs, us, vs, ws;
        linecolor=linecolor, arrowcolor=arrowcolor, linewidth=linewidth, arrowsize=arrowsize,
        show_axis=false, kwargs...    
    )
end

"Equivalent to above, but different arguments"
function plot_spins(lat::Lattice{3}, spins; linecolor=:grey, arrowcolor=:red, linewidth=0.1, arrowsize=0.3, kwargs...)
    sites = reinterpret(reshape, Float64, collect(lat))
    spins = reinterpret(reshape, Float64, spins.val)

    xs = vec(sites[1, 1:end, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end, 1:end])
    zs = vec(sites[3, 1:end, 1:end, 1:end, 1:end])
    us = vec(spins[1, 1:end, 1:end, 1:end, 1:end])
    vs = vec(spins[2, 1:end, 1:end, 1:end, 1:end])
    ws = vec(spins[3, 1:end, 1:end, 1:end, 1:end])

    GLMakie.arrows(
        xs, ys, zs, us, vs, ws;
        linecolor=linecolor, arrowcolor=arrowcolor, linewidth=linewidth, arrowsize=arrowsize,
        show_axis=false, kwargs...    
    )
end

# No support for higher than 3D visualization, sorry!

"Produce an animation of spin dynamics for the specified length of time"
function anim_integration(
    sys::SpinSystem{2}, fname, steps_per_frame, Δt, nframes;
    linecolor=:grey, arrowcolor=:red, linewidth=0.1, arrowsize=0.2, kwargs...
)
    sites = reinterpret(reshape, Float64, collect(sys.lattice))
    spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
    
    xs = vec(sites[1, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end])
    zs = zero(ys)
    us = GLMakie.Node(vec(spins[1, 1:end, 1:end, 1:end]))
    vs = GLMakie.Node(vec(spins[2, 1:end, 1:end, 1:end]))
    ws = GLMakie.Node(vec(spins[3, 1:end, 1:end, 1:end]))
    fig, ax, plot = GLMakie.arrows(
        xs, ys, zs, us, vs, ws;
        linecolor=linecolor, arrowcolor=arrowcolor, linewidth=linewidth, arrowsize=arrowsize,
        show_axis=false, kwargs...    
    )
    display(fig)

    framerate = 30
    integrator = HeunP(sys)

    GLMakie.record(fig, fname, 1:nframes; framerate=framerate) do frame
        for step in 1:steps_per_frame
            evolve!(integrator, Δt)
        end
        spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
        us[] = vec(spins[1, 1:end, 1:end, 1:end])
        vs[] = vec(spins[2, 1:end, 1:end, 1:end])
        ws[] = vec(spins[3, 1:end, 1:end, 1:end])
    end
end

"Produce an animation of spin dynamics for the specified length of time"
function anim_integration(
    sys::SpinSystem{3}, fname, steps_per_frame, Δt, nframes;
    linecolor=:grey, arrowcolor=:red, linewidth=0.1, arrowsize=0.2, kwargs...
)
    sites = reinterpret(reshape, Float64, collect(sys.lattice))
    spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
    
    xs = vec(sites[1, 1:end, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end, 1:end])
    zs = vec(sites[3, 1:end, 1:end, 1:end, 1:end])
    us = GLMakie.Node(vec(spins[1, 1:end, 1:end, 1:end, 1:end]))
    vs = GLMakie.Node(vec(spins[2, 1:end, 1:end, 1:end, 1:end]))
    ws = GLMakie.Node(vec(spins[3, 1:end, 1:end, 1:end, 1:end]))
    fig, ax, plot = GLMakie.arrows(
        xs, ys, zs, us, vs, ws;
        linecolor=linecolor, arrowcolor=arrowcolor, linewidth=linewidth, arrowsize=arrowsize,
        show_axis=false, kwargs...    
    )
    display(fig)

    framerate = 30
    integrator = HeunP(sys)

    GLMakie.record(fig, fname, 1:nframes; framerate=framerate) do frame
        for step in 1:steps_per_frame
            evolve!(integrator, Δt)
        end
        spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
        us[] = vec(spins[1, 1:end, 1:end, 1:end, 1:end])
        vs[] = vec(spins[2, 1:end, 1:end, 1:end, 1:end])
        ws[] = vec(spins[3, 1:end, 1:end, 1:end, 1:end])
    end
end

"Endless integration in a live window"
function live_integration(
    sys::SpinSystem{2}, steps_per_frame, Δt;
    linecolor=:grey, arrowcolor=:red, linewidth=0.1, arrowsize=0.2, kwargs...
)
    sites = reinterpret(reshape, Float64, collect(sys.lattice))
    spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
    
    xs = vec(sites[1, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end])
    zs = zero(ys)
    us = GLMakie.Node(vec(spins[1, 1:end, 1:end, 1:end]))
    vs = GLMakie.Node(vec(spins[2, 1:end, 1:end, 1:end]))
    ws = GLMakie.Node(vec(spins[3, 1:end, 1:end, 1:end]))
    fig, ax, plot = GLMakie.arrows(
        xs, ys, zs, us, vs, ws;
        linecolor=linecolor, arrowcolor=arrowcolor, linewidth=linewidth, arrowsize=arrowsize,
        show_axis=false, kwargs...    
    )
    display(fig)

    framerate = 30
    integrator = HeunP(sys)

    while true
        for step in 1:steps_per_frame
            evolve!(integrator, Δt)
        end
        spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
        us[] = vec(spins[1, 1:end, 1:end, 1:end])
        vs[] = vec(spins[2, 1:end, 1:end, 1:end])
        ws[] = vec(spins[3, 1:end, 1:end, 1:end])
        sleep(1/framerate)
    end
end

"Endless integration in a live window"
function live_integration(
    sys::SpinSystem{3}, steps_per_frame, Δt;
    linecolor=:grey, arrowcolor=:red, linewidth=0.1, arrowsize=0.2, kwargs...
)
    sites = reinterpret(reshape, Float64, collect(sys.lattice))
    spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
    
    xs = vec(sites[1, 1:end, 1:end, 1:end, 1:end])
    ys = vec(sites[2, 1:end, 1:end, 1:end, 1:end])
    zs = vec(sites[3, 1:end, 1:end, 1:end, 1:end])
    us = GLMakie.Node(vec(spins[1, 1:end, 1:end, 1:end, 1:end]))
    vs = GLMakie.Node(vec(spins[2, 1:end, 1:end, 1:end, 1:end]))
    ws = GLMakie.Node(vec(spins[3, 1:end, 1:end, 1:end, 1:end]))
    fig, ax, plot = GLMakie.arrows(
        xs, ys, zs, us, vs, ws;
        linecolor=linecolor, arrowcolor=arrowcolor, linewidth=linewidth, arrowsize=arrowsize,
        show_axis=false, kwargs...    
    )
    display(fig)

    framerate = 30
    integrator = HeunP(sys)

    while true
        for step in 1:steps_per_frame
            evolve!(integrator, Δt)
        end
        spins = 0.2 .* reinterpret(reshape, Float64, sys.sites)
        us[] = vec(spins[1, 1:end, 1:end, 1:end, 1:end])
        vs[] = vec(spins[2, 1:end, 1:end, 1:end, 1:end])
        ws[] = vec(spins[3, 1:end, 1:end, 1:end, 1:end])
        sleep(1/framerate)
    end
end