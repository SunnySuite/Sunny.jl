using LinearAlgebra
import ColorSchemes, ColorTypes


################################################################################
# Functions for plotting on triangular plaquettes, in particular Berry curvature
################################################################################

function vec3_triple_product(s₁, s₂, s₃)
    s₁, s₂, s₃ = normalize.((s₁, s₂, s₃))
    return s₁ ⋅ (s₂ × s₃)
end

function sun_berry_curvature(z₁, z₂, z₃)
    z₁, z₂, z₃ = normalize.((z₁, z₂, z₃))
    n₁ = z₁ ⋅ z₂
    n₂ = z₂ ⋅ z₃
    n₃ = z₃ ⋅ z₁
    return angle(n₁ * n₂ * n₃)
end


function plaquette_idcs(dims::Tuple{Int,Int,Int})
    dx, dy, dz = dims
    (dz != 1) && println("Warning: Ignoring lattice c vector.")
    Triple = Tuple{Int,Int,Int}
    idcs = Array{Tuple{Triple,Triple,Triple},4}(undef, 2, dx + 1, dy + 1, 1)
    for j ∈ 1:(dy+1)
        for i ∈ 1:(dx+1)
            idcs[1, i, j, 1] = (
                (mod1(i  , dx), mod1(j,   dy), 1),
                (mod1(i+1, dx), mod1(j+1, dy), 1),
                (mod1(i  , dx), mod1(j+1, dy), 1),
            )
            idcs[2, i, j, 1] = (
                (mod1(i,   dx), mod1(j,   dy), 1),
                (mod1(i+1, dx), mod1(j,   dy), 1),
                (mod1(i+1, dx), mod1(j+1, dy), 1),
            )
        end
    end
    return idcs
end

plaquette_idcs(x::Array{T,4}) where T = plaquette_idcs(size(x)[1:3])

function plaquette_map(f::Function, x::Array{T,4}) where T
    
    dims = size(x)
    @assert dims[3] == dims[4] == 1 "Multiple sites and multiple layers are not supported."
    dims = dims[1:3]
    out = zeros(Float64, 2, dims...)
    idcs_all = plaquette_idcs(x)
    for i ∈ CartesianIndices(dims)
        idcs = idcs_all[1, i]
        out[1, i] = f(x[idcs[1]..., 1], x[idcs[2]..., 1], x[idcs[3]..., 1])
        idcs = idcs_all[2, i]
        out[2, i] = f(x[idcs[1]..., 1], x[idcs[2]..., 1], x[idcs[3]..., 1])
    end
    out
end


################################################################################
# Plotting functions
################################################################################
function aspect_ratio(x_panel, y_panel, x_offset, y_offset, numrows, numcols)
    corners = [
        (0, 0),
        numrows * y_panel + (numrows - 1) * y_offset,
        numcols * x_panel + (numcols - 1) * x_offset,
        numcols * x_panel + numrows * y_panel + (numrows - 1) * y_offset + (numcols - 1) * x_offset,
    ]
    xs = [c[1] for c ∈ corners]
    ys = [c[2] for c ∈ corners]
    x1, x2 = minimum(xs), maximum(xs)
    y1, y2 = minimum(ys), maximum(ys)
    return abs((x2 - x1) / (y2 - y1))
end


function plot_chirality(xs, sys;
    colorscheme=ColorSchemes.RdBu, clims=(-0.5, 0.5), offset_spacing=1,
    numcols=nothing, texts=nothing, force_aspect=true, fig_kwargs...
)
    # Consolidate lattice info and panel layout
    numpanels = length(xs)
    isnothing(numcols) && (numcols = numpanels)
    numrows = floor(Int, (numpanels - 1) / numcols) + 1
    vecs = sys.crystal.latvecs
    v₁, v₂ = vecs[:, 1], vecs[:, 2]
    nx, ny = size(xs[1])[1:2]
    v₁, v₂ = Point3f(v₁), Point3f(v₂)
    x, y = [1.0, 0, 0], [0.0, 1, 0]
    x_offset = offset_spacing * (x ⋅ v₁) * x
    y_offset = -offset_spacing * (y ⋅ v₂) * y
    x_panel = nx * v₁
    y_panel = -ny * v₂
    aspect = aspect_ratio(x_panel, y_panel, x_offset, y_offset, numrows, numcols)

    # Set up figure
    fig = Figure(; fig_kwargs...)
    if force_aspect
        ax = Axis(fig[1, 1:length(xs)]; aspect)
    else
        ax = Axis(fig[1, 1:length(xs)])
    end
    hidespines!(ax)
    hidedecorations!(ax)

    # Plot panels
    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁+v₂, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂+v₁]))
    for (i, x) ∈ enumerate(xs)
        r, c = fldmod1(i, numcols)
        v₀ = (c - 1) * (x_panel + x_offset) + (r - 1) * (y_panel + y_offset)

        f = (eltype(x) == Sunny.Vec3) ? vec3_triple_product : sun_berry_curvature
        χ = plaquette_map(f, x)
        pgons = GLMakie.Polygon[]
        colors = ColorTypes.RGB{Float64}[]
        for r ∈ 1:nx
            for c ∈ 1:ny
                base = (r - 1) * v₁ + (c - 1) * v₂ + v₀
                push!(pgons, plaq1(base))
                push!(colors, get(colorscheme, χ[1, r, c, 1, 1], clims))
                push!(pgons, plaq2(base))
                push!(colors, get(colorscheme, χ[2, r, c, 1, 1], clims))
            end
        end
        poly!(ax, pgons; color=colors)
        if !isnothing(texts)
            text!(ax, v₀[1], v₀[2]; text=texts[i])
        end
    end
    return fig
end

function plot_chirality(sys;
    colorscheme=ColorSchemes.RdBu, clims=(-0.5, 0.5), offset=1,
    numcols=nothing, texts=nothing, fig_kwargs...
)
    x = (sys.mode == :dipole) ? sys.dipoles : sys.coherents
    return plot_chirality([x], sys; colorscheme, clims, offset, numcols, texts, fig_kwargs...)
end