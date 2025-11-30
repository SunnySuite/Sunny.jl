function get_unit_energy(units, into)
    if isnothing(units)
        isnothing(into) || error("`into` parameter requires `units` parameter")
        return (1.0, "Energy")
    else
        sym = @something into units.energy
        str = Sunny.unit_strs[sym]
        return (getproperty(units, sym), "Energy ($str)")
    end
end

function colorrange_from_data(; data, saturation, sensitivity, allpositive)
    max_data = ndims(data) == 1 ? data : vec(maximum(data; dims=1))
    min_data = ndims(data) == 1 ? data : vec(minimum(data; dims=1))
    cmax = Statistics.quantile(filter(!isnan, max_data), saturation)
    cmin = Statistics.quantile(filter(!isnan, min_data), 1 - saturation)

    # The returned pair (lo, hi) should be strictly ordered, lo < hi, for use in
    # Makie.Colorbar
    if allpositive
        # Intensities are supposed to be non-negative in this branch, but might
        # not be due to any of: User error, very small round-off around negative
        # zero, or Gibbs ringing in a KPM calculation.
        @assert 0 <= sensitivity < 1
        if cmax <= 0
            return (0, 1e-15)
        else
            (sensitivity, 1) .* cmax
        end
    else
        cscale = max(abs(cmax), abs(cmin))
        if iszero(cscale)
            return (-1e-15, 1e-15)
        else
            return (-cscale, cscale)
        end
    end
end

# Makie.heatmap can fail on large grids. A workaround is to wrap the data in a
# Makie.Resampler if any of the array dimensions is "large"
# https://github.com/MakieOrg/Makie.jl/issues/4950. For OpenGL backends, the
# cutoff 2,000 should (hopefully) be conservative for most systems. The Cairo
# backend only errors for sizes exceeding 32,745.
function is_texture_too_big(data)
    return any(>=(2000), size(data))
end

function heatmap_aux!(ax, x, y, data; opts...)
    if is_texture_too_big(data)
        pl = Makie.heatmap!(ax, x, y, Makie.Resampler(data); opts...)
        # Avoid transparency bug https://github.com/MakieOrg/Makie.jl/issues/5215
        pl.plots[1].visible[]
        pl.plots[1].visible[] = false
        return pl
    else
        return Makie.heatmap!(ax, x, y, data; opts...)
    end
end

"""
    plot_intensities!(panel, res; opts...)

Mutating variant of [`plot_intensities`](@ref) that allows drawing into a single
panel of a Makie figure. Returns a `Makie.Axis` object to allow for further
customization of labels, ticks, legends, etc.

# Example

```julia
fig = Figure()
plot_intensities!(fig[1, 1], res1; title="Panel 1")
plot_intensities!(fig[2, 1], res2; title="Panel 2")
display(fig)
```
"""
function Sunny.plot_intensities!(panel, res::Sunny.BandIntensities{Float64}; colormap=nothing,
                                 colorrange=nothing, saturation=0.9, sensitivity=0.0025,
                                 allpositive=true, interpolate=true, units=nothing, into=nothing,
                                 ylims=nothing, fwhm=nothing, title="", axis=NamedTuple())
    unit_energy, ylabel = get_unit_energy(units, into)
    axis = (; title, axis...)
 
    if res.qpts isa Sunny.QPath 
        mindisp, maxdisp = extrema(res.disp)
        ylims = @something ylims (min(0, mindisp), 1.1*maxdisp) ./ unit_energy
        ebounds = ylims .* unit_energy
        energies = range(ebounds[1], ebounds[2], 512)
        fwhm = @something fwhm 0.02*(ebounds[2]-ebounds[1])
        σ = fwhm/2√(2log(2))
        kernel = Sunny.Broadening(x -> exp(-x^2/2σ^2)) # Gaussian without normalization
        (; data) = Sunny.broaden(res; energies, kernel)

        colorrange_suggest = colorrange_from_data(; data, saturation, sensitivity, allpositive)
        colormap = @something colormap (allpositive ? reverse_thermal_fade : blue_white_red_fade)
        colorrange = @something colorrange colorrange_suggest

        xticklabelrotation = maximum(length.(res.qpts.xticks[2])) > 3 ? π/6 : 0.0
        ax = Makie.Axis(panel; xlabel="Momentum (r.l.u.)", ylabel, res.qpts.xticks, xticklabelrotation, limits=(nothing, ylims), axis...)
        heatmap_aux!(ax, (1, size(data, 2)), ylims, data'; colorrange, colormap, interpolate)
        for i in axes(res.disp, 1)
            Makie.lines!(ax, res.disp[i,:]/unit_energy; color=(:lightskyblue3, 0.75))
        end
        return ax
    else
        error("Cannot plot type $(typeof(res.qpts))")
    end
end

function grid_aspect_ratio(cryst::Crystal, grid::Sunny.QGrid{2})
    # Aspect ratio for global distances
    Δq_global = cryst.recipvecs * (grid.qs[end] - grid.qs[begin])
    e1, e2 = normalize.(Ref(cryst.recipvecs) .* grid.axes)
    abs(dot(e1, e2)) < 1e-12 || error("Cannot yet plot non-orthogonal grid")
    return (Δq_global ⋅ e1) / (Δq_global ⋅ e2)
end

function suggest_labels_for_grid(grid::Sunny.QGrid{N}) where N
    (; axes, offset) = grid

    varidxs = [findmax(abs.(a))[2] for a in axes]
    if varidxs[2] == varidxs[1]
        varidxs[2] = mod1(varidxs[2] + 1, 3)
    end
    varstrs = ("H", "K", "L")

    labels = map(axes, varstrs[varidxs]) do axis, c
        elems = map(axis) do x
            if abs(x) < 1e-12
                "0"
            else
                Sunny.coefficient_to_math_string(x)*c
            end
        end
        return "[" * join(elems, ", ") * "]"
    end

    if norm(offset) > 1e-12
        label1 = labels[begin] * " + " * Sunny.vec3_to_string(offset)
        labels = (label1, labels[2:end]...)
    end

    return labels
end


function Sunny.plot_intensities!(panel, res::Sunny.Intensities{Float64}; colormap=nothing, colorrange=nothing,
                                 saturation=0.9, allpositive=true, interpolate=false, units=nothing, into=nothing,
                                 ylims=nothing, title="", axis=NamedTuple())
    axis = (; title, axis...)
    (; crystal, qpts, data, energies) = res

    colorrange_suggest = colorrange_from_data(; data, saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)
    colorrange = @something colorrange colorrange_suggest

    if qpts isa Sunny.QPath
        unit_energy, ylabel = get_unit_energy(units, into)
        xticklabelrotation = maximum(length.(qpts.xticks[2])) > 3 ? π/6 : 0.0
        ax = Makie.Axis(panel[1, 1]; xlabel="Momentum (r.l.u.)", ylabel, qpts.xticks, xticklabelrotation, limits=(nothing, ylims), axis...)
        hm = heatmap_aux!(ax, (1, size(data, 2)), extrema(energies)./unit_energy, data'; colormap, colorrange, interpolate)
        Makie.Colorbar(panel[1, 2], hm)
        return ax
    elseif qpts isa Sunny.QGrid{2}
        if isone(length(energies))
            aspect = grid_aspect_ratio(crystal, qpts)
            xlabel, ylabel = suggest_labels_for_grid(qpts)
            (xbounds, ybounds) = zip(qpts.coefs_lo, qpts.coefs_hi)
            ax = Makie.Axis(panel[1, 1]; xlabel, ylabel, aspect, axis...)
            hm = heatmap_aux!(ax, xbounds, ybounds, dropdims(data; dims=1); colormap, colorrange, interpolate)
            Makie.Colorbar(panel[1, 2], hm)
            return ax
        else
            error("Cannot yet plot $(typeof(res))")
        end
    else
        error("Cannot yet plot $(typeof(res))")
    end
end

function Sunny.plot_intensities!(panel, res::Sunny.StaticIntensities{Float64}; colormap=nothing, colorrange=nothing,
                                 saturation=0.9, allpositive=true, interpolate=false, units=nothing, into=nothing,
                                 ylims=nothing, title="", axis=NamedTuple())
    axis = (; title, axis...)
    (; crystal, qpts, data) = res

    colorrange_suggest = colorrange_from_data(; data, saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)

    if qpts isa Sunny.QPath
        ylims = @something ylims colorrange (colorrange_suggest .* 1.1)
        xticklabelrotation = maximum(length.(qpts.xticks[2])) > 3 ? π/6 : 0.0
        ax = Makie.Axis(panel; xlabel="Momentum (r.l.u.)", ylabel="Intensity", qpts.xticks, xticklabelrotation, limits=(nothing, ylims), axis...)
        Makie.lines!(ax, data)
        return ax
    elseif qpts isa Sunny.QGrid{2}
        colorrange = @something colorrange colorrange_suggest
        aspect = grid_aspect_ratio(crystal, qpts)
        xlabel, ylabel = suggest_labels_for_grid(qpts)
        (xbounds, ybounds) = zip(qpts.coefs_lo, qpts.coefs_hi)
        ax = Makie.Axis(panel[1, 1]; xlabel, ylabel, aspect, axis...)
        hm = heatmap_aux!(ax, xbounds, ybounds, data; colormap, colorrange, interpolate)
        Makie.Colorbar(panel[1, 2], hm)
        return ax
    else
        error("Cannot yet plot $(typeof(res))")
    end
end

function Sunny.plot_intensities!(panel, res::Sunny.PowderIntensities{Float64}; colormap=nothing, colorrange=nothing,
                                 saturation=0.9, allpositive=true, interpolate=false, units=nothing, into=nothing,
                                 ylims=nothing, title="", axis=NamedTuple())
    axis = (; title, axis...)
    unit_energy, ylabel = get_unit_energy(units, into)
    xlabel = isnothing(units) ? "Momentum " : "Momentum ($(Sunny.unit_strs[units.length])⁻¹)" 
 
    colorrange_suggest = colorrange_from_data(; res.data, saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)
    colorrange = @something colorrange colorrange_suggest

    ax = Makie.Axis(panel[1, 1]; xlabel, ylabel, limits=(nothing, ylims), axis...)
    hm = heatmap_aux!(ax, extrema(res.radii), extrema(res.energies)./unit_energy, res.data'; colormap, colorrange, interpolate)
    Makie.Colorbar(panel[1, 2], hm)
    return ax
end

function Sunny.plot_intensities!(panel, res::Sunny.PowderStaticIntensities{Float64}; colorrange=nothing,
                                 saturation=1.0, allpositive=true, interpolate=false, units=nothing,
                                 ylims=nothing, title="", axis=NamedTuple())
    axis = (; title, axis...)
    xlabel = isnothing(units) ? "Momentum " : "Momentum ($(Sunny.unit_strs[units.length])⁻¹)" 
    ylabel = "Intensity"
 
    colorrange_suggest = colorrange_from_data(; res.data, saturation, sensitivity=0, allpositive)
    ylims = @something ylims colorrange (colorrange_suggest .* 1.1)
    ax = Makie.Axis(panel[1, 1]; xlabel, ylabel, limits=(nothing, ylims), axis...)
    Makie.lines!(ax, res.radii, res.data)
    return ax
end

# Axes will currently be labeled as a linear combination of crystal lattice
# vectors. See https://github.com/SunnySuite/Sunny.jl/pull/310 for details.
function Sunny.plot_intensities!(panel, res::Sunny.BinnedIntensities{Float64}; colormap=nothing, colorrange=nothing,
                                 saturation=0.9, allpositive=true, interpolate=false, units=nothing, into=nothing,
                                 title="", axis=NamedTuple(), divide_counts=true)
    axis = (; title, axis...)
    unit_energy, elabel = get_unit_energy(units, into)

    data = divide_counts ? (res.data ./ res.counts) : res.data

    colorrange_suggest = colorrange_from_data(; data, saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)
    colorrange = @something colorrange colorrange_suggest

    # Low-dimension cases
    n_dims_resolved = count(res.params.numbins .!= 1)

    if n_dims_resolved == 0
        # No resolved data: Just display the one value!
        ax = Makie.Axis(panel[1,1]; axis...)
        text = Sunny.number_to_simple_string(data[1]; digits=4)
        Makie.text!(ax, 0, 0; text)
        return ax
    elseif n_dims_resolved == 1
        # Only resolved on one axis!
        x_axis = findfirst(res.params.numbins .!= 1)
        xlabel = (x_axis == 4) ? elabel : Sunny.covector_name(res.params.covectors[x_axis, :])
        ax = Makie.Axis(panel[1,1]; xlabel, ylabel="Integrated Intensity", axis...)
        bcs = Sunny.axes_bincenters(res.params)
        bcs[4] /= unit_energy
        Makie.barplot!(ax, bcs[x_axis], data[:]; colormap, colorrange)
        return ax
    elseif n_dims_resolved == 2
        x_axis = findfirst(res.params.numbins .!= 1)
        y_axis = x_axis + findfirst(res.params.numbins[x_axis+1:end] .!= 1)
        xlabel = Sunny.covector_name(res.params.covectors[x_axis, :])
        ylabel = (y_axis == 4) ? elabel : Sunny.covector_name(res.params.covectors[y_axis, :])
        ax = Makie.Axis(panel[1,1]; xlabel, ylabel, axis...)
        bcs = Sunny.axes_bincenters(res.params)
        bcs[4] /= unit_energy
        data = reshape(data, size(data, x_axis), size(data, y_axis))
        hm = heatmap_aux!(ax, bcs[x_axis], bcs[y_axis], data; colormap, colorrange, interpolate)
        Makie.Colorbar(panel[1, 2], hm)
        return ax
    elseif n_dims_resolved > 2
        error("Higher-dimensional binned data visualization not yet implemented!")
    end
end

#=
  * `axis`: An additional collection of named arguments that will be passed to
    the `Makie.Axis` constructor. This allows to override the axis `xlabel`,
    `ylabel`, `xticks`, etc. See [Makie
    documentation](https://docs.makie.org/release/reference/blocks/axis#attributes)
    for details.
=#
"""
    plot_intensities(res; colormap=nothing, colorrange=nothing, allpositive=true,
                     saturation=0.9, units=nothing, into=nothing, ylims=nothing,
                     title="", interpolate=false)

Plot the intensities data in `res`, which may be obtained from
[`intensities_bands`](@ref), [`intensities`](@ref),
[`intensities_static`](@ref), [`powder_average`](@ref), or similar.

See also the mutating variant [`plot_intensities!`](@ref), which allows to draw
multiple panels into a `Makie.Figure`, and to further customize the
`Makie.Axis`.

Keyword arguments:

  * `colormap`: Color palette for plotting broadened intensities. See Makie docs
    for allowed values.
  * `colorrange`: A lower and upper bound on intensities. For heatmaps, these
    bounds define the intensity values that saturate the colormap.
  * `allpositive`: Should intensities be all positive, apart from numerical
    error? If true, the default colors will clip below zero intensity. If false,
    the default colors will be symmetric about zero intensity.
  * `saturation`: If `colorrange` is not explicitly set, this dimensionless
    parameter defines the upper saturated intensity value as a quantile of
    maximum intensities taken over wavevectors.
  * `interpolate`: Should 2D grid data be plotted with an interpolated heatmap?
    False by default, except when plotting discrete bands.
  * `units`: A [`Units`](@ref) instance for labeling axes and performing
    conversions.
  * `into`: A symbol for conversion into a new base energy unit (`:meV`, `:K`,
    etc.)
  * `ylims`: Limits of the y-axis.
  * `title`: An optional title for the plot.

When plotting discrete [`intensities_bands`](@ref) data, an additional argument
is:

  * `fwhm`: Apply Gaussian broadening with a specific full-width at
    half-maximum.

!!! tip "Exporting to PDF"

    To export the figure to a PDF file, use the CairoMakie backend and call

    ```julia
        fig = plot_intensities(..., interpolate=true)
        save("myfile.pdf", fig)
    ```

    The choice `interpolate=true` significantly reduces file sizes for plots
    with 2D heatmap data.
"""
function Sunny.plot_intensities(res::Sunny.AbstractIntensities; opts...)
    fig = Makie.Figure()
    Sunny.plot_intensities!(fig[1,1], res; opts...)
    return fig
end
