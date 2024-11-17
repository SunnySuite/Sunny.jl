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
    cmax = Statistics.quantile(filter(!isnan, vec(maximum(data; dims=1))), saturation)
    cmin = Statistics.quantile(filter(!isnan, vec(minimum(data; dims=1))), 1 - saturation)

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

"""
    plot_intensities!(panel, res::AbstractIntensities; opts...)

Mutating variant of [`plot_intensities`](@ref) that allows drawing into a single
panel of a Makie figure.

# Example

```julia
fig = Figure()
plot_intensities!(fig[1, 1], res1; title="Panel 1")
plot_intensities!(fig[2, 1], res2; title="Panel 2")
display(fig)
```
"""
function Sunny.plot_intensities!(panel, res::Sunny.BandIntensities{Float64}; colormap=nothing, colorrange=nothing,
                                 saturation=0.9, sensitivity=0.0025, allpositive=true, units=nothing, into=nothing,
                                 fwhm=nothing, ylims=nothing, title="", axisopts=NamedTuple())
    unit_energy, ylabel = get_unit_energy(units, into)
    axisopts = (; title, axisopts...)
 
    if res.qpts isa Sunny.QPath 
        mindisp, maxdisp = extrema(res.disp)
        ylims = @something ylims (min(0, mindisp), 1.2*maxdisp)
        energies = range(ylims[1], ylims[2], 512)
        fwhm = @something fwhm 0.02*(ylims[2]-ylims[1])
        σ = fwhm/2√(2log(2))
        kernel = Sunny.Broadening(x -> exp(-x^2/2σ^2)) # Gaussian without normalization
        broadened = Sunny.broaden(res; energies, kernel)

        colorrange_suggest = colorrange_from_data(; broadened.data, saturation, sensitivity, allpositive)
        colormap = @something colormap (allpositive ? Makie.Reverse(:thermal) : :bwr)
        colorrange = @something colorrange colorrange_suggest

        xticklabelrotation = maximum(length.(res.qpts.xticks[2])) > 3 ? π/6 : 0.0
        ax = Makie.Axis(panel; xlabel="Momentum (r.l.u.)", ylabel, limits=(nothing, ylims ./ unit_energy), res.qpts.xticks, xticklabelrotation, axisopts...)
        Makie.heatmap!(ax, axes(res.data, 2), collect(energies / unit_energy), broadened.data'; colorrange, colormap, lowclip=:white)
        for i in axes(res.disp, 1)
            Makie.lines!(ax, res.disp[i,:] / unit_energy; color=:lightskyblue3)
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
        label1 = labels[begin] * " + " * Sunny.fractional_vec3_to_string(offset)
        labels = (label1, labels[2:end]...)
    end

    return labels
end


function Sunny.plot_intensities!(panel, res::Sunny.Intensities{Float64}; colormap=nothing, colorrange=nothing, saturation=0.9, allpositive=true, units=nothing, into=nothing, title="", axisopts=NamedTuple())
    axisopts = (; title, axisopts...)
    (; crystal, qpts, data, energies) = res

    colorrange_suggest = colorrange_from_data(; data, saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)
    colorrange = @something colorrange colorrange_suggest
    
    if qpts isa Sunny.QPath
        unit_energy, ylabel = get_unit_energy(units, into)
        xticklabelrotation = maximum(length.(qpts.xticks[2])) > 3 ? π/6 : 0.0
        ax = Makie.Axis(panel[1, 1]; xlabel="Momentum (r.l.u.)", ylabel, qpts.xticks, xticklabelrotation, axisopts...)
        hm = Makie.heatmap!(ax, axes(data, 2), collect(energies/unit_energy), data'; colormap, colorrange)
        Makie.Colorbar(panel[1, 2], hm)
        return ax
    elseif qpts isa Sunny.QGrid{2}
        if isone(length(energies))
            aspect = grid_aspect_ratio(crystal, qpts)
            xlabel, ylabel = suggest_labels_for_grid(qpts)
            (xs, ys) = map(range, qpts.coefs_lo, qpts.coefs_hi, size(qpts.qs))
            ax = Makie.Axis(panel[1, 1]; xlabel, ylabel, aspect, axisopts...)
            hm = Makie.heatmap!(ax, xs, ys, dropdims(data; dims=1); colormap, colorrange)
            Makie.Colorbar(panel[1, 2], hm)
            return ax
        else
            error("Cannot yet plot $(typeof(res))")
        end
    else
        error("Cannot yet plot $(typeof(res))")
    end
end

function Sunny.plot_intensities!(panel, res::Sunny.StaticIntensities{Float64}; colormap=nothing, colorrange=nothing, saturation=0.9, allpositive=true, units=nothing, into=nothing, title="", axisopts=NamedTuple())
    axisopts = (; title, axisopts...)
    (; crystal, qpts, data) = res

    colorrange_suggest = colorrange_from_data(; data=reshape(data, 1, size(data)...), saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)
    colorrange = @something colorrange colorrange_suggest

    if qpts isa Sunny.QGrid{2}
        aspect = grid_aspect_ratio(crystal, qpts)
        xlabel, ylabel = suggest_labels_for_grid(qpts)
        (xs, ys) = map(range, qpts.coefs_lo, qpts.coefs_hi, size(qpts.qs))
        ax = Makie.Axis(panel[1, 1]; xlabel, ylabel, aspect, axisopts...)
        hm = Makie.heatmap!(ax, xs, ys, data; colormap, colorrange)
        Makie.Colorbar(panel[1, 2], hm)
        return ax
    elseif qpts isa Sunny.QPath
        xticklabelrotation = maximum(length.(qpts.xticks[2])) > 3 ? π/6 : 0.0
        ax = Makie.Axis(panel; xlabel="Momentum (r.l.u.)", ylabel="Intensity", qpts.xticks, xticklabelrotation, axisopts...)
        Makie.lines!(ax, data)
        Makie.ylims!(ax, colorrange)
        return ax
    else
        error("Cannot yet plot $(typeof(res))")
    end
end

function Sunny.plot_intensities!(panel, res::Sunny.PowderIntensities{Float64}; colormap=nothing, colorrange=nothing, saturation=0.9, allpositive=true, units=nothing, into=nothing, title="", axisopts=NamedTuple())
    axisopts = (; title, axisopts...)
    unit_energy, ylabel = get_unit_energy(units, into)
    xlabel = isnothing(units) ? "Momentum " : "Momentum ($(Sunny.unit_strs[units.length])⁻¹)" 
 
    colorrange_suggest = colorrange_from_data(; res.data, saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)
    colorrange = @something colorrange colorrange_suggest

    ax = Makie.Axis(panel[1, 1]; xlabel, ylabel, axisopts...)
    hm = Makie.heatmap!(ax, res.radii, collect(res.energies/unit_energy), res.data'; colormap, colorrange)
    Makie.Colorbar(panel[1, 2], hm)
    return ax
end

# Axes will currently be labeled as a linear combination of crystal lattice
# vectors. See https://github.com/SunnySuite/Sunny.jl/pull/310 for details.
function Sunny.plot_intensities!(panel, res::Sunny.BinnedIntensities{Float64}; colormap=nothing, colorrange=nothing, saturation=0.9, allpositive=true, units=nothing, into=nothing, title="", axisopts=NamedTuple(), divide_counts=true)
    axisopts = (; title, axisopts...)
    unit_energy, elabel = get_unit_energy(units, into)

    data = divide_counts ? (res.data ./ res.counts) : res.data

    colorrange_suggest = colorrange_from_data(; data, saturation, sensitivity=0, allpositive)
    colormap = @something colormap (allpositive ? :gnuplot2 : :bwr)
    colorrange = @something colorrange colorrange_suggest

    # Low-dimension cases
    n_dims_resolved = count(res.params.numbins .!= 1)

    if n_dims_resolved == 0
        # No resolved data: Just display the one value!
        ax = Makie.Axis(panel[1,1]; axisopts...)
        text = Sunny.number_to_simple_string(data[1]; digits=4)
        Makie.text!(ax, 0, 0; text)
        return ax
    elseif n_dims_resolved == 1
        # Only resolved on one axis!
        x_axis = findfirst(res.params.numbins .!= 1)
        xlabel = (x_axis == 4) ? elabel : Sunny.covector_name(res.params.covectors[x_axis, :])
        ax = Makie.Axis(panel[1,1]; xlabel, ylabel="Integrated Intensity", axisopts...)
        bcs = Sunny.axes_bincenters(res.params)
        bcs[4] /= unit_energy
        Makie.barplot!(ax, bcs[x_axis], data[:]; colormap, colorrange)
        return ax
    elseif n_dims_resolved == 2
        x_axis = findfirst(res.params.numbins .!= 1)
        y_axis = x_axis + findfirst(res.params.numbins[x_axis+1:end] .!= 1)
        xlabel = Sunny.covector_name(res.params.covectors[x_axis, :])
        ylabel = (y_axis == 4) ? elabel : Sunny.covector_name(res.params.covectors[y_axis, :])
        ax = Makie.Axis(panel[1,1]; xlabel, ylabel, axisopts...)
        bcs = Sunny.axes_bincenters(res.params)
        bcs[4] /= unit_energy
        data = reshape(data, size(data, x_axis), size(data, y_axis))
        hm = Makie.heatmap!(ax, bcs[x_axis], bcs[y_axis], data; colormap, colorrange)
        Makie.Colorbar(panel[1, 2], hm)
        return ax
    elseif n_dims_resolved > 2
        error("Higher-dimensional binned data visualization not yet implemented!")
    end
end

#=
  * `axisopts`: An additional collection of named arguments that will be passed
    to the `Makie.Axis` constructor. This allows to override the axis `xlabel`,
    `ylabel`, `xticks`, etc. See [Makie
    documentation](https://docs.makie.org/release/reference/blocks/axis#attributes)
    for details.
=#
"""
    plot_intensities(res::BandIntensities; colormap=nothing, colorrange=nothing, allpositive=true,
                     saturation=0.9, sensitivity=0.0025, fwhm=nothing, ylims=nothing, units=nothing,
                     into=nothing, title="")
    plot_intensities(res::Intensities; colormap=nothing, colorrange=nothing, allpositive=true, 
                     saturation=0.9, units=nothing, into=nothing, title="")
    plot_intensities(res::StaticIntensities; colormap=nothing, colorrange=nothing, allpositive=true, 
                     saturation=0.9, units=nothing, into=nothing, title="")
    plot_intensities(res::PowderIntensities; colormap=nothing, colorrange=nothing, allpositive=true, 
                     saturation=0.9, units=nothing, into=nothing, title="")

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
  * `sensitivity`: When plotting `BandIntensities`, this defines a lower bound
    on the visible intensities as a fraction of the upper saturated intensity
    value.
  * `fwhm`: When plotting `BandIntensities`, this overrides the full-width at
    half-maximum value used for Gaussian broadening.
  * `ylims`: Limits of the y-axis.
  * `units`: A [`Units`](@ref) instance for labeling axes and performing
    conversions.
  * `into`: A symbol for conversion into a new base energy unit (e.g. `:meV`,
    `:K`, etc.)
  * `title`: An optional title for the plot.
"""
function Sunny.plot_intensities(res::Sunny.AbstractIntensities; opts...)
    fig = Makie.Figure()
    Sunny.plot_intensities!(fig[1,1], res; opts...)
    return fig
end
