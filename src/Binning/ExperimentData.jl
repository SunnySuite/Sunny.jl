#=
"""
    generate_mantid_script_from_binning_parameters(params::BinningParameters)

Generate a Mantid script which bins data according to the given
[`BinningParameters`](@ref).

!!! warning "Units"  
    Take care to ensure the units are correct (R.L.U. or absolute). You may want
    to call `Sunny.bin_rlu_as_absolute_units!` or
    `Sunny.bin_absolute_units_as_rlu!` first.
"""
function generate_mantid_script_from_binning_parameters(params)
    covectorsK = params.covectors # Please call rlu_to_absolute_units! first if needed
    # function bin_string(k)
    #     if params.numbins[k] == 1
    #         return "$(params.binsstart[k]),$(params.binend[k])"
    #     else
    #         return "$(params.binsstart[k]),$(params.binend[k])"
    #     end
    # end

    # FIXME: These covectorsK need to be inverted before they are valid QDimensions
    return """MakeSlice(InputWorkspace="merged_mde_INPUT",
        QDimension0="$(covectorsK[1,1]),$(covectorsK[1,2]),$(covectorsK[1,3])",
        QDimension1="$(covectorsK[2,1]),$(covectorsK[2,2]),$(covectorsK[2,3])",
        QDimension2="$(covectorsK[3,1]),$(covectorsK[3,2]),$(covectorsK[3,3])",
        Dimension0Binning="$(params.binstart[1]),$(params.binwidth[1]),$(params.binend[1])",
        Dimension1Name="DeltaE",
        Dimension1Binning="$(params.binstart[2]),$(params.binwidth[2]),$(params.binend[2])",
        Dimension2Binning="$(params.binstart[3]),$(params.binwidth[3]),$(params.binend[3])",
        Dimension3Binning="$(params.binstart[4]),$(params.binwidth[4]),$(params.binend[4])",
        Dimension3Name="QDimension1",
        Smoothing="0",
        OutputWorkspace="Histogram_OUTPUT")
        """
end
=#

"""
    params, signal = load_nxs(filename; field="signal")

Given the name of a Mantid-exported `MDHistoWorkspace` file, load the
[`BinningParameters`](@ref) and the signal from that file.

To load another field instead of the signal, specify e.g.
`field="errors_squared"`. Typical fields include `errors_squared`, `mask`,
`num_events`, and `signal`.
"""
function load_nxs(filename; field="signal")
    JLD2.jldopen(filename, "r") do file
        # This variable is basically the "Mantid W-Matrix". See discussion on
        # Github: https://github.com/SunnySuite/Sunny.jl/pull/256.
        spatial_covectors = zeros(3, 3)

        signal = file["MDHistoWorkspace"]["data"][field]

        axes = JLD2.load_attributes(file, "MDHistoWorkspace/data/signal")["axes"]

        # Axes are just stored backwards in Mantid .nxs for some reason
        axes_names = reverse(split(axes,":"))

        data_dims = Vector{Vector{Float64}}(undef, length(axes_names))
        binwidth = zeros(4)
        binstart = zeros(4)
        binend = zeros(4)
        covectors = zeros(4, 4)
        spatial_covector_ixs = [0,0,0]
        std = x -> sqrt(sum((x .- sum(x) ./ length(x)).^2))
        for (i, name) in enumerate(axes_names)
            long_name = JLD2.load_attributes(file, "MDHistoWorkspace/data/$name")["long_name"]

            if long_name == "DeltaE"
                covectors[i,:] .= [0,0,0,1] # energy covector
            else # spatial covector case
                ix = findfirst(spatial_covector_ixs .== 0)
                spatial_covector_ixs[ix] = i
                lbl = parse_long_name(long_name)
                spatial_covectors[:,ix] .= lbl
            end

            data_dims[i] = file["MDHistoWorkspace"]["data"][name]
            binwidth[i] = sum(diff(data_dims[i])) / length(diff(data_dims[i]))
            if std(diff(data_dims[i])) > 1e-4 * binwidth[i]
              printstyled("Warning possible non-uniform binning: mean width = $(binwidth[i]),  std width = $(std(diff(data_dims[i])))"; color = :yellow)
            end

            binstart[i] = minimum(data_dims[i])

            # Place end of bin in center of last bin, according to Sunny convention
            binend[i] = maximum(data_dims[i]) - binwidth[i]/2
        end

        covectors[spatial_covector_ixs,1:3] .= inv(spatial_covectors)

        return BinningParameters(binstart, binend, binwidth, covectors), signal
    end
end

function Base.permutedims(params::BinningParameters, perm)
    out = copy(params)
    out.covectors .= params.covectors[perm,:]
    out.binwidth .= params.binwidth[perm]
    out.binstart .= params.binstart[perm]
    out.binend .= params.binend[perm]
    out
end

# Parse "[0.5H,0.3H,0.1H]" into the Julia vector [0.5, 0.3, 0.1]. These strings
# typically arises as the Mantid histogram axis label.
function parse_long_name(long_name)
    found = match(r"\[(|0|(?:-?[0-9]?(?:\.[0-9]+)?))[HKL]?,(|0|(?:-?[0-9]?(?:\.[0-9]+)?))[HKL]?,(|0|(?:-?[0-9]?(?:\.[0-9]+)?))[HKL]?\]", long_name)
    if isnothing(found)
        return nothing
    end
    return map(x -> isempty(x) ? 1. : x == "-" ? -1. : parse(Float64, x), found)
end

function quick_view_nxs(filename, keep_ax)
    integration_axes = setdiff(1:4, keep_ax)
    params, signal = load_nxs(filename)
    integrate_axes!(params, axes=integration_axes)
    int_signal = dropdims(sum(signal, dims=integration_axes); dims=Tuple(integration_axes))
    bcs = axes_bincenters(params)
    (bcs[keep_ax[1]], bcs[keep_ax[2]], int_signal)
end
