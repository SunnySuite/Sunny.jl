"""
    generate_mantid_script_from_binning_parameters(params::BinningParameters)

Generate a Mantid script which bins data according to the 
given [`BinningParameters`](@ref).

!!! warning "Units"

    Take care to ensure the units are correct (R.L.U. or absolute).
    You may want to call `Sunny.bin_rlu_as_absolute_units!` or `Sunny.bin_absolute_units_as_rlu!` first.
"""
function generate_mantid_script_from_binning_parameters(params)
    covectorsK = params.covectors # Please call rlu_to_absolute_units! first if needed
    #function bin_string(k)
        #if params.numbins[k] == 1
            #return "$(params.binsstart[k]),$(params.binend[k])"
        #else
            #return "$(params.binsstart[k]),$(params.binend[k])"
        #end
    #end
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

"""
    params, signal = load_nxs(filename)

Given the name of a Mantid-exported `MDHistoWorkspace` file, load the [`BinningParameters`](@ref) and the signal from that file.
"""
function load_nxs(filename)
    JLD2.jldopen(filename,"r") do file
        read_covectors_from_axes_labels = false
        covectors = Matrix{Float64}(undef,3,3)
        try
          try
            w_matrix = file["MDHistoWorkspace"]["experiment0"]["logs"]["W_MATRIX"]["value"]
            covectors .= reshape(w_matrix,3,3) # No transpose because stored as column
          catch e
            printstyled("Warning",color=:yellow)
            print(": failed to load W_MATRIX from Mantid file $filename due to:\n")
            println(e)
            println("Falling back to reading transform_from_orig")

            coordinate_system = file["MDHistoWorkspace"]["coordinate_system"][1]

            # Permalink to where this enum is defined:
            # https://github.com/mantidproject/mantid/blob/057df5b2de1c15b819c7dd06e50bdbf5461b09c6/Framework/Kernel/inc/MantidKernel/SpecialCoordinateSystem.h#L14
            system_name = ["None", "QLab", "QSample", "HKL"][coordinate_system + 1]

            # The matrix stored in the file is transpose of the actual
            # transform_from_orig matrix
            transform_from_orig = file["MDHistoWorkspace"]["transform_from_orig"]
            transform_from_orig = reshape(transform_from_orig,5,5)
            transform_from_orig = transform_from_orig'
            
            # U: Orthogonal rotation matrix
            # B: inv(lattice_vectors(...)) matrix, as in Sunny
            # The matrix stored in the file is transpose of U * B
            ub_matrix = file["MDHistoWorkspace"]["experiment0"]["sample"]["oriented_lattice"]["orientation_matrix"]
            ub_matrix = ub_matrix'

            # Following: https://docs.mantidproject.org/nightly/concepts/Lattice.html
            # It can be verified that the matrix G^* = (ub_matrix' * ub_matrix)
            # is equal to B' * B, where B = inv(lattice_vectors(...)), and the diagonal
            # entries of the inverse of this matrix are the lattice parameters squared
            #
            #abcMagnitude = sqrt.(diag(inv(ub_matrix' * ub_matrix)))
            #println("[a,b,c] = ",abcMagnitude)

            # This is how you extract the covectors from `transform_from_orig` and `ub_matrix`
            # TODO: Parse this from the `long_name` of the data_dims instead
            covectors .= 2π .* transform_from_orig[1:3,1:3] * ub_matrix
          end
        catch e
          printstyled("Warning",color=:yellow)
          print(": failed to read histogram axes from Mantid file $filename due to:\n")
          println(e)
          println("Defaulting to low-accuracy reading of axes labels!")
          read_covectors_from_axes_labels = true
        end

        signal = file["MDHistoWorkspace"]["data"]["signal"]

        axes = Dict(JLD2.load_attributes(file,"MDHistoWorkspace/data/signal"))[:axes]

        # Axes are just stored backwards in Mantid .nxs for some reason
        axes_names = reverse(split(axes,":"))

        data_dims = Vector{Vector{Float64}}(undef,length(axes_names))
        binwidth = Vector{Float64}(undef,4)
        binstart = Vector{Float64}(undef,4)
        binend = Vector{Float64}(undef,4)
        std = x -> sqrt(sum((x .- sum(x) ./ length(x)).^2))
        for (i,name) in enumerate(axes_names)
            if read_covectors_from_axes_labels
                long_name = Dict(JLD2.load_attributes(file,"MDHistoWorkspace/data/$name"))[:long_name]

                if i <= 3 # This is long_name contains a covector
                    covectors[i,:] .= parse_long_name(long_name)
                    covectors[i,:] .= transpose(pinv(covectors[i,:]))
                end
            end

            data_dims[i] = file["MDHistoWorkspace"]["data"][name]
            binwidth[i] = sum(diff(data_dims[i])) / length(diff(data_dims[i]))
            if std(diff(data_dims[i])) > 1e-4 * binwidth[i]
              printstyled("Warning possible non-uniform binning: mean width = $(binwidth[i]),  std width = $(std(diff(data_dims[i])))"; color = :yellow)
            end

            binstart[i] = minimum(data_dims[i])
            binend[i] = maximum(data_dims[i])
        end

        covectors4D = [covectors [0;0;0]; [0 0 0] 1]
        return BinningParameters(binstart,binend,binwidth,covectors4D), signal
    end
end

# Parse the "[0.5H,0.3H,0.1H]" type of Mantid string describing
# a histogram axis
function parse_long_name(long_name)
  # You're welcome
  found = match(r"\[(|0|(?:-?[0-9]?(?:\.[0-9]+)?))[HKL]?,(|0|(?:-?[0-9]?(?:\.[0-9]+)?))[HKL]?,(|0|(?:-?[0-9]?(?:\.[0-9]+)?))[HKL]?\]",long_name)
  if isnothing(found)
    return nothing
  end
  return map(x -> isempty(x) ? 1. : x == "-" ? -1. : parse(Float64,x),found)
end

function quick_view_nxs(filename,keep_ax)
    integration_axes = setdiff(1:4,keep_ax)
    params, signal = load_nxs(filename)
    integrate_axes!(params,axes = integration_axes)
    int_signal = dropdims(sum(signal,dims = integration_axes);dims = Tuple(integration_axes))
    bcs = axes_bincenters(params)
    (bcs[keep_ax[1]],bcs[keep_ax[2]],int_signal)
end


