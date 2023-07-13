
"""
    BinningParameters(binstart,binend,binwidth;covectors = I(4))
    BinningParameters(binstart,binend;numbins,covectors = I(4))

Describes a 4D parallelepided histogram in a format compatible with experimental Inelasitic Neutron Scattering data.
See [`generate_shiver_script`](@ref) to convert [`BinningParameters`](@ref) to a format understandable by the [Mantid software](https://www.mantidproject.org/), or [`load_nxs_binning_parameters`](@ref) to load a [`BinningParameters`](@ref) from a Mantid `.nxs` file.
 
The coordinates of the histogram axes are specified by multiplication 
of `(q,ω)` with each row of the `covectors` matrix.
Since the `covectors` matrix is the identity matrix by default, the default coordinates
are `(q,ω)` in Reciprocal Lattice Units (R.L.U.). 
To use physical units instead, see [`rlu_to_absolute_units!`](@ref).

The convention for the binning scheme is that:
- The left edge of the first bin starts at `binstart`
- The bin width is `binwidth`
- The last bin contains `binend`
- There are no "partial bins;" the last bin may contain values greater than `binend`. C.f. [`count_bins`](@ref).

A `value` can be binned by computing its bin index:

    coords = covectors * value
    bin_ix = 1 .+ floor.(Int64,(coords .- binstart) ./ binwidth)
"""
mutable struct BinningParameters
    binstart::MVector{4,Float64}
    binend::MVector{4,Float64}
    binwidth::MVector{4,Float64}
    covectors::MMatrix{4,4,Float64}
end

function Base.show(io::IO, ::MIME"text/plain", params::BinningParameters)
    printstyled(io, "Binning Parameters\n"; bold=true, color=:underline)
    nbin = params.numbins
    for k = 1:4
        if nbin[k] == 1
            printstyled(io, "∫ Integrated"; bold=true)
        else
            printstyled(io, @sprintf("⊡ %5d bins",nbin[k]); bold=true)
        end
        @printf(io," from %+.3f to %+.3f along [", params.binstart[k], params.binend[k])
        axes_names = ["x","y","z","E"]
        inMiddle = false
        for j = 1:4
            if params.covectors[k,j] != 0.
                if(inMiddle)
                    print(io," ")
                end
                @printf(io,"%+.2f d%s",params.covectors[k,j],axes_names[j])
                inMiddle = true
            end
        end
        @printf(io,"] (Δ = %.3f)", params.binwidth[k]/norm(params.covectors[k,:]))
        println(io)
    end
end

# Support numbins as a (virtual) property, even though only the binwidth is stored
Base.getproperty(params::BinningParameters, sym::Symbol) = sym == :numbins ? count_bins(params.binstart,params.binend,params.binwidth) : getfield(params,sym)

function Base.setproperty!(params::BinningParameters, sym::Symbol, numbins)
    if sym == :numbins
        binwidth = (params.binend .- params.binstart) ./ numbins

        # *Ensure* that the last bin contains params.binend
        binwidth .+= eps.(binwidth) 

        setfield!(params,:binwidth,binwidth)
    else
        setfield!(params,sym,numbins)
    end
end

"""
    count_bins(binstart,binend,binwidth)

Returns the number of bins in the binning scheme implied by `binstart`, `binend`, and `binwidth`.
To count the bins in a [`BinningParameters`](@ref), use `params.numbins`.

This function defines how partial bins are handled, so it should be used preferentially over
computing the number of bins manually.
"""
count_bins(bin_start,bin_end,bin_width) = ceil.(Int64,(bin_end .- bin_start) ./ bin_width)

function BinningParameters(binstart,binend,binwidth;covectors = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    return BinningParameters(binstart,binend,binwidth,covectors)
end

function BinningParameters(binstart,binend;numbins,kwargs...)
    binwidth = (binend .- binstart) ./ numbins
    binwidth .+= eps.(binwidth)
    return BinningParameters(binstart,binend,binwidth;kwargs...)
end

"""
    integrate_axes!(params::BinningParameters; axes)
Integrate over one or more axes of the histogram by setting the number of bins
in that axis to 1. Examples:

    integrate_axes!(params; axes = [2,3])
    integrate_axes!(params; axes = 2)
"""
function integrate_axes!(params::BinningParameters;axes)
    for k in axes
        nbins = [params.numbins.data...]
        nbins[k] = 1
        params.numbins = SVector{4}(nbins)
    end
    return params
end

# Find an axis-aligned bounding box containing the histogram
function binning_parameters_aabb(params)
    (; binstart, binend, covectors) = params
    bin_edges = [binstart binend]
    this_corner = MVector{4,Float64}(undef)
    q_corners = MMatrix{4,16,Float64}(undef)
    for j = 1:16 # The sixteen corners of a 4-cube
        for k = 1:4 # The four axes
            this_corner[k] = bin_edges[k,1 + (j >> (k-1) & 1)]
        end
        q_corners[:,j] = covectors \ this_corner
    end
    lower_aabb_q = minimum(q_corners,dims=2)[1:3]
    upper_aabb_q = maximum(q_corners,dims=2)[1:3]
    return lower_aabb_q, upper_aabb_q
end

"""
    rlu_to_absolute_units!(params::BinningParameters,[specification of reciprocal lattice vectors])

Converts the [`BinningParameters`](@ref) from Reciprocal Lattice Units R.L.U. to absolute (physical) units.
After calling `rlu_to_absolute_units!`, 

The second argument may be a 3x3 matrix specifying the reciprocal lattice vectors, or any of these objects:
- [`Crystal`](@ref)
- [`System`](@ref)
- [`StructureFactor`](@ref)
- [`SpinWaveTheory`](@ref)
"""
function rlu_to_absolute_units!(params::BinningParameters,recip_vecs)
    covectorsQ = params.covectors

    # covectorsQ * q = covectorsK * k = covectorsK * recip_vecs * q
    covectorsK = covectorsQ * [inv(recip_vecs) [0;0;0]; [0 0 0] 1]
    params.covectors = MMatrix{4,4}(covectorsK)
end
rlu_to_absolute_units!(crystal::Crystal,params::BinningParameters) = rlu_to_absolute_units!(2π*inv(crystal.latvecs)',params)
rlu_to_absolute_units!(sys::System,params::BinningParameters) = rlu_to_absolute_units!(sys.crystal,params)
rlu_to_absolute_units!(sf::StructureFactor,params::BinningParameters) = rlu_to_absolute_units!(sf.crystal,params)
rlu_to_absolute_units!(swt::SpinWaveTheory,params::BinningParameters) = rlu_to_absolute_units!(swt.recipvecs_chem,params)


"""
    generate_shiver_script(params::BinningParameters)

Generate a Mantid/Shiver script which bins data according to the 
given [`BinningParameters`](@ref). You may want to call [`rlu_to_absolute_units!`](@ref) first.
"""
function generate_shiver_script(params)
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
    params, signal = load_nxs_binning_parameters(filename)

Given the name of a Mantid-exported `MDHistoWorkspace` file, load the [`BinningParameters`](@ref) and the signal from that file.
"""
function load_nxs_binning_parameters(filename)
    JLD2.jldopen(filename) do file
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

        # This is how you extract the covectors
        covectors = 2π .* transform_from_orig[1:3,1:3] * ub_matrix

        signal = file["MDHistoWorkspace"]["data"]["signal"]

        axes = Dict(JLD2.load_attributes(file,"MDHistoWorkspace/data/signal"))[:axes]

        # Axes are just stored backwards in Mantid .nxs for some reason
        axes_names = reverse(split(axes,":"))

        data_dims = Vector{Vector{Float64}}(undef,length(axes_names))
        binwidth = Vector{Float64}(undef,length(axes_names))
        binstart = Vector{Float64}(undef,length(axes_names))
        binend = Vector{Float64}(undef,length(axes_names))
        std = x -> sqrt(sum((x .- sum(x) ./ length(x)).^2))
        for (i,name) in enumerate(axes_names)
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

function quick_view_nxs(filename,keep_ax)
    integration_axes = setdiff(1:4,keep_ax)
    params, signal = load_nxs_binning_parameters(filename)
    integrate_axes!(params,axes = integration_axes)
    int_signal = dropdims(sum(signal,dims = integration_axes);dims = Tuple(integration_axes))
    bcs = axes_bincenters(params)
    (bcs[keep_ax[1]],bcs[keep_ax[2]],int_signal)
end

"""
    unit_resolution_binning_parameters(sf::StructureFactor)

Create [`BinningParameters`](@ref) which place one histogram bin centered at each possible `(q,ω)` scattering vector of the crystal.
This is the finest possible binning without creating bins with zero scattering vectors in them.

This function can be used without reference to a [`StructureFactor`](@ref) using this alternate syntax to manually specify the bin centers for the energy axis and the lattice size:

    unit_resolution_binning_parameters(ω_bincenters = ωs(sf),latsize = sf.latsize)
"""
function unit_resolution_binning_parameters(ωvals,latsize) 
    numbins = (latsize...,length(ωvals))
    # Bin centers should be at Sunny scattering vectors
    maxQ = 1 .- (1 ./ numbins)
    total_size = (maxQ[1],maxQ[2],maxQ[3],maximum(ωvals)) .- (0.,0.,0.,minimum(ωvals))
    binwidth = total_size ./ (numbins .- 1)
    binwidth = binwidth .+ eps.(binwidth)
    binstart = (0.,0.,0.,minimum(ωvals)) .- (binwidth ./ 2)
    binend = (maxQ[1],maxQ[2],maxQ[3],maximum(ωvals)) .+ (binwidth ./ 2)

    return BinningParameters(binstart,binend,binwidth)
end

unit_resolution_binning_parameters(sf::StructureFactor) = unit_resolution_binning_parameters(ωs(sf),sf.latsize)

function unit_resolution_binning_parameters(ωvals::Vector{Float64})
    ωbinwidth = (maximum(ωvals) - minimum(ωvals)) / (length(ωvals) - 1)
    ωbinwidth += eps(ωbinwidth)
    ωstart = minimum(ωvals) - ωbinwidth / 2
    ωend = maximum(ωvals) + ωbinwidth / 2
    return ωstart, ωend, ωbinwidth
end

"""
    one_dimensional_cut_binning_parameter(sf::StructureFactor, cut_from_q, cut_to_q, cut_bins::Int64, cut_width::Float64; plane_normal = [0,0,1],cut_height = cutwidth)

Creates [`BinningParameters`](@ref) which make a 1D cut in Q-space.
 
The x-axis of the resulting histogram consists of `cut_bins`-many bins ranging
from `cut_from_q` to `cut_to_q`. The orientation of the binning in the transverse directions is
determined automatically using `plane_normal`.
The width of the bins in the transverse direciton is controlled by `cut_width` and `cut_height`.

If the cut is too narrow, there will be very few scattering vectors per bin, or
the number per bin will vary substantially along the cut.
If the output appears under-resolved, try increasing `cut_width`.

The four axes of the resulting histogram are:
  1. Along the cut
  2. Fist transverse Q direction
  3. Second transverse Q direction
  4. Energy

This function can be used without reference to a [`StructureFactor`](@ref) using this alternate syntax to manually specify the bin centers for the energy axis:

    one_dimensional_cut_binning_parameter(ω_bincenters = ωs(sf),...)
"""
function one_dimensional_cut_binning_parameters(ωvals::Vector{Float64},cut_from_q,cut_to_q,cut_bins::Int64,cut_width;plane_normal = [0,0,1],cut_height = cut_width)
    # This covector should measure progress along the cut in r.l.u.
    cut_covector = normalize(cut_to_q - cut_from_q)
    # These two covectors should be perpendicular to the cut, and to each other
    transverse_covector = normalize(plane_normal × cut_covector)
    cotransverse_covector = normalize(transverse_covector × cut_covector)

    start_x = cut_covector ⋅ cut_from_q
    end_x = cut_covector ⋅ cut_to_q

    transverse_center = transverse_covector ⋅ cut_from_q # Equal to using cut_to_q
    cotransverse_center = cotransverse_covector ⋅ cut_from_q

    ωstart, ωend, ωbinwidth = unit_resolution_binning_parameters(ωvals)


    binstart = [start_x,transverse_center - cut_width/2,cotransverse_center - cut_height/2,ωstart]
    binend = [end_x,transverse_center + cut_width/2,cotransverse_center + cut_height/2,ωend]
    numbins = [cut_bins,1,1,length(ωvals)]
    covectors = [cut_covector... 0; transverse_covector... 0; cotransverse_covector... 0; 0 0 0 1]

    return BinningParameters(binstart,binend;numbins = numbins, covectors = covectors)
end
one_dimensional_cut_binning_parameters(sf::StructureFactor,cut_from_q,cut_to_q,cut_bins,cut_width;kwargs...) = one_dimensional_cut_binning_parameters(ωs(sf),cut_from_q,cut_to_q,cut_bins,cut_width;kwargs...)

"""
    axes_bincenters(params::BinningParameters)

Returns tick marks which label the bins of the histogram described by [`BinningParameters`](@ref) by their bin centers.

The following alternative syntax can be used to compute bin centers for a single axis:

    axes_bincenters(binstart,binend,binwidth)
"""
function axes_bincenters(binstart,binend,binwidth)
    bincenters = []
    for k = 1:length(binstart)
        first_center = binstart[k] .+ binwidth[k] ./ 2
        nbin = count_bins(binstart[k],binend[k],binwidth[k])
        push!(bincenters,range(first_center,step = binwidth[k],length = nbin))
    end
    bincenters
end
axes_bincenters(params::BinningParameters) = axes_bincenters(params.binstart,params.binend,params.binwidth)

"""
    connected_path_bins(sf,qs,density,args...;kwargs...)

Takes a list of wave vectors, `qs`, and builds a series of histogram `BinningParameters`
whose first axis traces a path through the provided points.
The second and third axes are integrated over according to the `args` and `kwargs`,
which are passed through to `one_dimensional_cut_binning_parameters`.

Also returned is a list of marker indices corresponding to the input points, and
a list of ranges giving the indices of each histogram `x`-axis within a concatenated histogram.
The `density` parameter is given in samples per reciprocal lattice unit (R.L.U.).
"""
function connected_path_bins(recip_vecs,ωvals,qs,density,args...;kwargs...)
    nPts = length(qs)
    params = []
    markers = []
    ranges = []
    total_bins_so_far = 0
    push!(markers, total_bins_so_far+1)
    for k = 1:(nPts-1)
        startPt = qs[k]
        endPt = qs[k+1]
        dist = norm(recip_vecs*(endPt - startPt))
        nBins = round(Int64,density * norm(endPt-startPt))
        push!(params,one_dimensional_cut_binning_parameters(ωvals,startPt,endPt,nBins,args...;kwargs...))
        push!(ranges, total_bins_so_far .+ (1:nBins))
        total_bins_so_far = total_bins_so_far + nBins
        push!(markers, total_bins_so_far+1)
    end
    return params, markers, ranges
end
connected_path_bins(sf::StructureFactor, qs::Vector, density,args...;kwargs...) = connected_path_bins(2π*inv(sf.crystal.latvecs)', ωs(sf), qs, density,args...;kwargs...)
connected_path_bins(sw::SpinWaveTheory, ωvals, qs::Vector, density,args...;kwargs...) = connected_path_bins(sw.recipvecs_chem, ωvals, qs, density,args...;kwargs...)


"""
    intensity, counts = intensities_binned(sf::StructureFactor, params::BinningParameters; formula, integrated_kernel)

Given correlation data contained in a [`StructureFactor`](@ref) and [`BinningParameters`](@ref) describing the
shape of a histogram, compute the intensity and normalization for each histogram bin using a given [`intensity_fomula`](@ref), or `intensity_formula(sf,:perp)` by default.

This is an alternative to [`intensities_interpolated`](@ref) which bins the scattering intensities into a histogram
instead of interpolating between them at specified `qs` values. See [`unit_resolution_binning_parameters`](@ref)
for a reasonable default choice of [`BinningParameters`](@ref) which roughly emulates [`intensities_interpolated`](@ref) with `interpolation = :round`.

If a function `integrated_kernel(Δω)` is passed, it will be used as the CDF of a kernel function for energy broadening.
For example,
`integrated_kernel = Δω -> atan(Δω/η)/pi` (c.f. [`integrated_lorentzian`](@ref) implements Lorentzian broadening with parameter `η`.
Currently, energy broadening is only supported if the [`BinningParameters`](@ref) are such that the first three axes are purely spatial and the last (energy) axis is `[0,0,0,1]`.
"""
function intensities_binned(sf::StructureFactor, params::BinningParameters;
    integrated_kernel = nothing,formula = intensity_formula(sf,:perp)
)
    (; binwidth, binstart, binend, covectors, numbins) = params
    output_intensities = zeros(Float64,numbins...)
    output_counts = zeros(Float64,numbins...)
    ωvals = ωs(sf)
    recip_vecs = 2π*inv(sf.crystal.latvecs)'

    # Find an axis-aligned bounding box containing the histogram
    lower_aabb_q, upper_aabb_q = binning_parameters_aabb(params)

    # Round the axis-aligned bounding box *outwards* to lattice sites
    # SQTODO: are these bounds optimal?
    Ls = sf.latsize
    lower_aabb_cell = floor.(Int64,lower_aabb_q .* Ls .+ 1) 
    upper_aabb_cell = ceil.(Int64,upper_aabb_q .* Ls .+ 1)

    # Loop over every scattering vector in the bounding box
    for cell in CartesianIndices(Tuple(((:).(lower_aabb_cell,upper_aabb_cell))))
        base_cell = CartesianIndex(mod1.(cell.I,Ls)...)
        for (iω,ω) in enumerate(ωvals)

            # Compute intensity
            # [c.f. all_exact_wave_vectors, but we need `cell' index as well here]
            q = SVector((cell.I .- 1) ./ Ls) # q is in R.L.U.

            # Figure out which bin this goes in
            v = [q...,ω]
            coords = covectors * v
            xyztBin = 1 .+ floor.(Int64,(coords .- binstart) ./ binwidth)

            if isnothing(integrated_kernel) # `Delta-function energy' logic
                # Check this bin is within the 4D histogram bounds
                if all(xyztBin .<= numbins) && all(xyztBin .>= 1)
                    ci = CartesianIndex(xyztBin.data)
                    k = recip_vecs * q
                    intensity = formula.calc_intensity(sf,k,base_cell,iω)
                    output_intensities[ci] += intensity
                    output_counts[ci] += 1
                end
            else # `Energy broadening into bins' logic
                # For now, only support broadening for `simple' energy axes
                if covectors[4,:] == [0,0,0,1] && norm(covectors[1:3,:] * [0,0,0,1]) == 0

                    # Check this bin is within the *spatial* 3D histogram bounds
                    # If we are energy-broadening, then scattering vectors outside the histogram
                    # in the energy direction need to be considered
                    if all(xyztBin[1:3] .<= numbins[1:3]) &&  all(xyztBin[1:3] .>= 1)

                        # Calculate source scattering vector intensity only once
                        ci = CartesianIndex(xyztBin.data)
                        k = recip_vecs * q
                        intensity = formula.calc_intensity(sf,k,base_cell,iω)
                        # Broaden from the source scattering vector (k,ω) to
                        # each target bin ci_other
                        for iωother = 1:numbins[4]
                            ci_other = CartesianIndex(ci[1],ci[2],ci[3],iωother)
                            # Start and end points of the target bin
                            a = binstart[4] + (iωother - 1) * binwidth[4]
                            b = binstart[4] + iωother * binwidth[4]

                            # P(ω picked up in bin [a,b]) = ∫ₐᵇ Kernel(ω' - ω) dω'
                            fraction_in_bin = integrated_kernel(b - ω) - integrated_kernel(a - ω)
                            output_intensities[ci_other] += fraction_in_bin * intensity
                            output_counts[ci_other] += fraction_in_bin
                        end
                    end
                else
                    error("Energy broadening not yet implemented for histograms with complicated energy axes")
                end
            end
        end
    end
    return output_intensities, output_counts
end

"""
    intensities_bin_centers(swt::SpinWaveTheory, params::BinningParameters, η::Float64)

Computes the unpolarized inelastic neutron scattering intensities given a
`SpinWaveTheory`, histogram described by its `BinningParameters`, and
a Lorentzian broadening parameter `η`.

Note that this method only calculates the intensity at the bin centers--it doesn't
integrate over the bins in any way. The output will be the same shape as if it were
histogrammed data.
"""
function intensities_bin_centers(swt::SpinWaveTheory, params::BinningParameters, η::Float64)
    (; sys) = swt
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    nmodes  = Nf * Nm

    bin_centers = axes_bincenters(params)
    qs = []
    ωs = []
    # Record all histogram bin centers
    for ci in CartesianIndices(params.numbins.data)
        qx_center = bin_centers[1][ci[1]]
        qy_center = bin_centers[2][ci[2]]
        qz_center = bin_centers[3][ci[3]]

        q_center = [qx_center,qy_center,qz_center]
        push!(qs,q_center)
        ω_center = bin_centers[4][ci[4]]
        push!(ωs,ω_center)
    end

    # Compute SWT at bin center qs
    disp, Sαβs = dssf(swt, qs)

    is = zeros(Float64,params.numbins...)
    for (cii,ci) in enumerate(CartesianIndices(params.numbins.data))
        q = qs[cii]
        polar_mat = polarization_matrix(swt.recipvecs_chem * q)

        for band = 1:nmodes
            band_intensity = real(sum(polar_mat .* Sαβs[cii,band]))
            # At a Goldstone mode, where the intensity is divergent, use a delta-function for the intensity.
            if (disp[cii, band] < 1.0e-3) && (band_intensity > 1.0e3)
                is[ci] += band_intensity
            else
                #SQTODO: This calculation is fake. It needs to integrate over the bin.
                is[ci] += band_intensity * lorentzian(ωs[cii]-disp[cii,band], η)
            end
        end
    end
    return is
end


