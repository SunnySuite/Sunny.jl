
"""
    BinningParameters(binstart,binend,binwidth;covectors = I(4))
    BinningParameters(binstart,binend;numbins,covectors = I(4))

Describes a 4D parallelepided histogram in a format compatible with experimental Inelasitic Neutron Scattering data.
See [`generate_mantid_script_from_binning_parameters`](@ref) to convert [`BinningParameters`](@ref) to a format understandable by the [Mantid software](https://www.mantidproject.org/), or [`load_nxs`](@ref) to load [`BinningParameters`](@ref) from a Mantid `.nxs` file.
 
The coordinates of the histogram axes are specified by multiplication 
of `(k,ω)` with each row of the `covectors` matrix, with `k` given in absolute units.
Since the default `covectors` matrix is the identity matrix, the default axes are
`(kx,ky,kz,ω)` in absolute units.
To bin `(q,ω)` values given in reciprocal lattice units (R.L.U.) instead, see [`bin_rlu_as_absolute_units!`](@ref).

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
        bin_edges = axes_binedges(params)
        first_edges = map(x -> x[1],bin_edges)
        last_edges = map(x -> x[end],bin_edges)
        @printf(io," from %+.3f to %+.3f along [", first_edges[k], last_edges[k])
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

Base.copy(p::BinningParameters) = BinningParameters(copy(p.binstart),copy(p.binend),copy(p.binwidth),copy(p.covectors))

# Support numbins as a (virtual) property, even though only the binwidth is stored
Base.getproperty(params::BinningParameters, sym::Symbol) = sym == :numbins ? count_bins(params.binstart,params.binend,params.binwidth) : getfield(params,sym)

function Base.setproperty!(params::BinningParameters, sym::Symbol, numbins)
    if sym == :numbins
        # *Ensure* that the last bin contains params.binend
        params.binwidth .= (params.binend .- params.binstart) ./ (numbins .- 0.5)
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
    params = BinningParameters(binstart,binend,[0.,0,0,0];kwargs...)
    params.numbins = numbins # Use the setproperty to do it correctly
    params
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
    bin_edges = axes_binedges(params)
    first_edges = map(x -> x[1],bin_edges)
    last_edges = map(x -> x[end],bin_edges)
    bin_edges = [first_edges last_edges]
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
If `params` expects to bin values `(k,ω)` in absolute units, then calling

    bin_rlu_as_absolute_units!(params::BinningParameters,[reciprocal lattice vectors])

will modifiy the `covectors` in `params` so that they will accept `(q,ω)` in Reciprocal Lattice Units (R.L.U.) instead.
Conversly, if `params` expects `(q,ω)` R.L.U., calling

    bin_absolute_units_as_rlu!(params::BinningParameters,[reciprocal lattice vectors])

will adjust `params` to instead accept `(k,ω)` absolute units.

The second argument may be a 3x3 matrix specifying the reciprocal lattice vectors, or any of these objects:
- [`Crystal`](@ref)
- [`System`](@ref)
- [`SampledCorrelations`](@ref)
- [`SpinWaveTheory`](@ref)
"""
bin_absolute_units_as_rlu!, bin_rlu_as_absolute_units!

function bin_rlu_as_absolute_units!(params::BinningParameters,recip_vecs::AbstractMatrix)
    covectorsK = params.covectors

    # covectorsQ * q = covectorsK *  recip_vecs * q = covectorsK * k
    # covectorsQ     = covectorsK *  recip_vecs
    covectorsQ       = covectorsK * [recip_vecs [0;0;0]; [0 0 0] 1]
    params.covectors = MMatrix{4,4}(covectorsQ)
    params
end

# covectorsK * k = covectorsQ * inv(recip_vecs) * k = covectorsQ * q
# covectorsK     = covectorsQ * inv(recip_vecs)
bin_absolute_units_as_rlu!(params::BinningParameters,recip_vecs::AbstractMatrix) = bin_rlu_as_absolute_units!(params,inv(recip_vecs))

bin_absolute_units_as_rlu!(params::BinningParameters,crystal::Crystal) = bin_absolute_units_as_rlu!(params,2π*inv(crystal.latvecs)')
bin_absolute_units_as_rlu!(params::BinningParameters,sys::System) = bin_absolute_units_as_rlu!(params,sys.crystal)
bin_absolute_units_as_rlu!(params::BinningParameters,sc::SampledCorrelations) = bin_absolute_units_as_rlu!(params,sc.crystal)
bin_absolute_units_as_rlu!(params::BinningParameters,swt::SpinWaveTheory) = bin_absolute_units_as_rlu!(params,swt.sys)

bin_rlu_as_absolute_units!(params::BinningParameters,crystal::Crystal) = bin_absolute_units_as_rlu!(params,crystal.latvecs'/2π)
bin_rlu_as_absolute_units!(params::BinningParameters,sys::System) = bin_rlu_as_absolute_units!(params,sys.crystal)
bin_rlu_as_absolute_units!(params::BinningParameters,sc::SampledCorrelations) = bin_rlu_as_absolute_units!(params,sc.crystal)
bin_rlu_as_absolute_units!(params::BinningParameters,swt::SpinWaveTheory) = bin_rlu_as_absolute_units!(params,swt.sys)


"""
    unit_resolution_binning_parameters(sc::SampledCorrelations; units = :absolute)

Create [`BinningParameters`](@ref) which place one histogram bin centered at each possible `(q,ω)` scattering vector of the crystal.
This is the finest possible binning without creating bins with zero scattering vectors in them.

Setting `units = :rlu` returns [`BinningParameters`](@ref) which accept values in R.L.U. instead of absolute units.

This function can be used without reference to a [`SampledCorrelations`](@ref) using an alternate syntax to manually specify the bin centers for the energy axis and the lattice size:

    unit_resolution_binning_parameters(ω_bincenters,latsize,[reciprocal lattice vectors]; [units])

As in [`bin_absolute_units_as_rlu!`](@ref), the last argument may be a 3x3 matrix specifying the reciprocal lattice vectors, or any of these objects:
- [`Crystal`](@ref)
- [`System`](@ref)
- [`SampledCorrelations`](@ref)
- [`SpinWaveTheory`](@ref)

Lastly, binning parameters for a single axis may be specifed by their bin centers:

    (binstart,binend,binwidth) = unit_resolution_binning_parameters(bincenters::Vector{Float64})
"""
function unit_resolution_binning_parameters(ωvals,latsize,args...;units = :absolute)
    numbins = (latsize...,length(ωvals))
    # Bin centers should be at Sunny scattering vectors
    maxQ = 1 .- (1 ./ numbins)
    
    min_val = (0.,0.,0.,minimum(ωvals))
    max_val = (maxQ[1],maxQ[2],maxQ[3],maximum(ωvals))
    total_size = max_val .- min_val

    binwidth = total_size ./ (numbins .- 1)
    binstart = (0.,0.,0.,minimum(ωvals)) .- (binwidth ./ 2)
    binend = (maxQ[1],maxQ[2],maxQ[3],maximum(ωvals)) # bin end is well inside of last bin

    params = BinningParameters(binstart,binend,binwidth)

    # Special case for when there is only one bin in a direction
    for i = 1:4
        if numbins[i] == 1
            params.binwidth[i] = 1.
            params.binstart[i] = min_val[i] - (params.binwidth[i] ./ 2)
            params.binend[i] = min_val[i]
        end
    end

    if units == :absolute
        bin_absolute_units_as_rlu!(params,args...)
    elseif units == :rlu
        params
    else
        error("unit_resolution_binning_parameters: Unsupported units $units")
    end
end

unit_resolution_binning_parameters(sc::SampledCorrelations; kwargs...) = unit_resolution_binning_parameters(ωs(sc),sc.latsize,sc;kwargs...)

function unit_resolution_binning_parameters(ωvals::AbstractVector{Float64})
    if !all(abs.(diff(diff(ωvals))) .< 1e-12)
      @warn "Non-uniform bins will be re-spaced into uniform bins"
    end
    if length(ωvals) == 1
      error("Can not infer bin width given only one bin center")
    end
    ωbinwidth = (maximum(ωvals) - minimum(ωvals)) / (length(ωvals) - 1)
    ωstart = minimum(ωvals) - ωbinwidth / 2
    ωend = maximum(ωvals)

    return ωstart, ωend, ωbinwidth
end

"""
    slice_2D_binning_parameter(sc::SampledCorrelations, cut_from_q, cut_to_q, cut_bins::Int64, cut_width::Float64; plane_normal = [0,0,1],cut_height = cutwidth, units = :absolute)

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

Setting `units = :rlu` returns [`BinningParameters`](@ref) which accept values in R.L.U. instead of absolute units.

This function can be used without reference to a [`SampledCorrelations`](@ref) using this alternate syntax to manually specify the bin centers for the energy axis:

    slice_2D_binning_parameter(ω_bincenters, cut_from, cut_to,...)

where `ω_bincenters` specifies the energy axis, and both `cut_from` and `cut_to` are arbitrary covectors, in any units.
"""
function slice_2D_binning_parameters(ωvals::Vector{Float64},cut_from_q,cut_to_q,cut_bins::Int64,cut_width;plane_normal = [0,0,1],cut_height = cut_width)
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
    xstart, xend, xbinwidth = unit_resolution_binning_parameters(range(start_x,end_x,length = cut_bins))

    binstart = [xstart,transverse_center - cut_width/2,cotransverse_center - cut_height/2,ωstart]
    binend = [xend,transverse_center,cotransverse_center,ωend]
    numbins = [cut_bins,1,1,length(ωvals)]
    covectors = [cut_covector... 0; transverse_covector... 0; cotransverse_covector... 0; 0 0 0 1]

    BinningParameters(binstart,binend;numbins = numbins, covectors = covectors)
end

function slice_2D_binning_parameters(sc::SampledCorrelations,cut_from_q,cut_to_q,args...;units = :absolute,kwargs...)
    params = slice_2D_binning_parameters(ωs(sc),cut_from_q,cut_to_q,args...;kwargs...)

    if units == :absolute
        bin_absolute_units_as_rlu!(params,sc)
    elseif units == :rlu
        params
    else
        error("slice_2D_binning_parameters: Unsupported units $units")
    end
end

"""
    axes_bincenters(params::BinningParameters)

Returns tick marks which label the bins of the histogram described by [`BinningParameters`](@ref) by their bin centers.

The following alternative syntax can be used to compute bin centers for a single axis:

    axes_bincenters(binstart,binend,binwidth)
"""
function axes_bincenters(binstart,binend,binwidth)
    bincenters = Vector{AbstractRange{Float64}}(undef,0)
    for k = 1:length(binstart)
        first_center = binstart[k] .+ binwidth[k] ./ 2
        nbin = count_bins(binstart[k],binend[k],binwidth[k])
        push!(bincenters,range(first_center,step = binwidth[k],length = nbin))
    end
    bincenters
end
axes_bincenters(params::BinningParameters) = axes_bincenters(params.binstart,params.binend,params.binwidth)

function axes_binedges(binstart,binend,binwidth)
    binedges = Vector{AbstractRange{Float64}}(undef,0)
    for k = 1:length(binstart)
        nbin = count_bins(binstart[k],binend[k],binwidth[k])
        push!(binedges,range(binstart[k],step = binwidth[k],length = nbin + 1))
    end
    binedges
end
axes_binedges(params::BinningParameters) = axes_binedges(params.binstart,params.binend,params.binwidth)



"""
    connected_path_bins(sc,qs,density,args...;kwargs...)

Takes a list of wave vectors, `qs`, and builds a series of histogram [`BinningParameters`](@ref)
whose first axis traces a path through the provided points.
The second and third axes are integrated over according to the `args` and `kwargs`,
which are passed through to [`slice_2D_binning_parameters`](@ref).

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
        # Density is taken in R.L.U. since that's where the
        # scattering vectors are equally spaced!
        nBins = round(Int64,density * norm(endPt - startPt))
        # SQTODO: Automatic density that adjusts itself lower
        # if there are not enough (e.g. zero) counts in some bins

        param = slice_2D_binning_parameters(ωvals,startPt,endPt,nBins,args...;kwargs...)
        bin_absolute_units_as_rlu!(param,recip_vecs)
        push!(params,param)
        push!(ranges, total_bins_so_far .+ (1:nBins))
        total_bins_so_far = total_bins_so_far + nBins
        push!(markers, total_bins_so_far+1)
    end
    return params, markers, ranges
end
connected_path_bins(sc::SampledCorrelations, qs::Vector, density,args...;kwargs...) = connected_path_bins(2π*inv(sc.crystal.latvecs)', ωs(sc), qs, density,args...;kwargs...)
connected_path_bins(sw::SpinWaveTheory, ωvals, qs::Vector, density,args...;kwargs...) = connected_path_bins(sw.recipvecs_chem, ωvals, qs, density,args...;kwargs...)


"""
    intensity, counts = intensities_binned(sc::SampledCorrelations, params::BinningParameters; formula, integrated_kernel)

Given correlation data contained in a [`SampledCorrelations`](@ref) and [`BinningParameters`](@ref) describing the
shape of a histogram, compute the intensity and normalization for each histogram bin using a given [`intensity_formula`](@ref), or `intensity_formula(sc,:perp)` by default.

The [`BinningParameters`](@ref) are expected to accept `(k,ω)` in absolute units.

This is an alternative to [`intensities_interpolated`](@ref) which bins the scattering intensities into a histogram
instead of interpolating between them at specified `qs` values. See [`unit_resolution_binning_parameters`](@ref)
for a reasonable default choice of [`BinningParameters`](@ref) which roughly emulates [`intensities_interpolated`](@ref) with `interpolation = :round`.

If a function `integrated_kernel(Δω)` is passed, it will be used as the CDF of a kernel function for energy broadening.
For example,
`integrated_kernel = Δω -> atan(Δω/η)/pi` (c.f. [`integrated_lorentzian`](@ref) implements Lorentzian broadening with parameter `η`.
Currently, energy broadening is only supported if the [`BinningParameters`](@ref) are such that the first three axes are purely spatial and the last (energy) axis is `[0,0,0,1]`.
"""
function intensities_binned(sc::SampledCorrelations, params::BinningParameters;
    integrated_kernel = nothing,formula = intensity_formula(sc,:perp) :: ClassicalIntensityFormula
)
    (; binwidth, binstart, binend, covectors, numbins) = params
    output_intensities = zeros(Float64,numbins...)
    output_counts = zeros(Float64,numbins...)
    ωvals = ωs(sc)
    recip_vecs = 2π*inv(sc.crystal.latvecs)'

    # Find an axis-aligned bounding box containing the histogram.
    # The AABB needs to be in q-space because that's where we index
    # over the scattering modes.
    q_params = copy(params)
    bin_rlu_as_absolute_units!(q_params,sc)
    lower_aabb_q, upper_aabb_q = binning_parameters_aabb(q_params)

    # Round the axis-aligned bounding box *outwards* to lattice sites
    # SQTODO: are these bounds optimal?
    Ls = sc.latsize
    lower_aabb_cell = floor.(Int64,lower_aabb_q .* Ls .+ 1) 
    upper_aabb_cell = ceil.(Int64,upper_aabb_q .* Ls .+ 1)

    q = MVector{3,Float64}(undef)
    v = MVector{4,Float64}(undef)
    k = view(v,1:3)
    coords = MVector{4,Float64}(undef)
    xyztBin = MVector{4,Int64}(undef)
    xyzBin = view(xyztBin,1:3)

    # Pre-compute discrete broadening kernel from continuous one provided
    if !isnothing(integrated_kernel)
        fraction_in_bin = Vector{Vector{Float64}}(undef,length(ωvals))
        for (iω,ω) in enumerate(ωvals)
            fraction_in_bin[iω] = Vector{Float64}(undef,numbins[4])
            for iωother = 1:numbins[4]
                ci_other = CartesianIndex(xyzBin[1],xyzBin[2],xyzBin[3],iωother)
                # Start and end points of the target bin
                a = binstart[4] + (iωother - 1) * binwidth[4]
                b = binstart[4] + iωother * binwidth[4]

                # P(ω picked up in bin [a,b]) = ∫ₐᵇ Kernel(ω' - ω) dω'
                fraction_in_bin[iω][iωother] = integrated_kernel(b - ω) - integrated_kernel(a - ω)
            end
        end
    end

    # Loop over every scattering vector in the bounding box
    for cell in CartesianIndices(Tuple(((:).(lower_aabb_cell,upper_aabb_cell))))
        # Which is the analog of this scattering mode in the first BZ?
        base_cell = CartesianIndex(mod1.(cell.I,Ls)...)
        q .= ((cell.I .- 1) ./ Ls) # q is in R.L.U.
        k .= recip_vecs * q # But binning is done in absolute units
        for (iω,ω) in enumerate(ωvals)
            if isnothing(integrated_kernel) # `Delta-function energy' logic
                # Figure out which bin this goes in
                v[4] = ω
                mul!(coords,covectors,v)
                xyztBin .= 1 .+ floor.(Int64,(coords .- binstart) ./ binwidth)

                # Check this bin is within the 4D histogram bounds
                if all(xyztBin .<= numbins) && all(xyztBin .>= 1)
                    intensity = formula.calc_intensity(sc,SVector{3,Float64}(k),base_cell,iω)

                    ci = CartesianIndex(xyztBin.data)
                    output_intensities[ci] += intensity
                    output_counts[ci] += 1
                end
            else # `Energy broadening into bins' logic
                # For now, only support broadening for `simple' energy axes
                if covectors[4,:] == [0,0,0,1] && norm(covectors[1:3,:] * [0,0,0,1]) == 0

                    # Check this bin is within the *spatial* 3D histogram bounds
                    # If we are energy-broadening, then scattering vectors outside the histogram
                    # in the energy direction need to be considered
                    mul!(view(coords,1:3),view(covectors,1:3,1:3), view(v,1:3))
                    xyzBin .= 1 .+ floor.(Int64,(view(coords,1:3) .- view(binstart,1:3)) ./ view(binwidth,1:3))
                    if all(xyzBin .<= view(numbins,1:3)) &&  all(xyzBin .>= 1)

                        # Calculate source scattering vector intensity only once
                        intensity = formula.calc_intensity(sc,SVector{3,Float64}(k),base_cell,iω)

                        # Broaden from the source scattering vector (k,ω) to
                        # each target bin ci_other
                        ci_other = CartesianIndex(xyzBin[1],xyzBin[2],xyzBin[3])
                        view(output_intensities,ci_other,:) .+= fraction_in_bin[iω] .* intensity
                        view(output_counts,ci_other,:) .+= fraction_in_bin[iω]
                    end
                else
                    error("Energy broadening not yet implemented for histograms with complicated energy axes")
                end
            end
        end
    end
    return output_intensities, output_counts
end

