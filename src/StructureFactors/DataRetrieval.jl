################################################################################
# Basic functions for retrieving 𝒮(𝐪,ω) values
################################################################################

# Internal function for getting a single 𝒮(q, ω) intensity
function calc_intensity(sf::StructureFactor, k, cell, ω, iω, contractor, kT, ffdata, ::Val{NCorr}, ::Val{NAtoms}) where {NCorr, NAtoms}
    elems = phase_averaged_elements(view(sf.data,:,:,:,cell,iω), k, sf, ffdata, Val(NCorr), Val(NAtoms))
    intensity = contract(elems, k, contractor)
    return intensity * classical_to_quantum(ω, kT)  # DDTodo: Probably make this a post-processing step 
end

classical_to_quantum(ω, kT::Float64) = ω > 0 ? ω/(kT*(1 - exp(-ω/kT))) : iszero(ω) ? 1.0 : -ω*exp(ω/kT)/(kT*(1 - exp(ω/kT)))
classical_to_quantum(ω, ::Nothing) = 1.0


function prepare_form_factors(sf, formfactors)
    if isnothing(formfactors)
        cryst = isnothing(sf.origin_crystal) ? sf.crystal : sf.origin_crystal 
        class_indices = [findfirst(==(class_label), cryst.classes) for class_label in unique(cryst.classes)]
        formfactors = [FormFactor{Sunny.EMPTY_FF}(; atom) for atom in class_indices]
    end
    formfactors = upconvert_form_factors(formfactors) # Ensure formfactors have consistent type
    return propagate_form_factors(sf, formfactors)
end


# Describes a 4D parallelepided histogram in a format compatible with Mantid/Horace
# 
# The coordinates of the histogram axes are specified by multiplication 
# of `(q,ω)' with each row of the covectors matrix
#
# The convention is that:
# - The left edge of the first bin starts at `binstart`
# - The last bin contains `binend`
# - There are no ``partial bins;'' the last bin may contain values greater than `binend`
# - The bin width is `binwidth`
#
# A value can be binned by computing its bin number:
# 
#     bin = 1 + floor(Int64,(value - binstart) / binwidth)
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

# This function defines how partial bins are handled.
count_bins(bin_start,bin_end,bin_width) = ceil.(Int64,(bin_end .- bin_start) ./ bin_width)

# Default coordinates are (Qx,Qy,Qz,ω)
function BinningParameters(binstart,binend,binwidth;covectors = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    return BinningParameters(binstart,binend,binwidth,covectors)
end

function BinningParameters(binstart,binend;numbins,kwargs...)
    binwidth = (binend .- binstart) ./ numbins
    binwidth .+= eps.(binwidth)
    return BinningParameters(binstart,binend,binwidth;kwargs...)
end

# Integrate over one or more axes of the histogram by setting the number of bins
# in that axis to 1
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

function rlu_to_absolute_units!(sf,params::BinningParameters)
    recip_vecs = 2π*inv(sf.crystal.latvecs)'
    covectorsQ = params.covectors

    # covectorsQ * q = covectorsK * k = covectorsK * recip_vecs * q
    covectorsK = covectorsQ * [inv(recip_vecs) [0;0;0]; [0 0 0] 1]
    params.covectors = MMatrix{4,4}(covectorsK)
end


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

# This places one histogram bin around each possible Sunny scattering vector.
# This is the finest possible binning without creating bins with zero scattering vectors in them.
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

# Creates BinningParameters which implement a 1D cut in (Qx,Qy,Qz) space.
# 
# The x-axis of the resulting histogram consists of `cut_bins`-many bins ranging
# from `cut_from_q` to `cut_to_q`. The binning in the transverse directions is
# determined automatically using `plane_normal`, and has size controlled by `cut_width`.
#
# If the cut is too narrow, there will be very few scattering vectors per bin, or
# the number per bin will vary substantially along the cut.
#
# The four axes of the resulting histogram are:
#   1. Along the cut
#   2. Fist transverse Q direction
#   3. Second transverse Q direction
#   4. Energy
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

# Returns tick marks which label the histogram bins by their bin centers
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



# Given correlation data contained in `sf` and BinningParameters describing the
# shape of a histogram, compute the intensity and normalization for each histogram bin.
#
# This is an alternative to `intensities` which bins the scattering intensities into a histogram
# instead of interpolating between them at specified `qs` values. See `unit_resolution_binning_parameters`
# for a reasonable default choice of `BinningParameters` which roughly emulates `intensities` with `NoInterp`.
#
# If a function `integrated_kernel(Δω)` is passed, it will be used as the CDF of a kernel function for energy broadening.
# For example,
# `integrated_kernel = Δω -> atan(Δω/η)/pi` implements `lorentzian` broadening with parameter `η`.
# Currently, energy broadening is only supported if the `BinningParameters` are such that the first three axes are purely spatial and the last (energy) axis is `[0,0,0,1]`.
function intensities_binned(sf::StructureFactor, params::BinningParameters, mode;
    integrated_kernel = nothing,
    kT = nothing,
    formfactors = nothing,
)
    (; binwidth, binstart, binend, covectors, numbins) = params
    output_intensities = zeros(Float64,numbins...)
    output_counts = zeros(Float64,numbins...)
    ωvals = ωs(sf)
    recip_vecs = 2π*inv(sf.crystal.latvecs)'
    ffdata = prepare_form_factors(sf, formfactors)
    contractor = if mode == :perp
        DipoleFactor(sf)
    elseif typeof(mode) <: Tuple{Int, Int}
        Element(sf, mode)
    else
        Trace(sf)
    end

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
                    NCorr, NAtoms = size(sf.data)[1:2]
                    intensity = calc_intensity(sf,k,base_cell,ω,iω,contractor, kT, ffdata, Val(NCorr), Val(NAtoms))
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
                        NCorr, NAtoms = size(sf.data)[1:2]
                        intensity = calc_intensity(sf,k,base_cell,ω,iω,contractor, kT, ffdata, Val(NCorr), Val(NAtoms))
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


# Note that requests for intensities often come in lists of nearby q values.
# Since the data is inherently discretized, this often results in repeated calls
# for values at the same discrete points. Since basis reduction is done for each
# of this calls, this results in a large amount of repeated calculation. This
# function analyzes repetitions in advance and prunes them out. This is
# ugly, but the speedup when tested on a few simple, realistic examples was
# 3-5x.
function pruned_stencil_info(sf::StructureFactor, qs, interp::InterpolationScheme{N}) where N
    # Count the number of contiguous regions with unchanging values. If all
    # values are unique, returns the length of q_info. Note comparison is on m
    # values rather than index values and the m values are the first element of
    # the a tuple, that is, we're checking x[1] == y[1] in the map.
    m_info = map(q -> stencil_points(sf, q, interp), qs)
    numregions = sum(map((x,y) -> x[1] == y[1] ? 0 : 1, m_info[1:end-1], m_info[2:end])) + 1
    
    # Remove repeated stencil points and count number of instances of each
    ms_ref, idcs_ref = stencil_points(sf, qs[1], interp)
    ms_all  = fill(ntuple(x->zero(Vec3), N), numregions)
    ms_all[1] = ms_ref 
    idcs_all = fill(ntuple(x->CartesianIndex((-1,-1,-1)), N), numregions)
    idcs_all[1] = idcs_ref 
    counts = zeros(Int64, numregions)
    c = counts[1] = 1
    for q in qs[2:end] 
        ms, idcs = stencil_points(sf, q, interp)
        if ms != ms_ref
            ms_ref = ms 
            c += 1
            ms_all[c] =  ms
            idcs_all[c] = idcs 
        end
        counts[c] += 1
    end
    @assert sum(counts) == length(m_info)

    # Calculate corresponding q (RLU) and k (global) vectors
    recip_vecs = 2π*inv(sf.crystal.latvecs)'  # Note, qs will be in terms of sf.crystal by this point, not origin_crystal
    qs_all = map(ms_all) do ms
       map(m -> m ./ sf.latsize, ms) 
    end

    ks_all = map(qs_all) do qs
        map(q -> recip_vecs * q, qs)
    end
    
    return (; qs_all, ks_all, idcs_all, counts)
end


"""
    intensities(sf::StructureFactor, qs, mode; interpolation = nothing,
                    kT = nothing, formfactors = nothing, negative_energies = false)

The basic function for retrieving ``𝒮(𝐪,ω)`` information from a
`StructureFactor`. Maps an array of wave vectors `qs` to an array of structure
factor intensities, including an additional energy index. The values of ``ω``
associated with the energy index can be retrieved by calling [`ωs`](@ref). The
three coordinates of each wave vector are measured in reciprocal lattice units,
i.e., multiples of the reciprocal lattice vectors.

- `mode`: Should be one of `:trace`, `:perp`, or `:full`. Determines an optional
    contraction on the indices ``α`` and ``β`` of ``𝒮^{αβ}(q,ω)``. Setting
    `trace` yields ``∑_α 𝒮^{αα}(q,ω)``. Setting `perp` will contract
    ``𝒮^{αβ}(q,ω)`` with the dipole factor ``δ_{αβ} - q_{α}q_{β}``, returning
    the unpolarized intensity. Setting `full` will return all elements
    ``𝒮^{αβ}(𝐪,ω)`` without contraction.
- `interpolation`: Since ``𝒮(𝐪, ω)`` is calculated on a finite lattice, data
    is only available at discrete wave vectors. By default, Sunny will round a
    requested `q` to the nearest available wave vector. Linear interpolation can
    be applied by setting `interpolation=:linear`.
- `kT`: If a temperature is provided, the intensities will be rescaled by a
    temperature- and ω-dependent classical-to-quantum factor. `kT` should be
    specified when making comparisons with spin wave calculations or
    experimental data.
- `formfactors`: To apply form factor corrections, provide this keyword with a
    vector of `FormFactor`s, one for each unique site in the unit cell. Sunny
    will symmetry propagate the results to all equivalent sites.
- `negative_energies`: If set to `true`, Sunny will return the periodic
    extension of the energy axis. Most users will not want this.
"""
function intensities(sf::StructureFactor, qs, mode;
    interpolation = :none,
    kT = nothing,
    formfactors = nothing,
    negative_energies = false,
    static_warn = true
)
    qs = Vec3.(qs)
    NCorr, NAtoms = size(sf.data)[1:2]

    # If working on reshaped system, assume qs given as coordinates in terms of
    # reciprocal vectors of original crystal and convert them to qs in terms of
    # the reciprocal vectors of the reshaped crystal.
    if !isnothing(sf.origin_crystal)
        rvecs_reshaped = inv(sf.crystal.latvecs)'       # Note, leading 2π will cancel
        rvecs_origin = inv(sf.origin_crystal.latvecs)'
        qs = map(q -> rvecs_reshaped \ rvecs_origin * q, qs)
    end

    # Make sure it's a dynamical structure factor 
    if static_warn && size(sf.data, 7) == 1
        error("`intensities` given a StructureFactor with no dynamical information. Call `instant_intensities` to retrieve instantaneous (static) structure factor data.")
    end

    # If temperature given, ensure it's greater than 0.0
    if !isnothing(kT) && iszero(kT)
        error("`kT` must be greater than zero.")
    end

    # Set up interpolation scheme
    interp = if interpolation == :none
        NoInterp()
    elseif interpolation == :linear
        LinearInterp()
    end

    # Set up element contraction
    contractor = if mode == :trace
        Trace(sf)
    elseif mode == :perp
        DipoleFactor(sf)
    elseif mode == :full
        FullTensor(sf)
    elseif typeof(mode) <: Tuple{Int, Int}
        Element(sf, mode)
    end

    # Propagate form factor information (if any)
    ffdata = prepare_form_factors(sf, formfactors) 

    # Precompute index information and preallocate
    ωvals = ωs(sf; negative_energies)
    nω = length(ωvals) 
    stencil_info = pruned_stencil_info(sf, qs, interp) 
    intensities = zeros(contractor, size(qs)..., nω)
    
    # Call type stable version of the function
    intensities!(intensities, sf, qs, ωvals, interp, contractor, kT, ffdata, stencil_info, Val(NCorr), Val(NAtoms)) 

    return intensities
end


# Actual intensity calculation
function intensities!(intensities, sf::StructureFactor, q_targets::Array, ωvals, interp::InterpolationScheme, contraction::Contraction{T}, temp, ffdata, stencil_info, ::Val{NCorr}, ::Val{NAtoms}) where {T, NCorr, NAtoms}
    li_intensities = LinearIndices(intensities)
    ci_qs = CartesianIndices(q_targets)
    (; qs_all, ks_all, idcs_all, counts) = stencil_info 
    for (iω, ω) in enumerate(ωvals)
        iq = 0
        for (qs, ks, idcs, numrepeats) in zip(qs_all, ks_all, idcs_all, counts)
            local_intensities = stencil_intensities(sf, ks, idcs, ω, iω, interp, contraction, temp, ffdata, Val(NCorr), Val(NAtoms)) 
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(ci_qs[iq], iω)]
                intensities[idx] = interpolated_intensity(sf, q_targets[iq], qs, local_intensities, interp) 
            end
        end
    end
    return intensities
end



"""
    instant_intensities(sf::StructureFactor, qs, mode; kwargs...)

Return ``𝒮(𝐪)`` intensities at wave vectors `qs`. The functionality is very
similar to [`intensities`](@ref), except the returned array has dimensions
identical to `qs`. If called on a `StructureFactor` with dynamical information,
i.e., ``𝒮(𝐪,ω)``, the ``ω`` information is integrated out.
"""
function instant_intensities(sf::StructureFactor, qs, mode; kwargs...)
    datadims = size(qs)
    ndims = length(datadims)
    vals = intensities(sf, qs, mode; static_warn=false, kwargs...)
    static_vals = sum(vals, dims=(ndims+1,))
    return reshape(static_vals, datadims)
end


"""
    connected_path(recip_vecs, qs::Vector, density)

Takes a list of wave vectors, `qs`, and builds an expanded list of wave vectors
that traces a path through the provided points. Also returned is a list of
marker indices corresponding to the input points. The `density` parameter is
given in samples per inverse Å.

Instead of `recip_vecs`, the first argument may be either a `StructureFactor` or
a `SpinWaveTheory`.
"""
function connected_path(recip_vecs, qs::Vector, density)
    @assert length(qs) >= 2 "The list `qs` should include at least two wavevectors."
    qs = Vec3.(qs)

    path = Vec3[]
    markers = Int[]
    for i in 1:length(qs)-1
        push!(markers, length(path)+1)
        q1, q2 = qs[i], qs[i+1]
        dist = norm(recip_vecs*(q1 - q2))
        npoints = round(Int, dist*density)
        for n in 1:npoints
            push!(path, (1 - (n-1)/npoints)*q1 + (n-1)*q2/npoints)
        end
    end
    push!(markers, length(path)+1)
    push!(path, qs[end])

    return (path, markers)
end
connected_path(sf::StructureFactor, qs::Vector, density) = connected_path(2π*inv(sf.crystal.latvecs)', qs, density)
connected_path(sw::SpinWaveTheory, qs::Vector, density) = connected_path(sw.recipvecs_chem, qs, density)


"""
    lorentzian(x, η) 

Returns ``η/(π(x^2 + η^2))``.
"""
lorentzian(x, η) = η/(π*(x^2 + η^2))

"""
    broaden_energy(sf::StructureFactor, vals, kernel::Function; negative_energies=false)

Performs a real-space convolution along the energy axis of an array of
intensities. Assumes the format of the intensities array corresponds to what
would be returned by [`intensities`](@ref). `kernel` must be a function that
takes two numbers: `kernel(ω, ω₀)`, where `ω` is a frequency, and `ω₀` is the
center frequency of the kernel. Sunny provides [`lorentzian`](@ref)
for the most common use case:

```
newvals = broaden_energy(sf, vals, (ω, ω₀) -> lorentzian(ω-ω₀, 0.2))
```
"""
function broaden_energy(sf::StructureFactor, is, kernel::Function; negative_energies=false)
    dims = size(is)
    ωvals = ωs(sf; negative_energies)
    out = zero(is)
    for (ω₀i, ω₀) in enumerate(ωvals)
        for (ωi, ω) in enumerate(ωvals)
            for qi in CartesianIndices(dims[1:end-1])
                out[qi,ωi] += is[qi,ω₀i]*kernel(ω, ω₀)
            end
        end
    end
    return out
end
