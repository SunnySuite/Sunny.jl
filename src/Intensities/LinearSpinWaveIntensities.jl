"""
    intensities_broadened(swt::SpinWaveTheory, qs, ωs, formula)

Computes the scattering intensities at each `(q,ω)` according to Linear Spin
Wave Theory and the given intensity `formula`. The required `formula` must have
a non-delta-function kernel, e.g.:
    
    formula = intensity_formula(swt, :perp; kernel = lorentzian(0.05))

or else the intensity at `ωs` which are not exactly on the dispersion curve can
not be calculated.

The intensity is computed at each wave vector in `qs` and each energy in `ωs`.
The output will be an array with indices identical to `qs`, with the last index
matching `ωs`.

Note that `qs` is an array of wave vectors of arbitrary dimension. Each element
``q`` of `qs` must be a 3-wavevector in absolute units.
"""
function intensities_broadened(swt::SpinWaveTheory, qs, ωs, formula)
    qs = Vec3.(qs)
    num_ω = length(ωs)

    return_type = typeof(formula).parameters[1]
    if return_type <: BandStructure 
        # This only happens if the user sets `kernel = delta_function_kernel`
        error("""intensities_broadened: Can't compute broadened intensities without a finite-width kernel.
                 Try: intensity_formula(...; kernel=lorentzian(fwhm=...))
                 """)
    end

    is = zeros(return_type, size(qs)..., num_ω)

    # Compute the intensity at each (k,ω) pair
    for qidx in CartesianIndices(qs)
        intensity_as_function_of_ω = formula.calc_intensity(swt, qs[qidx])
        is[qidx,:] .= intensity_as_function_of_ω(ωs)
    end

    return is
end

"""
    dispersion, intensities = intensities_bands(swt::SpinWaveTheory, qs, formula::SpinWaveIntensityFormula)

Computes the scattering intensities at each energy band for each momentum
transfer `q` in `qs`, according to linear spin wave theory and the given
intensity `formula`. The `formula` must have a delta-function kernel, e.g.:
    
    formula = intensity_formula(swt, :perp, formula; kernel=delta_function_kernel)

or else the bands will be broadened, and their intensity can not be computed.

The outputs will be arrays with indices identical to `qs`, with the last index
giving the band index. `dispersions` reports the energy of each band, while
`intensities` reports the scattering intensity.
"""
function intensities_bands(swt::SpinWaveTheory, qs, formula::SpinWaveIntensityFormula)
    if !isnothing(formula.kernel)
        # This is only triggered if the user has explicitly specified a formula with e.g. kT
        # corrections applied, but has not disabled the broadening kernel.
        error("""intensities_bands: Can't compute band intensities if a broadening kernel is applied.
                 Try intensity_formula(...; kernel = delta_function_kernel)
                 """)
    end

    qs = Vec3.(qs)
    nmodes = nbands(swt)

    # Get the type parameter from the BandStructure
    return_type = typeof(formula).parameters[1].parameters[2]

    band_dispersions = zeros(Float64,size(qs)...,nmodes)
    band_intensities = zeros(return_type,size(qs)...,nmodes)
    for qidx in CartesianIndices(qs)
        band_structure = formula.calc_intensity(swt, qs[qidx])

        # Place the BandStructure at each point into its location in the array
        band_dispersions[qidx,:] .= band_structure.dispersion
        band_intensities[qidx,:] .= band_structure.intensity
    end
    return band_dispersions, band_intensities
end



"""
    intensities_bin_centers(swt::SpinWaveTheory, params::BinningParameters, formula)

Computes the unpolarized inelastic neutron scattering intensities given a
`SpinWaveTheory`, histogram described by its `BinningParameters`, and
an [`intensity_formula`](@ref) with finite kernel width, e.g.:

    formula = intensity_formula(swt, :perp; kernel = lorentzian(0.05))

Note that this method only calculates the intensity at the bin centers--it doesn't
integrate over the bins in any way. The output will be the same shape as if it were
histogrammed data.
"""
function intensities_bin_centers(swt::SpinWaveTheory, params::BinningParameters, formula)
    if any(params.covectors[1:3,4] .!= 0.) || any(params.covectors[4,1:3] .!= 0.)
      error("Complicated binning parameters not supported by intensities_bin_centers")
    end

    bin_centers = axes_bincenters(params)


    # coords = covectors * (q,ω)
    coords_to_q = inv(params.covectors[1:3,1:3])

    is = zeros(Float64,params.numbins...)

    # Loop over qs
    for ci in CartesianIndices(params.numbins.data[1:3])
        x_center = bin_centers[1][ci[1]]
        y_center = bin_centers[2][ci[2]]
        z_center = bin_centers[3][ci[3]]

        q = SVector{3}(coords_to_q * [x_center; y_center; z_center])
        ωs = bin_centers[4]

        intensity_as_function_of_ω = formula.calc_intensity(swt, q)
        is[ci,:] .= intensity_as_function_of_ω(ωs)
    end
    is
end

function intensities_bin_multisample(swt::SpinWaveTheory, hist_params::BinningParameters, msaa_strategy, energy_msaa_strategy, formula)
    if any(hist_params.covectors[1:3,4] .!= 0.) || any(hist_params.covectors[4,1:3] .!= 0.)
      error("Complicated binning parameters not supported by intensities_bin_centers")
    end

    bin_edges = axes_binedges(hist_params)
    bin_diagonal_vector = reshape(hist_params.binwidth[1:3],3,1)

    # coords = hist_covectors * (q,ω)
    coords_to_q = inv(hist_params.covectors[1:3,1:3])

    is = zeros(Float64, hist_params.numbins...)
    counts = zeros(Float64, hist_params.numbins...)

    # Visits each bin in hist_params
    for ci in CartesianIndices(hist_params.numbins.data[1:3])
        x_lower = bin_edges[1][ci[1]]
        y_lower = bin_edges[2][ci[2]]
        z_lower = bin_edges[3][ci[3]]

        for msaa_location in msaa_strategy
            coords_xyz = [x_lower; y_lower; z_lower] .+ bin_diagonal_vector .* msaa_location
            q = SVector{3}(coords_to_q * coords_xyz)

            if isnothing(formula.kernel) # Need to handle BandStructure
                if !isempty(energy_msaa_strategy)
                    error("Attempted to multisample in energy with delta function kernel")
                end
                band_structure = formula.calc_intensity(swt,q)
                energy_bins = 1 .+ floor.(Int64,(band_structure.dispersion .- hist_params.binstart[4]) ./ hist_params.binwidth[4])
                for (i,bin_i) in enumerate(energy_bins)
                    if 1 <= bin_i <= hist_params.numbins[4]
                        is[ci,bin_i] += band_structure.intensity[i]
                        counts[ci,bin_i] += 1
                    end
                end
            else # Broadening is done by the formula
                intensity_as_function_of_ω = formula.calc_intensity(swt, q)
                energy_lower_edges = bin_edges[4][1:end-1]
                for (i,energy_sample) in enumerate(energy_msaa_strategy)
                    energy_sample
                    bin_edges[4]
                    ω_this_sample = energy_lower_edges .+ hist_params.binwidth[4] .* energy_sample
                    view(is,ci,:) .+= intensity_as_function_of_ω(ω_this_sample)
                end
                view(counts,ci,:) .+= length(energy_msaa_strategy)
            end
        end
    end
    is, counts
end


