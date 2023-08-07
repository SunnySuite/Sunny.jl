"""
    intensities_broadened(swt::SpinWaveTheory, ks, ωvals, formula)

Computes the scattering intensities at each `(Q,ω)` according to Linear Spin Wave Theory
and the given intensity `formula`. The required `formula` must have a non-delta-function
kernel, e.g.:
    
    formula = intensity_formula(swt; kernel = lorentzian(0.05))

or else the intensity at `ωvals` which are not exactly on the dispersion curve can not
be calculated.

The intensity is computed at each wave vector in `ks` and each energy in `ωvals`.
The output will be an array with indices identical to `ks`, with the last index matching `ωvals`.

Note that `ks` is an array of wave vectors of arbitrary dimension. Each element
``k`` of `ks` must be a 3-wavevector in absolute units.
"""
function intensities_broadened(swt::SpinWaveTheory, ks, ωvals, formula)
    ks = Vec3.(ks)
    nmodes = num_bands(swt)
    num_ω = length(ωvals)

    return_type = typeof(formula).parameters[1]
    if return_type <: BandStructure 
        # This only happens if the user sets `kernel = delta_function_kernel`
        error("intensities_broadened: Can't compute broadened intensities without a finite-width kernel.\nTry: intensity_formula(...; kernel = lorentzian(0.05))")
    end

    is = zeros(return_type, size(ks)..., num_ω)

    # Compute the intensity at each (k,ω) pair
    for kidx in CartesianIndices(ks)
        intensity_as_function_of_ω = formula.calc_intensity(swt,ks[kidx])
        is[kidx,:] .= intensity_as_function_of_ω(ωvals)
    end

    return is
end

"""
    dispersion, intensities = intensities_bands(swt::SpinWaveTheory, ks; [formula])

Computes the scattering intensities at each energy band for each momentum transfer `k` in `ks`, according
to Linear Spin Wave Theory and the given intensity `formula`.
The optional `formula` must have a delta-function kernel, e.g.:
    
    formula = intensity_formula(swt; kernel = delta_function_kernel)

or else the bands will be broadened, and their intensity can not be computed.

The outputs will be arrays with indices identical to `ks`, with the last index
giving the band index. `dispersions` reports the energy of each band, while
`intensities` reports the scattering intensity.
"""
function intensities_bands(swt::SpinWaveTheory, ks; formula = intensity_formula(swt;kernel = nothing) :: SpinWaveIntensityFormula)
    if !isnothing(formula.kernel)
        # This is only triggered if the user has explicitly specified a formula with e.g. kT
        # corrections applied, but has not disabled the broadening kernel.
        error("intensities_bands: Can't compute band intensities if a broadening kernel is applied.\nTry intensity_formula(...; kernel = delta_function_kernel)")
    end

    ks = Vec3.(ks)
    nmodes = num_bands(swt)

    # Get the type parameter from the BandStructure
    return_type = typeof(formula).parameters[1].parameters[2]

    band_dispersions = zeros(Float64,length(ks),nmodes)
    band_intensities = zeros(return_type,length(ks),nmodes)
    for kidx in CartesianIndices(ks)
        band_structure = formula.calc_intensity(swt, ks[kidx])

        # Place the BandStructure at each point into its location in the array
        band_dispersions[kidx,:] .= band_structure.dispersion
        band_intensities[kidx,:] .= band_structure.intensity
    end
    return band_dispersions, band_intensities
end

"""
    intensities_bin_centers(swt::SpinWaveTheory, params::BinningParameters, formula)

Computes the unpolarized inelastic neutron scattering intensities given a
`SpinWaveTheory`, histogram described by its `BinningParameters`, and
an [`intensity_formula`](@ref) with finite kernel width, e.g.:

    formula = intensity_formula(swt; kernel = lorentzian(0.05))

Note that this method only calculates the intensity at the bin centers--it doesn't
integrate over the bins in any way. The output will be the same shape as if it were
histogrammed data.
"""
function intensities_bin_centers(swt::SpinWaveTheory, params::BinningParameters, formula::SpinWaveIntensityFormula)
    if any(params.covectors[1:3,4] .!= 0.) || any(params.covectors[4,1:3] .!= 0.)
      error("Complicated binning parameters not supported by intensities_bin_centers")
    end

    bin_centers = axes_bincenters(params)


    q_params = copy(params)
    bin_rlu_as_absolute_units!(q_params,swt)
    # coords = q_covectors * (q,ω)
    coords_to_q = inv(q_params.covectors[1:3,1:3])

    is = zeros(Float64,params.numbins...)

    # Loop over qs
    for ci in CartesianIndices(params.numbins.data[1:3])
        x_center = bin_centers[1][ci[1]]
        y_center = bin_centers[2][ci[2]]
        z_center = bin_centers[3][ci[3]]

        q = SVector{3}(coords_to_q * [x_center;y_center;z_center])
        ωvals = bin_centers[4]

        intensity_as_function_of_ω = formula.calc_intensity(swt,q)
        is[ci,:] .= intensity_as_function_of_ω(ωvals)
    end
    is
end


