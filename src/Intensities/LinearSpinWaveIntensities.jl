"""
    intensities_broadened(swt::SpinWaveTheory, qs, ωvals, formula)

Computes the scattering intensities at each `(Q,ω)` according to Linear Spin Wave Theory
and the given intensity `formula`. The required `formula` must have a non-delta-function
kernel, e.g.:
    
    formula = intensity_formula(swt; kernel = lorentzian(0.05))

or else the intensity at `ωvals` which are not exactly on the dispersion curve can not
be calculated.

The intensity is computed at each wave vector in `qs` and each energy in `ωvals`.
The output will be an array with indices identical to `qs`, with the last index matching `ωvals`.

Note that `qs` is an array of wave vectors of arbitrary dimension. Each element
``q`` of `qs` must be a 3-vector in reciprocal lattice units. I.e., ``qᵢ`` is
given in ``2π/|aᵢ|`` with ``|aᵢ|`` the lattice constant of the chemical lattice.
"""
function intensities_broadened(swt::SpinWaveTheory, qs, ωvals, formula)
    qs = Vec3.(qs)
    nmodes = num_bands(swt)
    num_ω = length(ωvals)

    return_type = typeof(formula).parameters[1]
    if return_type <: BandStructure 
        # This only happens if the user set `kernel = delta_function_kernel`
        error("intensities_broadened: Can't compute broadened intensities without a finite-width kernel.\nTry: intensity_formula(...; kernel = lorentzian(0.05))")
    end

    is = zeros(return_type, size(qs)..., num_ω)

    # Compute the intensity at each (q,ω) pair
    for qidx in CartesianIndices(qs)
        intensity_as_function_of_ω = formula.calc_intensity(swt,qs[qidx])
        is[qidx,:] .= intensity_as_function_of_ω(ωvals)
    end

    return is
end

"""
    dispersions, intensities = intensities_bands(swt::SpinWaveTheory, qs; [formula])

Computes the scattering intensities at each energy band for each `q` in `qs`, according
to Linear Spin Wave Theory and the given intensity `formula`.
The optional `formula` must have a delta-function kernel, e.g.:
    
    formula = intensity_formula(swt; kernel = delta_function_kernel)

or else the bands will be broadened, and their intensity can not be computed.

The outputs will be arrays with indices identical to `qs`, with the last index
giving the band index. `dispersions` reports the energy of each band, while
`intensities` reports the scattering intensity.
"""
function intensities_bands(swt::SpinWaveTheory, qs; formula = intensity_formula(swt;kernel = nothing) :: SpinWaveIntensityFormula)
    if !isnothing(formula.kernel)
        # This is only triggered if the user has explicitly specified a formula with e.g. kT
        # corrections applied, but has not disabled the broadening kernel.
        error("intensities_bands: Can't compute band intensities if a broadening kernel is applied.\nTry intensity_formula(...; kernel = delta_function_kernel)")
    end

    qs = Vec3.(qs)
    nmodes = num_bands(swt)

    # Get the type parameter from the BandStructure
    return_type = typeof(formula).parameters[1].parameters[2]

    band_dispersions = zeros(Float64,length(qs),nmodes)
    band_intensities = zeros(return_type,length(qs),nmodes)
    for qidx in CartesianIndices(qs)
        band_structure = formula.calc_intensity(swt, qs[qidx])

        # Place the BandStructure at each point into its location in the array
        band_dispersions[qidx,:] .= band_structure.dispersion
        band_intensities[qidx,:] .= band_structure.intensity
    end
    return band_dispersions, band_intensities
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


