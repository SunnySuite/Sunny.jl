function find_idx_of_nearest_fft_energy(ref, val)
    for i in axes(ref[1:end-1], 1)
        x1, x2 = ref[i], ref[i+1]
        if x1 <= val < x2 
            if abs(x1 - val) < abs(x2 - val)
                return i
            else
                return i+1
            end
        end
    end
    # Deal with edge case arising due to FFT index ordering
    if ref[end] <= val < 0.0
        if abs(ref[end] - val) <= abs(val)
            return length(ref)
        else
            return 1
        end
    end
    error("Value does not lie in bounds of reference list.")
end

# If the user specifies an energy list, round to the nearest available energies
# and give the corresponding indices into the raw data. This is fairly
# inefficient, though the cost is likely trivial next to the rest of the
# computation. Since this is an edge case (user typically expected to choose
# :available or :available_with_negative), not spending time on optimization
# now.
function rounded_energy_information(sc, energies)
    ωvals = available_energies(sc; negative_energies = true)
    energies_sorted = sort(energies)
    @assert all(x -> x ==true, energies .== energies_sorted) "Specified energies must be an ordered list."
    @assert minimum(energies) >= minimum(ωvals) && maximum(energies) <= maximum(ωvals) "Specified energies includes values for which there is no available data."
    ωidcs = map(val -> find_idx_of_nearest_fft_energy(ωvals, val), energies)
    return ωvals[ωidcs], ωidcs
end

# Documented in same location as LSWT `intensities` function.
function intensities(sc::SampledCorrelations, qpts; energies, kernel=nothing, formfactors=nothing, kT)
    if isnan(sc.Δω) && energies != :available
        error("`SampledCorrelations` only contains static information. No energy information is available.")
    end
    if !isnothing(kernel)
        error("Kernel post-processing not yet available for `SampledCorrelations`.")
    end

    (; measure) = sc
    interp = NoInterp()
    qpts = Base.convert(AbstractQPoints, qpts[:])
    ff_atoms = propagate_form_factors_to_atoms(formfactors, sc.crystal)
    IntensitiesType = eltype(measure)
    crystal = !isnothing(sc.origin_crystal) ? sc.origin_crystal : sc.crystal

    # Determine energy information
    (ωs, ωidcs) = if energies == :available
        ωs = available_energies(sc; negative_energies=false)
        (ωs, axes(ωs, 1))
    elseif energies == :available_with_negative
        ωs = available_energies(sc; negative_energies=true)
        (ωs, axes(ωs, 1))
    else
        rounded_energy_information(sc, energies)
    end
    

    # Interpret q points in terms of original crystal. 
    q_targets = if !isnothing(sc.origin_crystal)
        convert = sc.crystal.recipvecs \ sc.origin_crystal.recipvecs
        [convert * Vec3(q) for q in qpts.qs]
    else
        qpts.qs
    end

    # Preallocation
    intensities = zeros(IntensitiesType, isnan(sc.Δω) ? 1 : length(ωs), length(qpts.qs)) # N.B.: Inefficient indexing order to mimic LSWT
    local_intensities = zeros(IntensitiesType, ninterp(interp)) 

    # Stencil and interpolation precalculation for q-space
    li_intensities = LinearIndices(intensities)
    ci_targets = CartesianIndices(q_targets)
    m_targets = [mod.(sc.latsize .* q_target, 1) for q_target in q_targets]
    (; qabs_all, idcs_all, counts) = pruned_stencil_info(sc, qpts.qs, interp) 

    # Calculate the interpolated intensities, avoiding repeated calls to
    # phase_averaged_elements. This is the computational expensive portion.
    for (n, iω) in enumerate(ωidcs)
        iq = 0
        for (qabs, idcs, numrepeats) in zip(qabs_all, idcs_all, counts)

            # Pull out nearest intensities formling the "stencil" used for
            # interpolation.
            for n in 1:ninterp(interp)
                correlations = phase_averaged_elements(
                    view(sc.data, :, :, :, idcs[n], iω), 
                    qabs[n], 
                    sc.crystal, 
                    ff_atoms, 
                    Val(size(sc.data, 1)), 
                    Val(size(sc.data, 2))
                )
                local_intensities[n] = measure.combiner(qabs[n], correlations)
            end

            # Perform interpolations for all requested qs using stencil
            # intensities.
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(n, ci_targets[iq])]
                intensities[idx] = interpolated_intensity(sc, m_targets[iq], local_intensities, interp) 
            end
        end
    end

    # Processing steps that depend on whether instant or dynamic correlations.
    if !isnan(sc.Δω)
        # Convert time axis to a density.
        # TODO: Why not do this with the definition of the FFT normalization?
        n_all_ω = size(sc.samplebuf, 6)
        intensities ./= (n_all_ω * sc.Δω)

        # Apply classical-to-quantum correspondence factor if temperature given.
        if !isnothing(kT) 
            c2q = classical_to_quantum.(ωs, kT)
            for i in axes(intensities, 2)
                intensities[:,i] .*= c2q
            end
        end
    else
        # If temperature is given for a SampledCorrelations with only
        # instantaneous data, throw an error. 
        !isnothing(kT) && error("Given `SampledCorrelations` only contains instant correlation data. Temperature corrections only available with dynamical correlation info.")
    end

    return if !isnan(sc.Δω)
        BroadenedIntensities(crystal, qpts, collect(ωs), intensities)
    else
        InstantIntensities(crystal, qpts, reshape(intensities, size(intensities, 2)))
    end
end


"""
    intensities_instant(sc::SampledCorrelations, qpts; kernel=nothing, formfactors=nothing, kT)

Calculate the instant (equal-time) correlations for a set of q-points in
reciprocal space.

This calculation may be performed on a `SampledCorrelations` regardless of
whether it contains dynamic correlation information or not. If the
`SampledCorrelation` does contain correlations from time-evolved trajectories,
the `intensities_instant` returns the energy-integrated correlations. If, in
addition, `kT` is given a numerical value, the dynamic correlations will be
multiplied by the "classical-to-quantum" correspondence factor before energy
integration. If `kT` is set to `nothing`, this correction will not be applied.
"""
# ```math
# \frac{\omega}{k_{\rm{B}}T}\left[1 + n_{\rm{B}}\left(\omega/T\right)\right].
# ```
# 
# where $n_{\rm{B}} = \left(e^{\omega/k_{\rm{T}}T}\right - 1)^{-1}$ is the Bose function.
function intensities_instant(sc::SampledCorrelations, qpts; kernel=nothing, formfactors=nothing, kT)
    return if isnan(sc.Δω)
        if !isnothing(kT) 
            error("Temperature corrections unavailable if `SampledCorrelations` does not contain dynamic correlations. Do not set `kT` value.")
        end
        intensities(sc, qpts; kernel, formfactors, kT, energies=:available) # Returns an InstantIntensities
    else
        is = intensities(sc, qpts; kernel, formfactors, kT, energies=:available_with_negative) # Returns a BroadenedIntensities
        data_new = reshape(sum(is.data, dims=(1,)), size(is.data, 2))
        InstantIntensities(is.crystal, is.qpts, data_new)
    end
end


function classical_to_quantum(ω, kT)
    if kT == Inf
        return 1.0
    end
    if ω > 0
        ω/(kT*(1 - exp(-ω/kT)))
    elseif iszero(ω)
        1.0
    else
        -ω*exp(ω/kT)/(kT*(1 - exp(ω/kT)))
    end
end

"""
    gaussian(; {fwhm, σ})

Returns the function `exp(-x^2/2σ^2) / √(2π*σ^2)`. Exactly one of `fwhm` or `σ`
must be specified, where `fwhm = (2.355...) * σ` denotes the full width at half
maximum.
"""
function gaussian06(; fwhm=nothing, σ=nothing)
    if sum(.!isnothing.((fwhm, σ))) != 1
        error("Exactly one of `fwhm` and `σ` must be specified.")
    end
    σ = Float64(@something σ (fwhm/2√(2log(2))))
    return x -> exp(-x^2/2σ^2) / √(2π*σ^2)
end


"""
    integrated_gaussian(; {fwhm, σ}) 

Returns the function `erf(x/√2σ)/2`, which is the integral of [`gaussian`](@ref)
over the range ``[0, x]``. Exactly one of `fwhm` or `σ` must be specified, where
`fwhm = (2.355...) * σ` denotes the full width at half maximum. Intended for use
with [`intensities_binned`](@ref).
"""
function integrated_gaussian(; fwhm=nothing, σ=nothing)
    if sum(.!isnothing.((fwhm, σ))) != 1
        error("Exactly one of `fwhm` and `σ` must be specified.")
    end
    σ = Float64(@something σ (fwhm/2√(2log(2))))
    return x -> erf(x/√2σ)/2
end

"""
    lorentzian(; fwhm)

Returns the function `(Γ/2) / (π*(x^2+(Γ/2)^2))` where `Γ = fwhm` is the full
width at half maximum.
"""
function lorentzian06(; fwhm)
    Γ = fwhm
    return x -> (Γ/2) / (π*(x^2+(Γ/2)^2))
end

"""
    integrated_lorentzian(; fwhm) 

Returns the function `atan(2x/Γ)/π`, which is the integral of
[`lorentzian`](@ref) over the range ``[0, x]``, where `Γ = fwhm` is the full
width at half maximum. Intended for use with [`intensities_binned`](@ref).
"""
function integrated_lorentzian(; fwhm)
    Γ = fwhm
    return x -> atan(2x/Γ)/π
end


"""
    broaden_energy(sc::SampledCorrelations, vals, kernel::Function; negative_energies=false)

Performs a real-space convolution along the energy axis of an array of
intensities. Assumes the format of the intensities array corresponds to what
would be returned by [`intensities_interpolated`](@ref). `kernel` must be a function that
takes two numbers: `kernel(ω, ω₀)`, where `ω` is a frequency, and `ω₀` is the
center frequency of the kernel. Sunny provides [`lorentzian`](@ref)
for the most common use case:

```
newvals = broaden_energy(sc, vals, (ω, ω₀) -> lorentzian06(fwhm=0.2)(ω-ω₀))
```
"""
function broaden_energy(sc::SampledCorrelations, is, kernel::Function; negative_energies=false)
    dims = size(is)
    ωvals = available_energies(sc; negative_energies)
    out = zero(is)
    for (ω₀i, ω₀) in enumerate(ωvals)
        for (ωi, ω) in enumerate(ωvals)
            for qi in CartesianIndices(dims[1:end-1])
                out[qi,ωi] += is[qi,ω₀i]*kernel(ω, ω₀)*sc.Δω
            end
        end
    end
    return out
end
