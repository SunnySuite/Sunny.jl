## TODO: Cleanup and docstring
## TODO: Add proper treatment for negative qs?
function intensities(sc::SampledCorrelations, qpts; measure=nothing, negative_energies=false, formfactors=nothing, interp=NoInterp())
    measure = !isnothing(measure) ? measure : sc.measure # TODO: Add checks to see if override is legit
    qpts = Base.convert(AbstractQPoints, qpts)

    ff_atoms = propagate_form_factors_to_atoms(formfactors, sc.crystal)
    # corr_ix = 1:length(measure.corr_pairs)  # This assumes all correlations are kept

    # Type stability for phase_averaged_elements
    IntensitiesType = eltype(measure)
    NInterp = 1 # Generalize this
    NCorr = Val(size(sc.data, 1))
    NAtoms = Val(size(sc.data, 2))

    # Interpret q points in terms of original crystal. 
    q_targets = if !isnothing(sc.origin_crystal)
        convert = sc.crystal.recipvecs \ sc.origin_crystal.recipvecs
        [convert * Vec3(q) for q in qpts.qs]
    else
        qpts.qs
    end

    ωvals = available_energies(sc; negative_energies)
    intensities = zeros(IntensitiesType, length(ωvals), length(qpts.qs))

    li_intensities = LinearIndices(intensities)
    ci_targets = CartesianIndices(q_targets)
    m_targets = [mod.(sc.latsize .* q_target, 1) for q_target in q_targets]

    (; qabs_all, idcs_all, counts) = pruned_stencil_info(sc, qpts.qs, interp) 
    local_intensities = zeros(IntensitiesType, NInterp) 

    for iω in eachindex(ωvals)
        iq = 0
        for (qabs, idcs, numrepeats) in zip(qabs_all, idcs_all, counts)

            # Pull out nearest intensities that are necessary for any interpolation
            for n in 1:NInterp
                correlations = phase_averaged_elements(view(sc.data, :, :, :, idcs[n], iω), qabs[n], sc.crystal, ff_atoms, NCorr, NAtoms)
                local_intensities[n] = measure.combiner(qabs[n], correlations)
            end

            # Perform interpolations 
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(iω, ci_targets[iq])]
                intensities[idx] = interpolated_intensity(sc, m_targets[iq], local_intensities, interp) 
            end
        end
    end

    # This converts the time axis to a density. TODO: Why not do this with the
    # definition of the FFT normalization?
    if !isnan(sc.Δω)
        n_all_ω = size(sc.samplebuf, 6)
        intensities ./= (n_all_ω * sc.Δω)
    end 

    crystal = !isnothing(sc.origin_crystal) ? sc.origin_crystal : sc.crystal
    return BroadenedIntensities(crystal, qpts, ωvals, intensities)
end

## TODO: Uncomment after decision made about instant correlations.
# """
#     instant_intensities_interpolated(sc::SampledCorrelations, qs, formula::ClassicalIntensityFormula; kwargs...)
# 
# Return ``𝒮(𝐪)`` intensities at wave vectors `qs`. The functionality is very
# similar to [`intensities_interpolated`](@ref), except the returned array has dimensions
# identical to `qs`. If called on a `SampledCorrelations` with dynamical information,
# i.e., ``𝒮(𝐪,ω)``, the ``ω`` information is integrated out.
# """
# function instant_intensities_interpolated(sc::SampledCorrelations, qs, formula; kwargs...)
#     datadims = size(qs)
#     ndims = length(datadims)
#     vals = intensities_interpolated(sc, qs, formula; instantaneous_warning=false, kwargs...)
#     static_vals = sum(vals, dims=(ndims+1,))
#     return reshape(static_vals, datadims)
# end


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
