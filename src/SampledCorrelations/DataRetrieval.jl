function find_idx_of_nearest_fft_energy(ref, val)
    for i in axes(ref[1:end-1], 1)
        x1, x2 = ref[i], ref[i+1]
        if x1 <= val <= x2 
            if abs(x1 - val) < abs(x2 - val)
                return i
            else
                return i+1
            end
        end
    end
    # Deal with edge case arising due to FFT index ordering
    if ref[end] <= val <= 0.0
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
function intensities_old(sc::SampledCorrelations, qpts; energies, kernel=nothing, formfactors=nothing, kT)
    if isnan(sc.Δω) && energies != :available
        error("`SampledCorrelations` only contains static information. No energy information is available.")
    end
    if !isnothing(kernel)
        error("Kernel post-processing not yet available for `SampledCorrelations`.")
    end

    (; measure) = sc
    interp = NoInterp()
    qpts = Base.convert(AbstractQPoints, qpts)
    qs = view(qpts.qs, :)
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

    # q-points in RLU for the reshaped crystal
    q_targets = to_reshaped_rlu.(Ref(sc), qs)

    # Preallocation
    nω = isnan(sc.Δω) ? 1 : length(ωs)
    intensities = zeros(IntensitiesType, nω, length(qs)) # N.B.: Inefficient indexing order to mimic LSWT
    local_intensities = zeros(IntensitiesType, ninterp(interp)) 

    # Stencil and interpolation precalculation for q-space
    li_intensities = LinearIndices(intensities)
    ci_targets = eachindex(q_targets)
    m_targets = [mod.(sc.latsize .* q_target, 1) for q_target in q_targets]
    index_info = (; li_intensities, ci_targets, m_targets)
    stencil_info = pruned_stencil_info(sc, qpts.qs, interp)
    NCorr  = Val{size(sc.data, 1)}()
    NAtoms = Val{size(sc.data, 2)}()

    intensities_interpolated_old!(intensities, local_intensities, sc, measure.combiner, ff_atoms, ωidcs, stencil_info, index_info, interp, NCorr, NAtoms)

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

    intensities = reshape(intensities, nω, size(qpts.qs)...)
    return if !isnan(sc.Δω)
        Intensities(crystal, qpts, collect(ωs), intensities)
    else
        InstantIntensities(crystal, qpts, dropdims(intensities; dims=1))
    end
end

function intensities_interpolated_old!(intensities, local_intensities, sc, combiner, ff_atoms, ωidcs, stencil_info, index_info, interp, ::Val{NCorr}, ::Val{NAtoms}) where {NCorr, NAtoms}
    (; qabs_all, idcs_all, counts) = stencil_info 
    (; li_intensities, ci_targets, m_targets) = index_info

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
                    Val(NCorr),
                    Val(NAtoms)
                )
                local_intensities[n] = combiner(qabs[n], correlations)
            end

            # Perform interpolations for all requested qs using stencil
            # intensities.
            for _ in 1:numrepeats
                iq += 1
                idx = li_intensities[CartesianIndex(n, ci_targets[iq])]
                intensities[idx] = interpolate(sc, m_targets[iq], local_intensities, interp) 
            end
        end
    end
end


"""
    intensities_instant(sc::SampledCorrelations, qpts; kernel=nothing, formfactors=nothing, kT)

Calculate the instant (equal-time) correlations for a set of q-points in
reciprocal space.

This calculation may be performed on a [SampledCorrelations](@ref) regardless of
whether it contains only static correlations or dynamic correlations. If the
`SampledCorrelation` does contain correlations from time-evolved trajectories,
then `intensities_instant` returns the energy-integrated correlations. If, in
addition, `kT` is given a numerical value, the dynamic correlations will be
multiplied by the "classical-to-quantum" correspondence factor before energy
integration:

```math
βω [1 + n_B(ω)],
```

where ``n_B(ω) = 1/(exp(βω) - 1)`` is the Bose function and ``β=1/(k_B T)``.

If `kT` is set to `nothing` (the default behavior), this correction will not be
applied. Note that temperature-dependent corrections are only available when the
`SampledCorrelations` contains dynamic correlations.
"""
function intensities_instant(sc::SampledCorrelations, qpts; kernel=nothing, formfactors=nothing, kT=nothing)
    return if isnan(sc.Δω)
        if !isnothing(kT) 
            error("Temperature corrections unavailable if `SampledCorrelations` does not contain dynamic correlations. Do not set `kT` value.")
        end
        intensities(sc, qpts; kernel, formfactors, kT, energies=:available) # Returns an InstantIntensities
    else
        res = intensities(sc, qpts; kernel, formfactors, kT, energies=:available_with_negative) # Returns a Intensities
        data_new = dropdims(sum(res.data, dims=1), dims=1)
        InstantIntensities(res.crystal, res.qpts, data_new)
    end
end


function classical_to_quantum(ω, kT)
    if ω > 0
        ω/(kT*(1 - exp(-ω/kT)))
    elseif iszero(ω)
        1.0
    else
        -ω*exp(ω/kT)/(kT*(1 - exp(ω/kT)))
    end
end


#=
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
=#






################################################################################
function intensities(sc::SampledCorrelations, qpts; energies, kernel=nothing, formfactors=nothing, kT)
    if isnan(sc.Δω) && energies != :available
        error("`SampledCorrelations` only contains static information. No energy information is available.")
    end
    if !isnothing(kernel)
        error("Kernel post-processing not yet available for `SampledCorrelations`.")
    end

    (; measure) = sc
    interp = NoInterp()
    qpts = Base.convert(AbstractQPoints, qpts)
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

    # Stencil and interpolation precalculation for q-space
    stencil_info = full_stencil_info(sc, qpts.qs)
    NCorr  = Val{size(sc.data, 1)}()
    NAtoms = Val{size(sc.data, 2)}()

    intensities_interpolated!(intensities, sc.data, sc.crystal, measure, ff_atoms, ωidcs, stencil_info, NCorr, NAtoms)

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
        Intensities(crystal, qpts, collect(ωs), intensities)
    else
        InstantIntensities(crystal, qpts, dropdims(intensities; dims=1))
    end
end

function intensities_interpolated!(intensities, data, crystal, measure::MeasureSpec{Op, F, Ret}, ff_atoms, ωidcs, stencil_info, ::Val{NCorr}, ::Val{NAtoms}) where {Op, F, Ret, NCorr, NAtoms}
    (; qabs_rounded, idcs) = stencil_info 
    for (n, (qabs, idx)) in enumerate(zip(qabs_rounded, idcs))
        prefactors = prefactors_for_phase_averaging(qabs, crystal, ff_atoms, Val(NCorr), Val(NAtoms))
        for iω in ωidcs
            elems = zero(MVector{NCorr,ComplexF64})
            for j in 1:NAtoms, i in 1:NAtoms
                elems .+= (prefactors[i] * conj(prefactors[j])) .* view(data, :, i, j, idx, iω)
            end
            intensities[iω, n] = measure.combiner(qabs, SVector{NCorr, ComplexF64}(elems))
        end
    end
end



