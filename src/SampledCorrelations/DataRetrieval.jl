# Takes a list of q points, converts into SampledCorrelation.data indices and
# corresponding exact wave vectors, and eliminates repeated elements.
function pruned_wave_vector_info(sc::SampledCorrelations, qs)

    # Round to the nearest wavevector and wrapped index
    Ls = size(sc.samplebuf)[2:4]
    ms = map(qs) do q 
        round.(Int, Ls .* q)
    end
    idcs = map(ms) do m
        CartesianIndex{3}(map(i -> mod(m[i], Ls[i])+1, (1, 2, 3)))
    end

    # Convert to absolute units (for form factors)
    qabs_rounded = map(m -> sc.crystal.recipvecs * (m ./ sc.sys_dims), ms)

    # List of "starting" pointers i where qabs_rounded[i-1] != qabs_rounded[i],
    # i.e., indices where the desired wave vector is distinct from the previous
    # one.
    starts = findall(i -> i == 1 || !isapprox(qabs_rounded[i-1], qabs_rounded[i]), eachindex(qabs_rounded))

    # Length of each run of repeated values
    counts = starts[2:end] - starts[1:end-1]
    append!(counts, length(idcs) - starts[end] + 1)

    # Remove contiguous repetitions
    qabs = qabs_rounded[starts]
    idcs = idcs[starts]

    return (; qabs, idcs, counts)
end


# Crude slow way to find the energy axis index closest to some given energy.
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

contains_dynamic_correlations(sc) = !isnan(sc.Δω)

# Documented under intensities function for LSWT. TODO: As a hack, this function
# is also being used as the back-end to intensities_static.
function intensities(sc::SampledCorrelations, qpts; energies, kernel=nothing, kT)
    if !isnothing(kT) && kT <= 0
        error("Positive `kT` required for classical-to-quantum corrections, or set `kT=nothing` to disable.")
    end
    if !isnothing(kernel)
        error("Kernel post-processing not yet available for `SampledCorrelations`.")
    end

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

    # Prepare memory and configuration variables for actual calculation
    qpts = Base.convert(AbstractQPoints, qpts)
    qs_reshaped = [to_reshaped_rlu(sc, q) for q in qpts.qs]

    ffs = sc.measure.formfactors[1, :] # FIXME
    intensities = zeros(eltype(sc.measure), isnan(sc.Δω) ? 1 : length(ωs), length(qpts.qs)) # N.B.: Inefficient indexing order to mimic LSWT
    q_idx_info = pruned_wave_vector_info(sc, qs_reshaped)
    crystal = @something sc.origin_crystal sc.crystal
    NCorr  = Val{size(sc.data, 1)}()
    # NPos = Val{size(sc.data, 2)}()
    NPos = Val{length(sc.crystal.positions)}()

    # Intensities calculation
    intensities_aux!(intensities, sc.data, sc.crystal, sc.positions, sc.measure.combiner, ffs, q_idx_info, ωidcs, NCorr, NPos)

    # Convert to a q-space density in original (not reshaped) RLU.
    intensities .*= det(sc.crystal.recipvecs) / det(crystal.recipvecs)

    # Post-processing steps for dynamical correlations 
    if contains_dynamic_correlations(sc) 
        # Convert time axis to a density.
        n_all_ω = size(sc.samplebuf, 6)
        intensities ./= (n_all_ω * sc.Δω)

        # Apply classical-to-quantum correspondence factor for finite kT
        if !isnothing(kT)
            # Equivalent to abs(ω/kT) * thermal_prefactor(ω; kT)
            c2q = [iszero(ω) ? 1 : abs((ω/kT) / (1 - exp(-ω/kT))) for ω in ωs]
            for i in axes(intensities, 2)
                intensities[:, i] .*= c2q
            end
        end
    end

    intensities = reshape(intensities, length(ωs), size(qpts.qs)...)

    return if contains_dynamic_correlations(sc) 
        Intensities(crystal, qpts, collect(ωs), intensities)
    else
        StaticIntensities(crystal, qpts, dropdims(intensities; dims=1))
    end
end

function intensities_aux!(intensities, data, crystal, positions, combiner, ff_atoms, q_idx_info, ωidcs, ::Val{NCorr}, ::Val{NPos}) where {NCorr, NPos}
    (; qabs, idcs, counts) = q_idx_info 
    (; recipvecs) = crystal 
    qidx = 1
    for (qabs, idx, count) in zip(qabs, idcs, counts)
        prefactors = prefactors_for_phase_averaging(qabs, recipvecs, view(positions, idx, :), ff_atoms, Val{NCorr}(), Val{NPos}())

        # Perform phase-averaging over all omega
        for (n, iω) in enumerate(ωidcs)
            elems = zero(SVector{NCorr, ComplexF64})
            for j in 1:NPos, i in 1:NPos
                elems += (prefactors[i] * conj(prefactors[j])) * SVector{NCorr}(view(data, :, i, j, idx, iω))
            end
            val = combiner(qabs, elems)
            intensities[n, qidx] = val
        end

        # Copy for repeated q-values
        for idx in qidx+1:qidx+count-1, n in axes(ωidcs, 1)
            intensities[n, idx] = intensities[n, qidx]
        end

        qidx += count
    end
end


function intensities_static(sc::SampledCorrelations, qpts; bounds = (-Inf, Inf), kT)
    ωs = available_energies(sc; negative_energies=true)
    ωidcs = findall(x -> bounds[1] <= x < bounds[2], ωs)
    if iszero(length(ωidcs))
        error("No information available within specified energy `bounds`. Try a larger interval.")
    end
    energies = sort(ωs[ωidcs])
    res = intensities(sc, qpts; kT, energies)
    data_new = dropdims(sum(res.data, dims=1), dims=1) * sc.Δω
    StaticIntensities(res.crystal, res.qpts, data_new)
end

function intensities_static(sc::SampledCorrelationsStatic, qpts)
    intensities(sc.parent, qpts; kT=nothing, energies=:available)
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


