"""
    available_wave_vectors(sc::SampledCorrelations; counts=(1,1,1))

Returns the grid of wave vectors for which `sc` contains exact values.
Optionally extend by a given number of `counts` along each grid axis. If the
system was not reshaped, then the number of Brillouin zones included is
`prod(counts)`.
"""
function available_wave_vectors(sc::SampledCorrelations; counts=(1,1,1))
    Ls = sc.sys_dims
    offsets = map(L -> isodd(L) ? 1 : 0, Ls)
    up = Ls .* counts
    hi = map(L -> L - div(L, 2), up) .- offsets
    lo = map(L -> 1 - div(L, 2), up) .- offsets

    orig_crystal = @something sc.origin_crystal sc.crystal
    convert = orig_crystal.recipvecs \ sc.crystal.recipvecs
    return [convert * Vec3(lx/Ls[1], ly/Ls[2], lz/Ls[3]) for lx in lo[1]:hi[1], ly in lo[2]:hi[2], lz in lo[3]:hi[3]]
end

"""
    available_energies(sc::SampledCorrelations; negative_energies=false)

Return the ω values for the energy index of a `SampledCorrelations`. By default,
only returns values for non-negative energies, which corresponds to the default
output of `intensities`. Set `negative_energies` to true to retrieve all ω
values.
"""
function available_energies(sc::SampledCorrelations; negative_energies=false)
    isnan(sc.Δω) && return NaN

    n_all_ω = size(sc.data, 7)
    n_non_neg_ω = div(n_all_ω, 2) + 1
    ωvals = collect(FFTW.fftfreq(n_all_ω, n_all_ω * sc.Δω))
    ωvals[n_non_neg_ω] *= -1  # Adjust for FFTW convention (which is largest frequency negative)
    return negative_energies ? ωvals : ωvals[1:n_non_neg_ω]
end
