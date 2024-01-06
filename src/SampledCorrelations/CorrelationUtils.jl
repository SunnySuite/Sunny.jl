"""
    available_wave_vectors(sc::SampledCorrelations; bzsize=(1,1,1))

Returns all wave vectors for which `sc` contains exact values. `bsize` specifies
the number of Brillouin zones to be included.
"""
function available_wave_vectors(sc::SampledCorrelations; bzsize=(1,1,1))
    Ls = size(sc.samplebuf)[2:4]  # If we had a sys, would use latsize
    offsets = map(L -> isodd(L) ? 1 : 0, Ls)
    up = Ls .* bzsize
    hi = map(L -> L - div(L, 2), up) .- offsets
    lo = map(L -> 1 - div(L, 2), up) .- offsets
    qs = zeros(Vec3, up...)
    for (k, lz) in enumerate(lo[3]:hi[3]), (j, ly) in enumerate(lo[2]:hi[2]), (i, lx) in enumerate(lo[1]:hi[1])
        qs[i,j,k] = Vec3(lx/Ls[1], ly/Ls[2], lz/Ls[3])

        # If the crystal has been reshaped, convert all wavevectors from RLU in the
        # the reshaped crystal to RLU in the original crystal
        if !isnothing(sc.origin_crystal)
            convert = sc.origin_crystal.recipvecs \ sc.crystal.recipvecs
            qs = [convert * q for q in qs]
        end
    end
    return qs
end

"""
    available_energies(sc::SampledCorrelations; negative_energies=false)

Return the ω values for the energy index of a `SampledCorrelations`. By default,
only returns values for non-negative energies, which corresponds to the default
output of `intensities`. Set `negative_energies` to true to retrieve all ω
values.
"""
function available_energies(sc::SampledCorrelations; negative_energies=false)
    isnan(sc.Δω) && (return NaN)

    nω = size(sc.data, 7)
    hω = div(nω, 2) + 1
    ωvals = FFTW.fftfreq(nω, nω * sc.Δω)
    return negative_energies ? ωvals : ωvals[1:hω]
end