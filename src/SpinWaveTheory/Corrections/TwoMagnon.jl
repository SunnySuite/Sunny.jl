"""
    intensities_two_magnon(swt::SpinWaveTheory, qpts; energies, kernel, grid, bin_width=nothing)

Calculates the two-magnon contribution to the dynamical spin structure factor at
temperature ``T = 0``. This intensity appears at sub-leading order in ``1/s``,
and is not captured by [`intensities`](@ref).

Whereas linear spin wave theory yields sharp bands from the part of the spin
operator that is _linear_ in Holstein-Primakoff bosons, the _quadratic_ part
creates two quasi-particles at once. Momentum conservation leaves their
individual momenta free, so the result is a continuum,

```math
𝒮(𝐪, ω) = \\frac{1}{2} ∫d𝐤 \\, Σ_{n₁ n₂} |A(𝐤, n₁; 𝐪-𝐤, n₂)|^2
          δ(ω - ε_{𝐤 n₁} - ε_{𝐪-𝐤, n₂}),
```

where the ``δ``-function is replaced by the provided broadening `kernel`. The
integral runs over the first magnetic Brillouin zone, and is performed on a
uniform `grid` of the given dimensions.

Rather than broadening each sampled pair energy ``ε_{𝐤n₁} + ε_{𝐪-𝐤,n₂}``
separately, the weight ``|A|^2`` is accumulated into bins of that energy of width
`bin_width`, and the `kernel` is applied to the bins. Total weight is preserved
exactly, and the shape carries an error of order `(bin_width / kernel width)²`.
By default `bin_width` is a small fraction of the full width at half maximum of
the `kernel`, putting that error well below the error of the `grid`.
"""
function intensities_two_magnon(swt::SpinWaveTheory, qpts; energies, kernel::AbstractBroadening, grid, bin_width=nothing)
    check_corrections_supported(swt)

    (; sys, measure, data) = swt
    num_observables(measure) == 0 && error("No observables! Construct SpinWaveTheory with a `measure` argument.")
    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")
    if isnothing(bin_width)
        isa(kernel, Broadening) && !isnan(kernel.fwhm) || error("Keyword `bin_width` is required for a kernel of unknown width.")
        bin_width = kernel.fwhm / 32
    end

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)

    @assert sys.dims == (1, 1, 1)
    Na = nsites(sys)
    Ncells = Na / natoms(cryst)
    L = nbands(swt)
    Nobs = num_observables(measure)
    Ncorr = num_correlations(measure)

    H = zeros(ComplexF64, 2L, 2L)
    A = zeros(ComplexF64, 2L, 2L)
    B = zeros(ComplexF64, 2L, 2L)
    Avec = zeros(ComplexF64, Nobs)
    corr = zeros(ComplexF64, Ncorr)

    # Only the component of each observable along the local quantization axis
    # couples to the two-magnon channel, because S^z = s - b†b while the
    # transverse components are linear in b at leading order.
    pref = zeros(ComplexF64, Nobs, Na)

    # Masses of the binned pair-energy measure, ρs[iq][b] sitting at energy (b-1)*bin_width
    ρs = [eltype(measure)[] for _ in qpts.qs]

    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        ps = loop_wavevectors(grid, q_reshaped)
        for μ in 1:Nobs, i in 1:Na
            O = (data::SWTDataDipole).observables[μ, i]
            pref[μ, i] = conj(observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)) * O[3]
        end
        ρ = ρs[iq]

        for k_reshaped in ps
            # Columns L+1:2L of `A` create a quasi-particle at +k, and those of
            # `B` create one at q-k. Their momenta sum to the momentum transfer.
            dynamical_matrix!(H, swt, -k_reshaped)
            εA = bogoliubov!(A, H)
            dynamical_matrix!(H, swt, k_reshaped - q_reshaped)
            εB = bogoliubov!(B, H)

            for n₁ in 1:L, n₂ in 1:L
                ϵ = -εA[L+n₁] - εB[L+n₂]
                for μ in 1:Nobs
                    Avec[μ] = sum(1:Na) do i
                        # Symmetrized amplitude for creating the unordered pair
                        pref[μ, i] * conj(A[L+i, L+n₁]*B[i, L+n₂] + A[i, L+n₁]*B[L+i, L+n₂])
                    end
                end
                map!(corr, measure.corr_pairs) do (μ, ν)
                    # The factor 1/2 avoids double counting pair (1, 2) as (2, 1)
                    Avec[μ] * conj(Avec[ν]) / 2Ncells
                end
                val = measure.combiner(q_global, corr) / length(ps)
                (bin, f) = bin_index(ϵ, bin_width)
                while length(ρ) < bin + 1
                    push!(ρ, zero(eltype(measure)))
                end
                ρ[bin] += (1 - f) * val
                ρ[bin+1] += f * val
            end
        end
    end

    # Broadening of the binned measure, deferred to here so that the wavevector loop
    # above costs nothing per energy. The bins hold masses rather than densities, so
    # they are convolved with no measure of their own.
    nbins = maximum(length, ρs)
    bins = zeros(eltype(measure), nbins, length(qpts.qs))
    for (iq, ρ) in enumerate(ρs)
        view(bins, eachindex(ρ), iq) .= ρ
    end
    ret = zeros(eltype(measure), length(energies), length(qpts.qs))
    broaden!(ret, (0:nbins-1) * bin_width, bins; energies, kernel, Δω=1)

    ret = reshape(ret, length(energies), size(qpts.qs)...)
    return Intensities(cryst, qpts, energies, ret)
end
