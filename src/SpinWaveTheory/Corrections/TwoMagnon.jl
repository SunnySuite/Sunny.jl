"""
    intensities_two_magnon(swt::SpinWaveTheory, qpts; energies, kernel, opts...)

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
integral runs over the first magnetic Brillouin zone.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration. See the
[HCubature](https://github.com/JuliaMath/HCubature.jl) documentation for
details. Because the integrand becomes singular as the `kernel` width shrinks,
a very narrow `kernel` will be expensive.
"""
function intensities_two_magnon(swt::SpinWaveTheory, qpts; energies, kernel::AbstractBroadening, opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")
    check_corrections_supported(swt)

    (; sys, measure, data) = swt
    num_observables(measure) == 0 && error("No observables! Construct SpinWaveTheory with a `measure` argument.")
    energies = collect(Float64, energies)
    issorted(energies) || error("energies must be sorted")

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

    ret = zeros(eltype(measure), length(energies), length(qpts.qs))

    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        for μ in 1:Nobs, i in 1:Na
            O = (data::SWTDataDipole).observables[μ, i]
            pref[μ, i] = conj(observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)) * O[3]
        end

        acc = hcubature((0, 0, 0), (1, 1, 1); opts...) do k_reshaped
            # Columns L+1:2L of `A` create a quasi-particle at +k, and those of
            # `B` create one at q-k. Their momenta sum to the momentum transfer.
            dynamical_matrix!(H, swt, -Vec3(k_reshaped))
            εA = bogoliubov!(A, H)
            dynamical_matrix!(H, swt, Vec3(k_reshaped) - q_reshaped)
            εB = bogoliubov!(B, H)

            out = zeros(eltype(measure), length(energies))
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
                val = measure.combiner(q_global, corr)
                for (iω, ω) in enumerate(energies)
                    out[iω] += kernel(ϵ, ω) * val
                end
            end
            return out
        end

        # Error bars in acc[2] are discarded
        view(ret, :, iq) .= acc[1]
    end

    ret = reshape(ret, length(energies), size(qpts.qs)...)
    return Intensities(cryst, qpts, energies, ret)
end
