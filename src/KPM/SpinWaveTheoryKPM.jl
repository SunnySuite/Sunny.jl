"""
    SpinWaveTheoryKPM(sys::System; measure, regularization=1e-8, tol)

A variant of [`SpinWaveTheory`](@ref) that uses the kernel polynomial method
(KPM) to calculate [`intensities`](@ref) [1]. This method avoids direct matrix
diagonalization, which scales cubically in the system size ``N``. Instead, using
``M`` iterative matrix-vector multiplications, the cost becomes ``𝒪(N M +
M^2)``. This provides a very significant acceleration when ``N`` is large.

The number of iterations is typically not too large: `M ≈ -2 log10(tol) Δϵ /
fwhm` where `Δϵ` is the estimated spectral bandwidth of excitations, `fwhm` is
the full width at half maximum of the broadening kernel, and `tol` is a
dimensionless tolerance parameter. Good choices are `0.05` (more speed) or
`0.01` (more accuracy).

## References

1. [H. Lane et al., _Kernel Polynomial Method for Linear Spin Wave Theory_
   (2023) [arXiv:2312.08349]](https://arxiv.org/abs/2312.08349).
"""
struct SpinWaveTheoryKPM
    swt :: SpinWaveTheory
    tol :: Float64
    lanczos_iters :: Int

    function SpinWaveTheoryKPM(sys::System; measure::Union{Nothing, MeasureSpec}, regularization=1e-8, tol, lanczos_iters=15)
        return new(SpinWaveTheory(sys; measure, regularization), tol, lanczos_iters)
    end
end


function mul_Ĩ!(y, x)
    L = size(y, 2) ÷ 2
    view(y, :, 1:L)    .= .+view(x, :, 1:L)
    view(y, :, L+1:2L) .= .-view(x, :, L+1:2L)
end

function mul_A!(swt, y, x, qs_reshaped, γ)
    L = size(y, 2) ÷ 2
    mul_dynamical_matrix!(swt, y, x, qs_reshaped)
    view(y, :, 1:L)    .*= +1/γ
    view(y, :, L+1:2L) .*= -1/γ
end

function set_moments!(moments, measure, u, α)
    map!(moments, measure.corr_pairs) do (μ, ν)
        dot(view(u, μ, :), view(α, ν, :))
    end
end


function intensities!(data, swt_kpm::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, kT=0.0, verbose=false)
    iszero(kT) || error("KPM does not yet support finite kT")
    qpts = convert(AbstractQPoints, qpts)

    (; swt, tol, lanczos_iters) = swt_kpm
    (; sys, measure) = swt
    cryst = orig_crystal(sys)

    isnothing(kernel.fwhm) && error("Cannot determine the kernel fwhm")

    @assert eltype(data) == eltype(measure)
    @assert size(data) == (length(energies), length(qpts.qs))

    Na = nsites(sys)
    Ncells = Na / natoms(cryst)
    Nf = nflavors(swt)
    L = Nf*Na
    Avec_pref = zeros(ComplexF64, Na) # initialize array of some prefactors

    Nobs = size(measure.observables, 1)
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)
    moments = ElasticArray{ComplexF64}(undef, Ncorr, 0)

    u = zeros(ComplexF64, Nobs, 2L)
    α0 = zeros(ComplexF64, Nobs, 2L)
    α1 = zeros(ComplexF64, Nobs, 2L)
    α2 = zeros(ComplexF64, Nobs, 2L)

    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q

        # Represent each local observable A(q) as a complex vector u(q) that
        # denotes a linear combination of HP bosons.

        for i in 1:Na
            r = sys.crystal.positions[i]
            ff = get_swt_formfactor(measure, 1, i)
            Avec_pref[i] = exp(2π*im * dot(q_reshaped, r))
            Avec_pref[i] *= compute_form_factor(ff, norm2(q_global))
        end

        if sys.mode == :SUN
            (; observables_localized) = swt.data::SWTDataSUN
            N = sys.Ns[1]
            for i in 1:Na, μ in 1:Nobs
                O = observables_localized[μ, i]
                for f in 1:Nf
                    u[μ, f + (i-1)*Nf]     = Avec_pref[i] * O[f, N]
                    u[μ, f + (i-1)*Nf + L] = Avec_pref[i] * O[N, f]
                end
            end
        else
            @assert sys.mode in (:dipole, :dipole_uncorrected)
            (; sqrtS, observables_localized) = swt.data::SWTDataDipole
            for i in 1:Na
                for μ in 1:Nobs
                    O = observables_localized[μ, i]
                    u[μ, i]   = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] + im*O[2])
                    u[μ, i+L] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] - im*O[2])
                end
            end
        end

        # Bound eigenvalue magnitudes and determine order of polynomial
        # expansion

        lo, hi = eigbounds(swt, q_reshaped, lanczos_iters)
        γ = 1.1 * max(abs(lo), hi)
        accuracy_factor = max(-3*log10(tol), 1)
        M = round(Int, accuracy_factor * max(2γ / kernel.fwhm, 3))
        resize!(moments, Ncorr, M)

        if verbose
            println("Bounds=", (lo, hi), " M=", M)
        end

        # Perform Chebyshev recursion

        q_repeated = fill(q_reshaped, Nobs)
        mul_Ĩ!(α0, u)
        mul_A!(swt, α1, α0, q_repeated, γ)
        set_moments!(view(moments, :, 1), measure, u, α0)
        set_moments!(view(moments, :, 2), measure, u, α1)
        for m in 3:M
            mul_A!(swt, α2, α1, q_repeated, γ)
            @. α2 = 2*α2 - α0
            set_moments!(view(moments, :, m), measure, u, α2)
            (α0, α1, α2) = (α1, α2, α0)
        end

        # Transform Chebyshev moments to intensities for each ω

        buf = zeros(2M)
        plan = FFTW.plan_r2r!(buf, FFTW.REDFT10)

        for (iω, ω) in enumerate(energies)
            # Ideally we would use thermal_prefactor instead of
            # thermal_prefactor_zero to allow for finite temperature effects.
            # Unfortunately, the Bose function's 1/x singularity introduces
            # divergence of the Chebyshev expansion integrals, and is tricky to
            # regularize. At kT=0, the occupation is a Heaviside step function.
            # To mitigate ringing artifacts associated with truncated Chebyshev
            # approximation, introduce smoothing on the energy scale σ. This is
            # the polynomial resolution scale times a prefactor that grows like
            # sqrt(accuracy) to reduce lingering ringing artifacts. See "AFM
            # KPM" for a test case where the smoothing degrades accuracy, and
            # "Disordered system with KPM" for an illustration of how smoothing
            # affects intensities at small ω.
            σ = sqrt(accuracy_factor) * (γ / M)
            thermal_prefactor_zero(x) = (tanh(x / σ) + 1) / 2
            f(x) = kernel(x, ω) * thermal_prefactor_zero(x)
            coefs = cheb_coefs!(M, f, (-γ, γ); buf, plan)
            # apply_jackson_kernel!(coefs)
            for i in 1:Ncorr
                corrbuf[i] = dot(coefs, view(moments, i, :)) / Ncells
            end
            data[iω, iq] = measure.combiner(q_global, corrbuf)
        end
    end

    return Intensities(cryst, qpts, collect(energies), data)
end

function intensities(swt_kpm::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, kT=0.0, verbose=false)
    qpts = convert(AbstractQPoints, qpts)
    data = zeros(eltype(swt_kpm.swt.measure), length(energies), length(qpts.qs))
    # return intensities!(data, swt_kpm, qpts; energies, kernel, kT, verbose)
    return intensities2!(data, swt_kpm, qpts; energies, kernel, kT, verbose)
end


function intensities2!(data, swt_kpm::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, kT=0.0, verbose=false)
    qpts = convert(AbstractQPoints, qpts)

    (; swt, tol, lanczos_iters) = swt_kpm
    (; sys, measure) = swt
    cryst = orig_crystal(sys)

    isnothing(kernel.fwhm) && error("Cannot determine the kernel fwhm")

    @assert eltype(data) == eltype(measure)
    @assert size(data) == (length(energies), length(qpts.qs))
    fill!(data, zero(eltype(data)))

    Na = nsites(sys)
    Ncells = Na / natoms(cryst)
    Nf = nflavors(swt)
    L = Nf*Na
    Avec_pref = zeros(ComplexF64, Na)

    Nobs = size(measure.observables, 1)
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)

    u = zeros(ComplexF64, 2L, Nobs)
    v = zeros(ComplexF64, 2L)
    Sv = zeros(ComplexF64, 2L)

    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q

        # Represent each local observable A(q) as a complex vector u(q) that
        # denotes a linear combination of HP bosons.

        for i in 1:Na
            r = sys.crystal.positions[i]
            ff = get_swt_formfactor(measure, 1, i)
            Avec_pref[i] = exp(2π*im * dot(q_reshaped, r))
            Avec_pref[i] *= compute_form_factor(ff, norm2(q_global))
        end

        if sys.mode == :SUN
            (; observables_localized) = swt.data::SWTDataSUN
            N = sys.Ns[1]
            for μ in 1:Nobs, i in 1:Na
                O = observables_localized[μ, i]
                for f in 1:Nf
                    u[f + (i-1)*Nf, μ]     = Avec_pref[i] * O[f, N]
                    u[f + (i-1)*Nf + L, μ] = Avec_pref[i] * O[N, f]
                end
            end
        else
            @assert sys.mode in (:dipole, :dipole_uncorrected)
            (; sqrtS, observables_localized) = swt.data::SWTDataDipole
            for μ in 1:Nobs, i in 1:Na
                O = observables_localized[μ, i]
                u[i, μ]   = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] + im*O[2])
                u[i+L, μ] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] - im*O[2])
            end
        end

        # Perform Lanczos calculation
 
        # w = Ĩ v
        function mulA!(w, v)
            @views w[1:L]    = +v[1:L]
            @views w[L+1:2L] = -v[L+1:2L]
            return w
        end

        # w = D v
        function mulS!(w, v)
            mul_dynamical_matrix!(swt, reshape(w, 1, :), reshape(v, 1, :), [q_reshaped])
            return w
        end

        resolution = (kernel.fwhm/2) / max(-log10(tol), 0)

        for ξ in 1:Nobs
            mulA!(v, view(u, :, ξ))
            mulS!(Sv, v)
            c = sqrt(Sv' * v)
            v ./= c
            tridiag, lhs_adj_Q = try
                lanczos(mulA!, mulS!, v; lhs=u, min_iters=lanczos_iters, resolution, verbose)
            catch e
                if e.msg == "S is not a positive definite measure"
                    rethrow(ErrorException("Not an energy-minimum; wavevector q = $q unstable."))
                else
                    rethrow()
                end
            end

            (; values, vectors) = eigen(tridiag)

            for (iω, ω) in enumerate(energies)
                f(x) = kernel(x, ω) * thermal_prefactor(x; kT)

                corr_ξ = c * lhs_adj_Q * vectors * Diagonal(f.(values)) * (vectors'[:, 1])

                # This step assumes that each local observable in the
                # correlation is Hermitian. In this case, bare correlations
                # should be symmetric, C[μ, ν] = C[ν, μ]*. The Lanczos
                # approximation C̃ breaks this symmetry. Restore it by looping
                # over ξ in 1:Nobs and accumulate Lanczos data C̃[:, ξ] in a
                # symmetric way. Accumulate C̃[μ, ξ] into C[μ, ν] if ξ = ν. Also
                # accumulate C̃[ν, ξ]* into C[μ, ν] if ξ = μ. A factor of 1/2
                # avoids double counting. In the special case that μ = ν, this
                # assigns real(C̃[μ, μ]) to C[μ, μ] only once.
                corrbuf .= 0
                for (i, (μ, ν)) in enumerate(measure.corr_pairs)
                    ξ == ν && (corrbuf[i] += (1/2) *     (corr_ξ[μ] / Ncells))
                    ξ == μ && (corrbuf[i] += (1/2) * conj(corr_ξ[ν] / Ncells))
                end

                # This step assumes that combiner is linear, so that it is valid
                # to move the ξ loop outside the data accumulation. One could
                # relax this assumption by preallocating an array of size (Nω,
                # Ncorr) to accumulate into corrbuf prior to calling combiner.
                data[iω, iq] += measure.combiner(q_global, corrbuf)
            end
        end
    end

    return Intensities(cryst, qpts, collect(energies), data)
end
