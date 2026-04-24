"""
    SpinWaveTheoryKPM(sys::System; measure, regularization=1e-8, tol=nothing,
                      niters=nothing, method=:lanczos)

A variant of [`SpinWaveTheory`](@ref) that estimates [`intensities`](@ref) using
iterated matrix-vector products. By avoiding direct matrix diagonalization, this
method reduces computational cost from cubic to linear-scaling in the system
size ``N``. Large system sizes can arise, e.g., for models of quenched disorder
or models with nearly incommensurate ordering wavevectors.

Energy resolution is controlled by the dimensionless `tol` parameter. Common
choices are `tol=0.05` (more speed) or `0.01` (more accuracy). This will
determine the number of iterations as `M ≈ -2 log10(tol) Δϵ / fwhm`, where `Δϵ`
is the estimated spectral bandwidth of excitations and `fwhm` is the full width
at half maximum of the user-supplied broadening `kernel`. Computational cost
scales like ``𝒪(N M + M^2)``. Use `niters` instead of `tol` to directly specify
``M``.

!!! warning "Intensity loss at low-energy excitations"

    Not all numerical artifacts can be resolved by reducing `tol`. In particular,
    there may be unavoidable intensity loss at low-energy excitations, e.g., near
    Goldstone modes. This type of error originates from finite numerical precision
    and ill-conditioning of the dynamical matrix.

!!! tip "Consider `SampledCorrelations` when calculating powder averages"

    Spin wave theory requires an independent calculation for each ``𝐪`` point of
    interest. Consequently, it can be very slow to sample a 3D volume of
    ``𝐪``-space, e.g., as required for a [`powder_average`](@ref). A compelling
    alternative may be [`SampledCorrelations`](@ref). It uses real-time spin
    dynamics to calculate structure factor data over the entire 3D grid of
    commensurate ``𝐪``-vectors in one shot. This may provide a considerable speedup
    at the cost of: limited ``𝐪``-space resolution and stochastic error due to
    statistical sampling.

Two choices of `method` are possible. Lanczos is the default because it achieves
appears near-optimal accuracy at fixed iterations ``M`` [1, 2] and can detect
energetic instabilities. The alternative, `method=:kpm`, implements the Kernel
Polynomial Method as described in Ref. [3], which may be of historical interest.

## References

1. [T. Chen, _The Lanczos algorithm for matrix functions: a handbook for
   scientists_ (2024) [arXiv:2410.11090]](https://arxiv.org/abs/2410.11090).
2. [N. Amsel, T. Chen, A. Greenbaum, C. Musco, C. Musco, _Near-Optimal
   Approximation of Matrix Functions by the Lanczos Method_ (2023)
   [arXiv:2303.03358]](https://arxiv.org/abs/2303.03358).
3. [H. Lane et al., _Kernel Polynomial Method for Linear Spin Wave Theory_
   (2023) [arXiv:2312.08349]](https://arxiv.org/abs/2312.08349).
"""
struct SpinWaveTheoryKPM
    swt :: SpinWaveTheory
    tol :: Float64
    niters :: Int
    niters_bounds :: Int # Number of Lanczos iterations to bound spectrum
    method :: Symbol

    function SpinWaveTheoryKPM(sys::System; measure::Union{Nothing, MeasureSpec}, regularization=1e-8,
                               tol=nothing, niters=nothing, niters_bounds=10, method=:lanczos)
        xor(isnothing(tol), isnothing(niters)) || error("Exactly one of `tol` or `niters` must be specified.")
        method in (:lanczos, :kpm) || error("The method must be one of :lanczos or :kpm")
        tol = @something tol 1.0
        niters = @something niters 0
        0 < tol <= 1 || error("Require 0 < tol <= 1")
        niters >= 0 || error("Require niters >= 0")
        return new(SpinWaveTheory(sys; measure, regularization), tol, niters, niters_bounds, method)
    end
end

function Base.show(io::IO, ::MIME"text/plain", swt_kry::SpinWaveTheoryKPM)
    (; swt) = swt_kry
    printstyled(io, "SpinWaveTheoryKPM ", mode_to_str(swt.sys), "\n"; bold=true, color=:underline)
    println(io, "  ", natoms(swt.sys.crystal), " atoms")
end

function intensities(swt_kry::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, kT=0.0, verbose=false)
    qpts = convert(AbstractQPoints, qpts)
    data = zeros(eltype(swt_kry.swt.measure), length(energies), size(qpts.qs)...)
    return intensities!(data, swt_kry, qpts; energies, kernel, kT, verbose)
end

function intensities!(data, swt_kry::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, kT=0.0, verbose=false)
    qpts = convert(AbstractQPoints, qpts)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    (; method) = swt_kry
    if method == :lanczos
        intensities_lanczos!(data, swt_kry, qpts; energies, kernel, kT, verbose)
    else
        @assert method == :kpm
        intensities_kpm!(data, swt_kry, qpts; energies, kernel, kT, verbose)
    end
end


function mul_Ĩ!(y, x)
    L = size(y, 1) ÷ 2
    view(y, 1:L, :)    .= .+view(x, 1:L, :)
    view(y, L+1:2L, :) .= .-view(x, L+1:2L, :)
end

function mul_A!(swt, y, x, qs_reshaped, γ)
    L = size(y, 1) ÷ 2
    mul_dynamical_matrix!(swt, transpose(y), transpose(x), qs_reshaped)
    view(y, 1:L, :)    .*= +1/γ
    view(y, L+1:2L, :) .*= -1/γ
end

function set_moments!(moments, measure, u, α)
    map!(moments, measure.corr_pairs) do (μ, ν)
        dot(view(u, :, μ), view(α, :, ν))
    end
end

function intensities_kpm!(data, swt_kry, qpts; energies, kernel, kT, verbose)
    iszero(kT) || error("The :kpm backend does not support finite kT")

    (; swt, tol, niters, niters_bounds) = swt_kry
    (; sys, measure) = swt
    cryst = orig_crystal(sys)

    isnothing(kernel.fwhm) && error("Cannot determine the kernel fwhm")

    @assert eltype(data) == eltype(measure)
    @assert size(data) == (length(energies), size(qpts.qs)...)

    Na = nsites(sys)
    Ncells = Na / natoms(cryst)
    Nf = nflavors(swt)
    L = Nf*Na

    Nobs = size(measure.observables, 1)
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)
    moments = ElasticArray{ComplexF64}(undef, Ncorr, 0)

    u = zeros(ComplexF64, 2L, Nobs)
    α0 = zeros(ComplexF64, 2L, Nobs)
    α1 = zeros(ComplexF64, 2L, Nobs)
    α2 = zeros(ComplexF64, 2L, Nobs)

    for iq in CartesianIndices(qpts.qs)
        q = qpts.qs[iq]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q

        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)

        # Find extreme eigenvalues and rescaling factor
        lo, hi = eigbounds(swt, q_reshaped, niters_bounds)
        γ = 1.1 * max(abs(lo), hi)
        Δϵ = hi - lo

        # Determine order of polynomial expansion
        if niters > 0
            @assert tol == 1
            M = niters
        else
            @assert 0 < tol <= 1
            resolution = (kernel.fwhm/2) / (-log10(tol))
            M = max(round(Int, Δϵ/resolution), 2)
        end

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
            # Unlike Lanczos, KPM requires a lot of hacks here. Restrict to kT=0
            # to avoid the 1/x divergence of thermal_prefactor(x). To mitigate
            # ringing artifacts, introduce a smoothing energy scale σ. This is
            # the energy resolution (Δϵ / M) times a prefactor that grows slowly
            # with decreasing error tolerance, to help control artifacts.
            accuracy_prefactor = 2 * max(sqrt(-log10(tol)), 1)
            σ = accuracy_prefactor * (Δϵ / M)
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

function intensities_lanczos!(data, swt_kry, qpts; energies, kernel, kT, verbose)
    (; swt, tol, niters, niters_bounds) = swt_kry
    (; sys, measure) = swt
    cryst = orig_crystal(sys)

    isnothing(kernel.fwhm) && error("Cannot determine the kernel fwhm")

    @assert eltype(data) == eltype(measure)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    fill!(data, zero(eltype(data)))

    Na = nsites(sys)
    Ncells = Na / natoms(cryst)
    Nf = nflavors(swt)
    L = Nf*Na

    Nobs = size(measure.observables, 1)
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)

    u = zeros(ComplexF64, 2L, Nobs)
    v = zeros(ComplexF64, 2L)
    Sv = zeros(ComplexF64, 2L)

    for iq in CartesianIndices(qpts.qs)
        q = qpts.qs[iq]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q

        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)

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

        # Determine bounds on either the number of iterations or the resolution
        if niters > 0
            @assert tol == 1
            min_iters = niters
            resolution = Inf
        else
            @assert 0.0 < tol <= 1
            min_iters = niters_bounds
            resolution = (kernel.fwhm/2) / (-log10(tol))
        end

        for ξ in 1:Nobs
            # Don't accumulate observables that are zero
            iszero(view(u, :, ξ)) && continue

            mulA!(v, view(u, :, ξ))
            mulS!(Sv, v)
            c = sqrt(real(Sv' * v))
            v ./= c
            tridiag, lhs_adj_Q = try
                lanczos(mulA!, mulS!, v; lhs=u, min_iters, resolution, verbose)
            catch e
                if e.msg == "S is not a positive definite measure"
                    rethrow(ErrorException("Not an energy-minimum; wavevector q = $q unstable."))
                else
                    rethrow()
                end
            end

            (; values, vectors) = try
                eigen(tridiag)
            catch e
                # Fallback for https://github.com/JuliaLang/LinearAlgebra.jl/issues/1491
                eigen(Hermitian(collect(tridiag)))
            end

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
