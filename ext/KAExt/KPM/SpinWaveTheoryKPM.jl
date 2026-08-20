using LinearAlgebra: SymTridiagonal, Diagonal, Hermitian, eigen, dot

struct SpinWaveTheoryKPMDeviceKA{TSwtD}
    swt_d         :: TSwtD
    swt_host      :: Sunny.SpinWaveTheory
    tol           :: Float64
    niters        :: Int
    niters_bounds :: Int
    method        :: Symbol
end

function Adapt.adapt_structure(to, swt_kry_d::SpinWaveTheoryKPMDeviceKA)
    swt_d_a = Adapt.adapt_structure(to, swt_kry_d.swt_d)
    SpinWaveTheoryKPMDeviceKA(
        swt_d_a,
        swt_kry_d.swt_host,
        swt_kry_d.tol,
        swt_kry_d.niters,
        swt_kry_d.niters_bounds,
        swt_kry_d.method,
    )
end

function Sunny.to_device(swt_kry::Sunny.SpinWaveTheoryKPM, backend::KernelAbstractions.Backend)
    swt_kry.method == :lanczos || error(
        "GPU device path only supports `method=:lanczos`. " *
        "Got method=$(swt_kry.method). The `:kpm` Chebyshev method is " *
        "not supported in the GPU device path."
    )
    swt_d = Sunny.to_device(swt_kry.swt, backend)
    return SpinWaveTheoryKPMDeviceKA(
        swt_d,
        swt_kry.swt,
        swt_kry.tol,
        swt_kry.niters,
        swt_kry.niters_bounds,
        swt_kry.method,
    )
end

function Sunny.intensities(swt_kry_d::SpinWaveTheoryKPMDeviceKA, qpts;
                           energies, kernel::Sunny.AbstractBroadening,
                           kT=0.0, verbose=false)
    qpts = convert(Sunny.AbstractQPoints, qpts)
    measure = swt_kry_d.swt_host.measure
    data = zeros(eltype(measure), length(energies), size(qpts.qs)...)
    return Sunny.intensities!(data, swt_kry_d, qpts;
                              energies, kernel, kT, verbose)
end

function Sunny.intensities!(data, swt_kry_d::SpinWaveTheoryKPMDeviceKA, qpts;
                            energies, kernel::Sunny.AbstractBroadening,
                            kT=0.0, verbose=false)
    qpts = convert(Sunny.AbstractQPoints, qpts)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    return intensities_lanczos_device_ka!(data, swt_kry_d, qpts;
                                          energies, kernel, kT, verbose)
end

function intensities_lanczos_device_ka!(data, swt_kry_d::SpinWaveTheoryKPMDeviceKA, qpts;
                                        energies, kernel, kT, verbose)
    swt = swt_kry_d.swt_host
    sys = swt.sys
    measure = swt.measure
    cryst = Sunny.orig_crystal(sys)

    isnothing(kernel.fwhm) && error("Cannot determine the kernel fwhm")

    @assert eltype(data) == eltype(measure)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    fill!(data, zero(eltype(data)))

    Na      = Sunny.nsites(sys)
    Ncells  = Na / Sunny.natoms(cryst)
    Nf      = Sunny.nflavors(swt)
    L       = Nf * Na
    twoL    = 2 * L

    Nobs    = size(measure.observables, 1)
    Ncorr   = length(measure.corr_pairs)

    tol           = swt_kry_d.tol
    niters_user   = swt_kry_d.niters
    niters_bounds = swt_kry_d.niters_bounds

    backend = KernelAbstractions.get_backend(swt_kry_d.swt_d.sys.extfield)

    q_reshaped_host = Sunny.Vec3[
        Sunny.to_reshaped_rlu(sys, q) for q in vec(qpts.qs)
    ]
    q_reshaped_all_d = _backend_array(backend, q_reshaped_host)

    v   = KernelAbstractions.zeros(backend, ComplexF64, twoL)
    vp  = KernelAbstractions.zeros(backend, ComplexF64, twoL)
    Sv  = KernelAbstractions.zeros(backend, ComplexF64, twoL)
    Svp = KernelAbstractions.zeros(backend, ComplexF64, twoL)
    w   = KernelAbstractions.zeros(backend, ComplexF64, twoL)
    Sw  = KernelAbstractions.zeros(backend, ComplexF64, twoL)

    lhs_dev = KernelAbstractions.zeros(backend, ComplexF64, twoL, Nobs)

    corr_dev  = KernelAbstractions.zeros(backend, ComplexF64, Nobs)
    corr_host = Vector{ComplexF64}(undef, Nobs)

    max_iters_global = max(twoL, 1)
    lhs_adj_Q_buf = Matrix{ComplexF64}(undef, Nobs, max_iters_global)

    Avec_pref = zeros(ComplexF64, Na)
    u_host    = zeros(ComplexF64, twoL, Nobs)
    corrbuf   = zeros(ComplexF64, Ncorr)

    for (chain_idx, iq) in enumerate(CartesianIndices(qpts.qs))
        q          = qpts.qs[iq]
        q_reshaped = q_reshaped_host[chain_idx]
        q_global   = cryst.recipvecs * q

        q_d_view = view(q_reshaped_all_d, chain_idx:chain_idx)

        for i in 1:Na
            r = sys.crystal.positions[i]
            ff = Sunny.get_swt_formfactor(measure, 1, i)
            Avec_pref[i] = exp(2π*im * dot(q_reshaped, r))
            Avec_pref[i] *= Sunny.compute_form_factor(ff, Sunny.norm2(q_global))
        end

        @assert sys.mode in (:dipole, :dipole_uncorrected) "GPU device path is dipole-only"
        let
            (; sqrtS, observables_localized) = swt.data::Sunny.SWTDataDipole
            for μ in 1:Nobs, i in 1:Na
                O = observables_localized[μ, i]
                u_host[i,   μ] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] + im*O[2])
                u_host[i+L, μ] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] - im*O[2])
            end
        end

        copyto!(lhs_dev, u_host)

        if niters_user > 0
            @assert tol == 1
            min_iters_chain = niters_user
            resolution_chain = Inf
        else
            @assert 0.0 < tol <= 1
            min_iters_chain = niters_bounds
            resolution_chain = (kernel.fwhm/2) / (-log10(tol))
        end
        max_iters_chain = max_iters_global

        for ξ in 1:Nobs
            iszero(view(u_host, :, ξ)) && continue

            mulA! = let L = L, twoL = twoL
                (w, v) -> begin
                    @views @. w[1:L]      = +v[1:L]
                    @views @. w[L+1:twoL] = -v[L+1:twoL]
                    return w
                end
            end

            mulS! = let swt_d = swt_kry_d.swt_d, q_d_view = q_d_view
                (Sw, w) -> begin
                    Sunny.mul_dynamical_matrix!(swt_d,
                        reshape(Sw, 1, :), reshape(w, 1, :), q_d_view)
                    return Sw
                end
            end

            @views mulA!(v, lhs_dev[:, ξ])

            mulS!(Sv, v)
            c = sqrt(real(dot(Sv, v)))
            @. v /= c

            tridiag, lhs_adj_Q = try
                lanczos_device_ka(
                    mulA!, mulS!,
                    v, vp, Sv, Svp, w, Sw,
                    lhs_dev, corr_dev, corr_host, lhs_adj_Q_buf;
                    min_iters = min_iters_chain,
                    max_iters = max_iters_chain,
                    resolution = resolution_chain,
                    verbose = verbose,
                )
            catch e
                if e isa ErrorException && occursin("not a positive definite measure", e.msg)
                    rethrow(ErrorException("Not an energy-minimum; wavevector q = $q unstable."))
                else
                    rethrow()
                end
            end

            (; values, vectors) = try
                eigen(tridiag)
            catch e
                eigen(Hermitian(collect(tridiag)))
            end

            for (iω, ω) in enumerate(energies)
                f(x) = kernel(x, ω) * Sunny.thermal_prefactor(x; kT)

                corr_ξ = c * lhs_adj_Q * vectors * Diagonal(f.(values)) * (vectors'[:, 1])

                corrbuf .= 0
                for (i, (μ, ν)) in enumerate(measure.corr_pairs)
                    ξ == ν && (corrbuf[i] += (1/2) *     (corr_ξ[μ] / Ncells))
                    ξ == μ && (corrbuf[i] += (1/2) * conj(corr_ξ[ν] / Ncells))
                end

                data[iω, iq] += measure.combiner(q_global, corrbuf)
            end
        end
    end

    return Sunny.Intensities(cryst, qpts, collect(energies), data)
end
