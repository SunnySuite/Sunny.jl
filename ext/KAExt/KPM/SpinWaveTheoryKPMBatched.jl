using LinearAlgebra: SymTridiagonal, Diagonal, Hermitian, eigen, dot

struct SpinWaveTheoryKPMBatchedDeviceKA{T<:AbstractFloat, TSwtD, TIncoming}
    swt_d         :: TSwtD
    swt_host      :: Sunny.SpinWaveTheory
    incoming      :: TIncoming
    tol           :: Float64
    niters        :: Int
    niters_bounds :: Int
    method        :: Symbol
    precision     :: Type{T}
end

function Adapt.adapt_structure(to, swt_kry_d::SpinWaveTheoryKPMBatchedDeviceKA{T}) where T
    swt_d_a    = Adapt.adapt_structure(to, swt_kry_d.swt_d)
    incoming_a = Adapt.adapt_structure(to, swt_kry_d.incoming)
    SpinWaveTheoryKPMBatchedDeviceKA{T, typeof(swt_d_a), typeof(incoming_a)}(
        swt_d_a,
        swt_kry_d.swt_host,
        incoming_a,
        swt_kry_d.tol,
        swt_kry_d.niters,
        swt_kry_d.niters_bounds,
        swt_kry_d.method,
        T,
    )
end

function Sunny.to_device_batched(swt_kry::Sunny.SpinWaveTheoryKPM,
                                 backend::KernelAbstractions.Backend;
                                 precision::Type{<:AbstractFloat}=Float64)
    swt_kry.method == :lanczos || error(
        "Unsupported method=$(swt_kry.method). Only :lanczos is supported on GPU."
    )

    precision in (Float32, Float64) || error(
        "Only Float32 and Float64 are supported for `precision`; got $precision."
    )

    sys = swt_kry.swt.sys
    if sys.mode in (:dipole, :dipole_uncorrected)
        swt_d    = Sunny.to_device(swt_kry.swt, backend, precision)
        incoming = build_incoming_bond_csr_ka(swt_kry.swt.sys, backend, precision)
    elseif sys.mode == :SUN
        swt_d = Sunny.to_device(swt_kry.swt, backend, precision)

        if precision === Float32
            incoming = build_empty_incoming_csr(backend)
        else
            incoming = build_incoming_bond_csr_sun_ka(swt_kry.swt.sys, backend)
        end
    else
        error("Unsupported sys.mode=$(sys.mode).")
    end

    return SpinWaveTheoryKPMBatchedDeviceKA{precision, typeof(swt_d), typeof(incoming)}(
        swt_d,
        swt_kry.swt,
        incoming,
        swt_kry.tol,
        swt_kry.niters,
        swt_kry.niters_bounds,
        swt_kry.method,
        precision,
    )
end

function Sunny.intensities(swt_kry_d::SpinWaveTheoryKPMBatchedDeviceKA, qpts;
                           energies, kernel::Sunny.AbstractBroadening,
                           kT=0.0, verbose=false)
    qpts = convert(Sunny.AbstractQPoints, qpts)
    measure = swt_kry_d.swt_host.measure
    data = zeros(eltype(measure), length(energies), size(qpts.qs)...)
    return Sunny.intensities!(data, swt_kry_d, qpts;
                              energies, kernel, kT, verbose)
end

function Sunny.intensities!(data, swt_kry_d::SpinWaveTheoryKPMBatchedDeviceKA{T}, qpts;
                            energies, kernel::Sunny.AbstractBroadening,
                            kT=0.0, verbose=false) where T
    qpts = convert(Sunny.AbstractQPoints, qpts)
    @assert size(data) == (length(energies), size(qpts.qs)...)

    if T === Float64
        return intensities_lanczos_device_batched_ka!(data, swt_kry_d, qpts;
                                                      energies, kernel, kT, verbose)
    else
        return intensities_lanczos_device_batched_fp32_ka!(data, swt_kry_d, qpts;
                                                           energies, kernel, kT, verbose)
    end
end

function intensities_lanczos_device_batched_ka!(data, swt_kry_d::SpinWaveTheoryKPMBatchedDeviceKA{T}, qpts;
                                                energies, kernel, kT, verbose) where T
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
    Nq      = length(qpts.qs)
    N_chains = Nq * Nobs

    tol           = swt_kry_d.tol
    niters_user   = swt_kry_d.niters
    niters_bounds = swt_kry_d.niters_bounds

    CT = Complex{T}
    backend = KernelAbstractions.get_backend(swt_kry_d.swt_d.sys.extfield)

    resolution = Inf
    if niters_user > 0
        @assert tol == 1
        max_iters_global = min(niters_user, twoL)
        max_iters_global += mod(max_iters_global, 2)
    else
        @assert 0.0 < tol <= 1
        q_idx_rep = max(1, length(qpts.qs) ÷ 2)
        q_rep = qpts.qs[q_idx_rep]
        q_rep_reshaped = Sunny.to_reshaped_rlu(sys, q_rep)
        lo, hi = Sunny.eigbounds(swt, q_rep_reshaped, niters_bounds)
        Δϵ = hi - lo
        resolution = (kernel.fwhm/2) / (-log10(tol))
        max_iters_unbounded = max(niters_bounds, fld(Δϵ, resolution))
        max_iters_with_safety = ceil(Int, 4.0 * max_iters_unbounded)
        max_iters_with_safety += mod(max_iters_with_safety, 2)
        max_iters_global = min(max_iters_with_safety, twoL)
        if verbose
            println("Eigbounds at q_idx=$q_idx_rep: Δϵ=$Δϵ")
            println("max_iters_global=$max_iters_global (4× safety on single-q estimate, buffer upper bound only)")
            println("Per-chain niters_eff determined adaptively at iteration $niters_bounds")
        end
    end
    @assert max_iters_global >= 1 "max_iters_global must be ≥ 1"

    q_reshaped_host = Sunny.Vec3[
        Sunny.to_reshaped_rlu(sys, q) for q in vec(qpts.qs)
    ]
    q_reshaped_chain_host = Vector{SVector{3, T}}(undef, N_chains)
    @inbounds for iq in 1:Nq
        for ξ in 1:Nobs
            c = (iq - 1) * Nobs + ξ
            q_reshaped_chain_host[c] = SVector{3, T}(q_reshaped_host[iq])
        end
    end
    q_reshaped_chain_dev = _backend_array(backend, q_reshaped_chain_host)

    v_batch   = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    vp_batch  = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    Sv_batch  = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    Svp_batch = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    w_batch   = KernelAbstractions.zeros(backend, CT, N_chains, twoL)
    Sw_batch  = KernelAbstractions.zeros(backend, CT, N_chains, twoL)

    lhs_dev   = KernelAbstractions.zeros(backend, CT,        twoL, Nobs, Nq)
    corr_dev  = KernelAbstractions.zeros(backend, ComplexF64, Nobs, Nobs, Nq)
    corr_host = Array{ComplexF64}(undef, Nobs, Nobs, Nq)

    α_chains_dev   = KernelAbstractions.zeros(backend, ComplexF64, N_chains)
    β²_chains_dev  = KernelAbstractions.zeros(backend, ComplexF64, N_chains)
    β_chains_dev   = KernelAbstractions.zeros(backend, T,          N_chains)
    α_chains_host  = Vector{ComplexF64}(undef, N_chains)
    β²_chains_host = Vector{ComplexF64}(undef, N_chains)

    lhs_adj_Q_buf = Array{ComplexF64}(undef, Nobs, max_iters_global, N_chains)

    u_host_per_chain = zeros(ComplexF64, twoL, N_chains)
    Avec_pref = zeros(ComplexF64, Na)

    if sys.mode == :SUN
        let
            (; observables_localized) = swt.data::Sunny.SWTDataSUN
            N_hilbert = sys.Ns[1]
            for iq in 1:Nq
                q        = qpts.qs[iq]
                q_r      = q_reshaped_host[iq]
                q_global = cryst.recipvecs * q
                for i in 1:Na
                    r  = sys.crystal.positions[i]
                    ff = Sunny.get_swt_formfactor(measure, 1, i)
                    Avec_pref[i]  = exp(2π*im * dot(q_r, r))
                    Avec_pref[i] *= Sunny.compute_form_factor(ff, Sunny.norm2(q_global))
                end
                for ξ in 1:Nobs
                    c = (iq - 1) * Nobs + ξ
                    for i in 1:Na
                        O = observables_localized[ξ, i]
                        for f in 1:Nf
                            l = f + (i - 1) * Nf
                            u_host_per_chain[l,   c] = Avec_pref[i] * O[f, N_hilbert]
                            u_host_per_chain[l+L, c] = Avec_pref[i] * O[N_hilbert, f]
                        end
                    end
                end
            end
        end
    else
        let
            (; sqrtS, observables_localized) = swt.data::Sunny.SWTDataDipole
            for iq in 1:Nq
                q          = qpts.qs[iq]
                q_reshaped = q_reshaped_host[iq]
                q_global   = cryst.recipvecs * q
                for i in 1:Na
                    r  = sys.crystal.positions[i]
                    ff = Sunny.get_swt_formfactor(measure, 1, i)
                    Avec_pref[i]  = exp(2π*im * dot(q_reshaped, r))
                    Avec_pref[i] *= Sunny.compute_form_factor(ff, Sunny.norm2(q_global))
                end
                for ξ in 1:Nobs
                    c = (iq - 1) * Nobs + ξ
                    for i in 1:Na
                        O = observables_localized[ξ, i]
                        u_host_per_chain[i,   c] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] + im*O[2])
                        u_host_per_chain[i+L, c] = Avec_pref[i] * (sqrtS[i] / √2) * (O[1] - im*O[2])
                    end
                end
            end
        end
    end

    chain_skip = falses(N_chains)
    @inbounds for c in 1:N_chains
        if iszero(view(u_host_per_chain, :, c))
            chain_skip[c] = true
        end
    end

    copyto!(reshape(lhs_dev, twoL, N_chains), convert(Array{CT}, u_host_per_chain))

    u_host_dummy = copy(u_host_per_chain)
    @inbounds for c in 1:N_chains
        if chain_skip[c]
            u_host_dummy[1, c] = ComplexF64(1, 0)
        end
    end

    v_batch_host = permutedims(u_host_dummy, (2, 1))
    @views v_batch_host[:, L+1:twoL] .*= -1
    copyto!(v_batch, convert(Array{CT}, v_batch_host))

    mulA! = let L = L, twoL = twoL
        (out, in_) -> begin
            @views @. out[:, 1:L]      = +in_[:, 1:L]
            @views @. out[:, L+1:twoL] = -in_[:, L+1:twoL]
            return out
        end
    end

    mulS! = let swt_d = swt_kry_d.swt_d,
                incoming = swt_kry_d.incoming,
                q_chain = q_reshaped_chain_dev
        (out, in_) -> begin
            mul_dynamical_matrix_batched_ka!(swt_d, incoming, out, in_, q_chain)
            return out
        end
    end

    αs_matrix, βs_matrix, lhs_adj_Q, n_iters_done_per_chain, c_norm_per_chain = try
        lanczos_device_batched_ka(
            mulA!, mulS!,
            v_batch, vp_batch, Sv_batch, Svp_batch, w_batch, Sw_batch,
            lhs_dev, corr_dev, corr_host,
            α_chains_dev, β²_chains_dev, β_chains_dev,
            α_chains_host, β²_chains_host,
            lhs_adj_Q_buf;
            max_iters = max_iters_global,
            min_iters = niters_bounds,
            resolution = resolution,
            verbose = verbose,
        )
    catch e
        if e isa ErrorException && occursin("not a positive definite measure", e.msg)
            rethrow(ErrorException("GPU batched Lanczos: not an energy-minimum; some q wavevector unstable. " *
                                   "Original error: $(e.msg)"))
        else
            rethrow()
        end
    end

    Threads.@threads for iq_linear in 1:Nq
        iq = CartesianIndices(qpts.qs)[iq_linear]
        q = qpts.qs[iq]
        q_global = cryst.recipvecs * q

        corrbuf = zeros(ComplexF64, Ncorr)

        for ξ in 1:Nobs
            c = (iq_linear - 1) * Nobs + ξ
            chain_skip[c] && continue
            c_norm_per_chain[c] == 0 && continue

            n = n_iters_done_per_chain[c]

            αs_chain = view(αs_matrix, 1:n, c)
            βs_chain = n > 1 ? view(βs_matrix, 1:n-1, c) : Float64[]

            tridiag = SymTridiagonal(collect(αs_chain), collect(βs_chain))

            (; values, vectors) = try
                eigen(tridiag)
            catch
                eigen(Hermitian(collect(tridiag)))
            end

            c_norm = c_norm_per_chain[c]

            lhs_adj_Q_chain = view(lhs_adj_Q, :, 1:n, c)

            for (iω, ω) in enumerate(energies)
                f(x) = kernel(x, ω) * Sunny.thermal_prefactor(x; kT)

                corr_ξ = c_norm * lhs_adj_Q_chain * vectors * Diagonal(f.(values)) * (vectors'[:, 1])

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
