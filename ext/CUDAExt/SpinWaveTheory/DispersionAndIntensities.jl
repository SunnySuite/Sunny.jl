using LinearAlgebra

function _set_identity(a)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if i > size(a, 2)
        return
    end

    iq = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    if iq > size(a, 3)
        return
    end

    L = div(size(a, 1), 2)
    a[i,i,iq] = i > L ? -1. : 1.
    return
end

function _frequencies(H, evalues)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if i > size(H, 2)
        return
    end

    j = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    if j > size(H, 3)
        return
    end

    λ = evalues[i, j]
    # Normalize columns of T so that para-unitarity holds, T† Ĩ T = Ĩ.
    view(H,:,i,j) .*= CUDA.rsqrt(abs(λ))

    # Inverse of λ are eigenvalues of Ĩ H, or equivalently, of √H Ĩ √H.
    evalues[i, j] = inv(λ)
    return
end

function _intensities(swt, qs, L, Ncells, H, Nobs, Na, Ncorr, recipvecs, intensity, kT, disp)
    iq = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if iq > size(H,3)
        return
    end
    q = qs[iq]
    Hq = view(H,:,:,iq)
    corrbuf = CuDynamicSharedArray(ComplexF64, (Ncorr, blockDim().x))
    corrbufq = view(corrbuf,:,threadIdx().x)
    Avec_pref = CuDynamicSharedArray(ComplexF64, (Nobs, Na, blockDim().x), (Ncorr + Nobs) * blockDim().x * sizeof(ComplexF64))
    Avec_prefq = view(Avec_pref,:,:,threadIdx().x)

    (; sys, data, measure) = swt
    q_global = recipvecs * q
    for i in 1:Na, μ in 1:Nobs
        r_global = global_position(sys, (1, 1, 1, i)) # + offsets[μ, i]
        ff = Sunny.get_swt_formfactor(measure, μ, i)
        Avec_prefq[μ, i] = exp(- im * dot(q_global, r_global))
        Avec_prefq[μ, i] *= Sunny.compute_form_factor(ff, Sunny.norm2(q_global))
    end

    Avec = CuDynamicSharedArray(ComplexF64, (Nobs, blockDim().x), Ncorr * blockDim().x * sizeof(ComplexF64))
    Avecq = view(Avec, :, threadIdx().x)
    # Fill `intensity` array
    for band in 1:L
        for idx in eachindex(Avecq)
            Avecq[idx] = 0.
        end
        t = view(Hq, :, band+L)
        for i in 1:Na, μ in 1:Nobs
            O = data.observables_localized[μ, i]
            # This is the Avec of the two transverse and one
            # longitudinal directions in the local frame. (In the
            # local frame, z is longitudinal, and we are computing
            # the transverse part only, so the last entry is zero)
            displacement_local_frame = Sunny.SA[t[i + Na] + t[i], im * (t[i + Na] - t[i]), 0.0]
            Avecq[μ] += Avec_prefq[μ, i] * (data.sqrtS[i]/√2) * (O' * displacement_local_frame)[1]
        end
        for idx in eachindex(corrbufq)
            (μ, ν) = measure.corr_pairs[idx]
            corrbufq[idx] = Avecq[μ] * conj(Avecq[ν]) / Ncells
        end

        intensity[band, iq] = Sunny.thermal_prefactor(disp[band, iq]; kT) * measure.combiner(q_global, corrbufq)
    end
    return
end

"""
    intensities_bands(swt::SpinWaveTheory, qpts; kT=0)

Calculate spin wave excitation bands for a set of q-points in reciprocal space.
This calculation is analogous to [`intensities`](@ref), but does not perform
line broadening of the bands.
"""
function intensities_bands(swt::SpinWaveTheoryDevice, qpts; kT=0, with_negative=false)
    (; sys, measure) = swt
    isempty(measure.observables) && error("No observables! Construct SpinWaveTheory with a `measure` argument.")
    with_negative && error("Option `with_negative=true` not yet supported.")
    #@assert sys.mode in (:dipole, :dipole_uncorrected)
    #@assert isnothing(sys.ewald)

    qpts = convert(Sunny.AbstractQPoints, qpts)

    cryst = sys.original_crystal

    # Number of (magnetic) atoms in magnetic cell
    @assert sys.dims == (1,1,1)
    Na = Sunny.nsites(sys)
    # Number of chemical cells in magnetic cell
    Ncells = Na / Sunny.natoms(cryst)
    # Number of quasiparticle modes
    L = Sunny.nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)
    
    # Temporary storage for pair correlations
    Nobs = Sunny.num_observables(measure)
    Ncorr = Sunny.num_correlations(measure)

    qs_d = qpts.qs

    # Given q in reciprocal lattice units (RLU) for the original crystal, return a
    # q_reshaped in RLU for the possibly-reshaped crystal.
    reshaped_rlu = inv(2π) * sys.crystal.latvecs' * cryst.recipvecs
    H_d = CUDA.zeros(ComplexF64, 2L, 2L, Nq)
    Sunny.dynamical_matrix!(H_d, swt, reshaped_rlu, qs_d)

    H_dp = [view(H_d,:,:,i) for i in 1:Nq]
    CUSOLVER.potrfBatched!('L', H_dp)

    I_d = CUDA.zeros(ComplexF64, 2L, 2L, Nq)
    kernel = @cuda launch=false _set_identity(I_d)
    config = launch_configuration(kernel.fun)
    optimal_threads_1d = config.threads

    threads_x = Base.min(2L, optimal_threads_1d)
    threads_y = Base.min(Nq, optimal_threads_1d ÷ threads_x)
    threads = (threads_x, threads_y) # e.g., (16, 32) or similar, max product is 1024

    blocks_x = cld(2L, threads_x)
    blocks_y = cld(Nq, threads_y)
    blocks = (blocks_x, blocks_y)

    kernel(I_d; threads=threads, blocks=blocks)

    I_dp = [view(I_d,:,:,i) for i in 1:Nq]
    CUBLAS.trsm_batched!('R', 'L', 'C', 'N', ComplexF64(1.), H_dp, I_dp)
    CUBLAS.trsm_batched!('L', 'L', 'N', 'N', ComplexF64(1.), H_dp, I_dp)
    #evalues_d , _ = CUSOLVER.heevjBatched!('V', 'L', I_d)
    evalues_d , _ = CUSOLVER.XsyevBatched!('V', 'L', I_d)
    CUBLAS.trsm_batched!('L', 'L', 'C', 'N', ComplexF64(1.), H_dp, I_dp)

    kernel = @cuda launch=false _frequencies(I_d, evalues_d)
    config = launch_configuration(kernel.fun)
    optimal_threads_1d = config.threads

    threads_x = Base.min(2L, optimal_threads_1d)
    threads_y = optimal_threads_1d ÷ threads_x
    threads_y = Base.min(Nq, threads_y)
    threads = (threads_x, threads_y) # e.g., (16, 32) or similar, max product is 1024

    blocks_x = cld(2L, threads_x)
    blocks_y = cld(Nq, threads_y)
    blocks = (blocks_x, blocks_y)

    kernel(I_d, evalues_d; threads=threads, blocks=blocks)
    disp_d = view(evalues_d,L+1:2L, :)

    intensity_d = CUDA.zeros(eltype(measure), L, Nq)
    kernel = @cuda launch=false _intensities(swt, qs_d, L, Ncells, I_d, Nobs, Na, Ncorr, cryst.recipvecs, intensity_d, kT, disp_d)
    get_shmem(threads; Nobs=Nobs, Na=Na, Ncorr=Ncorr) = threads * sizeof(ComplexF64) * (Nobs * (1 + Na) + Ncorr)
    config = launch_configuration(kernel.fun, shmem=threads->get_shmem(threads))
    threads = Base.min(Nq, config.threads)
    blocks = cld(Nq, threads)
    kernel(swt, qs_d, L, Ncells, I_d, Nobs, Na, Ncorr, cryst.recipvecs, intensity_d, kT, disp_d; threads=threads, blocks=blocks, shmem=get_shmem(threads))

    disp_d = reshape(CuArray(disp_d), L, size(qpts.qs)...)
    intensity_d = reshape(intensity_d, L, size(qpts.qs)...)
    return BandIntensitiesDevice(cryst, qpts, disp_d, intensity_d)
end
