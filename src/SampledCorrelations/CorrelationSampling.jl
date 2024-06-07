function observable_values!(buf, sys::System{N}, ops) where N
    if N == 0
        for (i, op) in enumerate(ops)
            for site in eachsite(sys)
                A = observable_at_site(op,site)
                dipole = sys.dipoles[site]
                buf[i,site] = A * dipole
            end
        end
    else
        Zs = sys.coherents
        for (i, op) in enumerate(ops)
            for site in eachsite(sys)
                A = observable_at_site(op,site)
                buf[i,site] = dot(Zs[site], A, Zs[site])
            end
        end
    end

    return nothing
end

function trajectory(sys::System{N}, dt, nsnaps, ops; kwargs...) where N
    num_ops = length(ops)

    traj_buf = zeros(N == 0 ? Float64 : ComplexF64, num_ops, sys.latsize..., natoms(sys.crystal), nsnaps)
    trajectory!(traj_buf, sys, dt, nsnaps, ops; kwargs...)

    return traj_buf
end

function trajectory!(buf, sys, dt, nsnaps, ops; measperiod = 1)
    @assert length(ops) == size(buf, 1)
    integrator = ImplicitMidpoint(dt)

    observable_values!(@view(buf[:,:,:,:,:,1]), sys, ops)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(sys, integrator)
        end
        observable_values!(@view(buf[:,:,:,:,:,n]), sys, ops)
    end

    return nothing
end

function new_sample!(sc::SampledCorrelations, sys::System)
    (; dt, samplebuf, measperiod, observables, processtraj!) = sc

    # Only fill the sample buffer half way; the rest is zero-padding
    buf_size = size(samplebuf, 6)
    nsnaps = (buf_size÷2) + 1
    samplebuf[:,:,:,:,:,(nsnaps+1):end] .= 0

    @assert size(sys.dipoles) == size(samplebuf)[2:5] "`System` size not compatible with given `SampledCorrelations`"

    trajectory!(samplebuf, sys, dt, nsnaps, observables.observables; measperiod)
    processtraj!(sc)

    return nothing
end

function no_processing(::SampledCorrelations)
    nothing
end

function accum_sample!(sc::SampledCorrelations; window)
    (; data, M, observables, samplebuf, corrbuf, nsamples, space_fft!, time_fft!, corr_fft!, corr_ifft!) = sc
    natoms = size(samplebuf)[5]

    num_time_offsets = size(samplebuf,6)
    T = (num_time_offsets÷2) + 1 # Duration that each signal was recorded for

    # Time offsets (in samples) Δt = [0,1,...,(T-1),-(T-1),...,-1] produced by 
    # the cross-correlation between two length-T signals
    time_offsets = FFTW.fftfreq(num_time_offsets,num_time_offsets)

    # Transform A(q) = ∑ exp(iqr) A(r).
    # This is opposite to the FFTW convention, so we must conjugate
    # the fft by a complex conjugation to get the correct sign.
    samplebuf .= conj.(samplebuf)
    space_fft! * samplebuf
    samplebuf .= conj.(samplebuf)

    # Transform A(ω) = ∑ exp(-iωt) A(t)
    # In samplebuf, the original signal is from 1:T, and the rest
    # is zero-padding, from (T+1):num_time_offsets. This allows a
    # usual FFT in the time direction, even though the signal isn't periodic.
    time_fft! * samplebuf

    # Number of contributions to the DFT sum (non-constant due to zero-padding).
    # Equivalently, this is the number of estimates of the correlation with
    # each offset Δt that need to be averaged over.
    n_contrib = reshape(T .- abs.(time_offsets),1,1,1,num_time_offsets)
    
    # As long as `num_time_offsets` is odd, there will be a non-zero number of
    # contributions, so we don't need this line
    @assert isodd(num_time_offsets)
    #n_contrib[n_contrib .== 0] .= Inf

    count = nsamples[1] += 1

    # Note that iterating over the `correlations` (a SortedDict) causes
    # allocations here. The contents of the loop contains no allocations. There
    # does not seem to be a big performance penalty associated with these
    # allocations.
    for j in 1:natoms, i in 1:natoms, (ci, c) in observables.correlations  
        α, β = ci.I

        sample_α = @view samplebuf[α,:,:,:,i,:]
        sample_β = @view samplebuf[β,:,:,:,j,:]
        databuf  = @view data[c,i,j,:,:,:,:]

        # According to Sunny convention, the correlation is between
        # α† and β. This conjugation implements both the dagger on the α
        # as well as the appropriate spacetime offsets of the correlation.
        @. corrbuf = conj(sample_α) * sample_β
        corr_ifft! * corrbuf
        corrbuf ./= n_contrib

        if window == :cosine
          # Apply a cosine windowing to force the correlation at Δt=±(T-1) to be zero
          # to force periodicity. In terms of the spectrum S(ω), this applys a smoothing
          # with a characteristic lengthscale of O(1) frequency bins.
          window_func = cos.(range(0,π,length = num_time_offsets + 1)[1:end-1]).^2
          corrbuf .*= reshape(window_func,1,1,1,num_time_offsets)
        end

        corr_fft! * corrbuf

        if isnothing(M)
            for k in eachindex(databuf)
                # Store the diff for one complex number on the stack.
                diff = corrbuf[k] - databuf[k]

                # Accumulate into running average
                databuf[k] += diff * (1/count)
            end
        else
            Mbuf = @view M[c,i,j,:,:,:,:]
            for k in eachindex(databuf)
                # Store old (complex) mean on stack.
                μ_old = databuf[k]

                # Update running mean.
                databuf[k] += (corrbuf[k] - databuf[k]) / count
                μ = databuf[k]

                # Update variance estimate.
                # Note that the first term of `diff` is real by construction
                # (despite appearances), but `real` is explicitly called to
                # avoid automatic typecasting errors caused by roundoff.
                Mbuf[k] += real((corrbuf[k] - μ_old)*conj(corrbuf[k] - μ))
            end
        end
    end

    return nothing
end


"""
    add_sample!(sc::SampledCorrelations, sys::System)

`add_trajectory` uses the spin configuration contained in the `System` to
generate a correlation data and accumulate it into `sc`. For static structure
factors, this involves analyzing the spin-spin correlations of the spin
configuration provided. For a dynamic structure factor, a trajectory is
calculated using the given spin configuration as an initial condition. The
spin-spin correlations are then calculated in time and accumulated into `sc`. 

This function will change the state of `sys` when calculating dynamical
structure factor data. To preserve the initial state of `sys`, it must be saved
separately prior to calling `add_sample!`. Alternatively, the initial spin
configuration may be copied into a new `System` and this new `System` can be
passed to `add_sample!`.
"""
function add_sample!(sc::SampledCorrelations, sys::System; window = :cosine)
    # Sunny now estimates the dynamical structure factor in two steps. First, it
    # estimates real-time correlations C(t) = ⟨S(t)S(0)⟩ from classical
    # dynamics. Second, it takes the Fourier transform of C(t) to get the
    # structure factor in energy space ω. Because time-correlations are
    # estimated from dynamical trajectories of finite duration T, there is no
    # data for C(t) when |t| > T. The inverse trajectory length sets a limit on
    # the energy resolution in the structure factor (dω ≳ 1/T). In practice,
    # rather than imposing a sharp cutoff on C(t) (a rectangular window), it is
    # favorable to multiply C(t) by a window function that goes smoothly to zero
    # as |t| -> T. For concreteness, Sunny selects a cosine window, but any
    # other smooth window would likely work similarly well in mitigating
    # artifacts. See https://github.com/SunnySuite/Sunny.jl/pull/246 for more
    # discussion about the calculation of intensities from classical dynamics. 
    # 
    # The `window` parameter to this function is *not* part of Sunny's public
    # API, and is subject to change at any time. Passing an alternative value
    # for `window` will replace the smooth cosine window with a rectangular
    # window (sharp truncation at |t| = T). This experimental feature is
    # provided so that expert users have the ability to extract real-time
    # dynamical correlations. In the future, a public API will be designed to
    # give more direct access to the real-time correlations.

    new_sample!(sc, sys)
    accum_sample!(sc; window)
end
