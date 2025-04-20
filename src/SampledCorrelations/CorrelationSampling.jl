function observable_values!(buf, sys::System{N}, observables, atom_idcs) where N
    if N == 0
        for i in axes(observables, 1)
            for site in eachsite(sys)
                obs = observables[i, site]
                dipole = sys.dipoles[site]
                buf[i, site] = obs ⋅ dipole
            end
        end
    else
        Zs = sys.coherents
        for idx in CartesianIndices(observables)
            _, la, lb, lc, pos = idx.I
            atom = atom_idcs[la, lb, lc, pos]
            buf[idx] = dot(Zs[la, lb, lc, atom], observables[idx], Zs[la, lb, lc, atom])
        end
    end
    return nothing
end

function trajectory!(buf, sys, integrator, nsnaps, observables, atom_idcs; measperiod=1)
    @assert size(observables, 1) == size(buf, 1)
    observable_values!(@view(buf[:,:,:,:,:,1]), sys, observables, atom_idcs)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(sys, integrator)
        end
        observable_values!(@view(buf[:,:,:,:,:,n]), sys, observables, atom_idcs)
    end
    return nothing
end

function new_sample!(sc::SampledCorrelations, sys::System)
    (; integrator, samplebuf, measperiod, observables, atom_idcs) = sc

    # Only fill the sample buffer half way; the rest is zero-padding
    buf_size = size(samplebuf, 6)
    nsnaps = (buf_size÷2) + 1
    samplebuf[:,:,:,:,:,(nsnaps+1):end] .= 0

    # @assert size(sys.dipoles) == size(samplebuf)[2:5] "`System` size not compatible with given `SampledCorrelations`"

    trajectory!(samplebuf, sys, integrator, nsnaps, observables, atom_idcs; measperiod)

    return nothing
end

function accum_sample!(sc::SampledCorrelations; window)
    (; data, M, corr_pairs, samplebuf, corrbuf, space_fft!, time_fft!, corr_fft!, corr_ifft!) = sc
    npos = size(samplebuf)[5]
    num_time_offsets = size(samplebuf, 6)
    T = (num_time_offsets÷2) + 1 # Duration that each signal was recorded for

    # Time offsets (in samples) Δt = [0,1,...,(T-1),-(T-1),...,-1] produced by 
    # the cross-correlation between two length-T signals
    time_offsets = FFTW.fftfreq(num_time_offsets, num_time_offsets)

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
    n_contrib = reshape(T .- abs.(time_offsets), 1, 1, 1, num_time_offsets)
    
    # As long as `num_time_offsets` is odd, there will be a non-zero number of
    # contributions, so we don't need this line
    #$ @assert isodd(num_time_offsets)
    n_contrib[n_contrib .== 0] .= Inf

    count = sc.nsamples += 1

    for j in 1:npos, i in 1:npos, (c, (α, β)) in enumerate(corr_pairs)
        # α, β = ci.I

        sample_α = @view samplebuf[α,:,:,:,i,:]
        sample_β = @view samplebuf[β,:,:,:,j,:]
        databuf  = @view data[c,i,j,:,:,:,:]

        # According to Sunny convention, the correlation is between
        # α† and β. This conjugation implements both the dagger on the α
        # as well as the appropriate spacetime offsets of the correlation.
        @. corrbuf = conj(sample_α) * sample_β
        corr_ifft! * corrbuf
        corrbuf ./= n_contrib

        @assert window in (:cosine, :rectangular)
        if window == :cosine
            # Multiply the real-time correlation data by a cosine window that
            # smoothly goes to zero at offsets approaching the trajectory
            # length, Δt → T. This smooth windowing mitigates ringing artifacts
            # that appear when imposing periodicity on the real-space
            # trajectory. Note, however, that windowing also broadens the signal
            # S(ω) on the characteristic scale of one frequency bin Δω = 2π/T.
            window_func = cos.(range(0, π, length=num_time_offsets+1)[1:end-1]).^2
            corrbuf .*= reshape(window_func, 1, 1, 1, num_time_offsets)
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
    add_sample!(sc::SampledCorrelationsStatic, sys::System)

Measure pair correlation data for the spin configuration in `sys`, and
accumulate these statistics into `sc`. For a dynamical
[`SampledCorrelations`](@ref), this involves time-integration of the provided
spin trajectory, recording correlations in both space and time. Conversely,
[`SampledCorrelationsStatic`](@ref), will record only spatial correlations for
the single spin configuration that is provided.

Time-integration will update the spin configuration of `sys` in-place. To avoid
this mutation, consider calling [`clone_system`](@ref) prior to `add_sample!`.
"""
function add_sample!(sc::SampledCorrelations, sys::System; window=:cosine)
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
    # The hidden option `window=:rectangular` will disable smooth windowing.
    # This may be of interest for extracting real-time dynamical correlations.

    new_sample!(sc, sys)
    accum_sample!(sc; window)
end

function add_sample!(sc::SampledCorrelationsStatic, sys::System; window=:cosine)
    add_sample!(sc.parent, sys; window)
end
