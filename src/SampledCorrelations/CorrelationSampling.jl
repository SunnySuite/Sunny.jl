function observable_values!(buf, sys::System{N}, ops; apply_g = true) where N
    if N == 0
        for site in eachsite(sys), (i, op) in enumerate(ops)
            dipole = sys.dipoles[site]
            if apply_g
              dipole = sys.gs[site] * dipole
            end
            buf[i,site] = op * dipole
        end
    else
        Zs = sys.coherents
        #num_ops =  size(ops′, 3)
        #ops = reinterpret(SMatrix{N, N, ComplexF64, N*N}, reshape(ops′, N*N, num_ops))

        # SQTODO: This allocates :(
        for (i, op) in enumerate(ops)
          matrix_operator = convert(Matrix{ComplexF64},op)
            for site in eachsite(sys)
                buf[i,site] = dot(Zs[site], matrix_operator, Zs[site])
            end
        end
    end

    return nothing
end

function trajectory(sys::System{N}, Δt, nsnaps, ops; kwargs...) where N
    num_ops = length(ops)

    traj_buf = zeros(N == 0 ? Float64 : ComplexF64, num_ops, sys.latsize..., natoms(sys.crystal), nsnaps)
    trajectory!(traj_buf, sys, Δt, nsnaps, ops; kwargs...)

    return traj_buf
end

function trajectory!(buf, sys, Δt, nsnaps, ops; measperiod = 1, apply_g = true)
    @assert length(ops) == size(buf, 1)
    integrator = ImplicitMidpoint(Δt)

    observable_values!(@view(buf[:,:,:,:,:,1]), sys, ops; apply_g)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(sys, integrator)
        end
        observable_values!(@view(buf[:,:,:,:,:,n]), sys, ops; apply_g)
    end

    return nothing
end

function new_sample!(sc::SampledCorrelations, sys::System)
    (; Δt, samplebuf, measperiod, apply_g, observables, processtraj!) = sc
    nsnaps = size(samplebuf, 6)
    @assert size(sys.dipoles) == size(samplebuf)[2:5] "`System` size not compatible with given `SampledCorrelations`"

    trajectory!(samplebuf, sys, Δt, nsnaps, observables.observables; measperiod, apply_g)
    processtraj!(sc)

    return nothing
end

# At the sacrifice of code modularity, this processing step could be effected
# more efficiently by simply taking the real part of the trajectory after the
# Fourier transform
function symmetrize!(sc::SampledCorrelations)
    (; samplebuf) = sc
    nsteps = size(samplebuf, 6)
    mid = floor(Int, nsteps/2)
    for t in 1:mid, idx in CartesianIndices(size(samplebuf)[1:5])
        samplebuf[idx, t] = samplebuf[idx, nsteps-t+1] = 0.5*(samplebuf[idx, t] + samplebuf[idx, nsteps-t+1])
    end
end

function subtract_mean!(sc::SampledCorrelations)
    (; samplebuf) = sc
    nsteps = size(samplebuf, 6)
    meanvals = sum(samplebuf, dims=6) ./ nsteps
    samplebuf .-= meanvals
end

function no_processing(::SampledCorrelations)
    nothing
end

function accum_sample!(sc::SampledCorrelations;alg = :no_window)
    (; data, variance, observables, samplebuf, nsamples, fft!) = sc
    natoms = size(samplebuf,5)

    time_T = size(samplebuf,6)
    time_2T = 2time_T
    left_zero_ix = 1:time_T
    right_zero_ix = (time_T + 1):time_2T

    # Zero-padded extension of the samplebuf
    right_zero_buffer = zeros(ComplexF64,size(samplebuf)[1:5]...,time_2T)
    right_zero_buffer[:,:,:,:,:,left_zero_ix] .= samplebuf
    right_zero_buffer[:,:,:,:,:,right_zero_ix] .= 0

    # Number of terms contributing to each auto-correlation sum due to zero-padding
    statistical_power = reshape(time_T .- abs.(FFTW.fftfreq(time_2T,time_2T)),1,1,1,time_2T)
    statistical_power[statistical_power .== 0] .= Inf

    #fft! * samplebuf # Apply pre-planned and pre-normalized FFT
    #fft! * left_zero_buffer
    fft! * right_zero_buffer
    count = nsamples[1] += 1

    # Note that iterating over the `correlations` (a SortedDict) causes
    # allocations here. The contents of the loop contains no allocations. There
    # does not seem to be a big performance penalty associated with these
    # allocations.
    for j in 1:natoms, i in 1:natoms, (ci, c) in observables.correlations  
        α, β = ci.I

        # Convention is: S{α,β} means <α(t)> <β(0)>, i.e. α is delayed
        traj_α = @view right_zero_buffer[α,:,:,:,i,:]
        traj_β = @view right_zero_buffer[β,:,:,:,j,:]

        # FFT-accelerated cross correlation. Since both signals are zero-padded
        # in real time, all parts of the result contain meaningful correlations.
        correlation = FFTW.ifft(traj_α .* conj.(traj_β),4)

        correlation .*= 2 # Cancels the factor of two in the denominator of ifft

        correlation ./= statistical_power # Sam N.B.: I have this commented out locally

        if alg == :window
          correlation .*= reshape(cos.(range(0,π,length = time_2T)).^2,(1,1,1,time_2T))
        elseif alg == :chop
          correlation[:,:,:,1] .= 0
        end

        FFTW.fft!(correlation,4)

        # Remove overlapping highest frequency (only required if window doesn't already set this to zero)
        #correlation[:,:,:,longest_correlation_delay + 1] .= 0

        databuf = @view data[c,i,j,:,:,:,:]

        if isnothing(M)
            for k in eachindex(databuf)
                # Store the diff for one complex number on the stack.
                diff = correlation[k] - databuf[k]

                # Accumulate into running average
                databuf[k] += diff * (1/count)
            end
        else
            Mbuf = @view M[c,i,j,:,:,:,:]
            for k in eachindex(databuf)
                # Store old (complex) mean on stack.
                μ_old = databuf[k]

                # Update running mean.
                matrixelem = correlation[k]
                databuf[k] += (matrixelem - databuf[k]) * (1/count)
                μ = databuf[k]

                # Update variance estimate.
                # Note that the first term of `diff` is real by construction
                # (despite appearances), but `real` is explicitly called to
                # avoid automatic typecasting errors caused by roundoff.
                Mbuf[k] += real((matrixelem - μ_old)*conj(matrixelem - μ))
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
spin-spin correlations are then calculating in time and accumulated into `sc`. 

This function will change the state of `sys` when calculating dynamical
structure factor data. To preserve the initial state of `sys`, it must be saved
separately prior to calling `add_sample!`. Alternatively, the initial spin
configuration may be copied into a new `System` and this new `System` can be
passed to `add_sample!`.
"""
function add_sample!(sc::SampledCorrelations, sys::System; alg = :no_window, processtraj! = no_processing) 
    new_sample!(sc, sys; processtraj!)
    accum_sample!(sc;alg)
end
