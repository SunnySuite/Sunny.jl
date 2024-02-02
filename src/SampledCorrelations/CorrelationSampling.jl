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
                #matrix_operator = convert(Matrix{ComplexF64},A)
                buf[i,site] = dot(Zs[site], A, Zs[site])
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

function trajectory!(buf, sys, Δt, nsnaps, ops; measperiod = 1)
    @assert length(ops) == size(buf, 1)
    integrator = ImplicitMidpoint(Δt)

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
    (; Δt, samplebuf, measperiod, observables, processtraj!) = sc
    nsnaps = size(samplebuf, 6)
    @assert size(sys.dipoles) == size(samplebuf)[2:5] "`System` size not compatible with given `SampledCorrelations`"

    trajectory!(samplebuf, sys, Δt, nsnaps, observables.observables; measperiod)
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

function accum_sample!(sc::SampledCorrelations)
    (; data, M, observables, samplebuf, nsamples, fft!) = sc
    natoms = size(samplebuf)[5]

    fft! * samplebuf # Apply pre-planned and pre-normalized FFT
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

        if isnothing(M)
            for k in eachindex(databuf)
                # Store the diff for one complex number on the stack.
                diff = sample_α[k] * conj(sample_β[k]) - databuf[k]

                # Accumulate into running average
                databuf[k] += diff * (1/count)
            end
        else
            Mbuf = @view M[c,i,j,:,:,:,:]
            for k in eachindex(databuf)
                # Store old (complex) mean on stack.
                μ_old = databuf[k]

                # Update running mean.
                matrixelem = sample_α[k] * conj(sample_β[k])
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
function add_sample!(sc::SampledCorrelations, sys::System) 
    new_sample!(sc, sys)
    accum_sample!(sc)
end
