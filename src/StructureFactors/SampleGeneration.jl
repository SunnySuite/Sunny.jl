function observable_expectations!(buf, sys::System{N}, ops′::Array{ComplexF64, 3}) where N
    Zs = sys.coherents
    num_ops =  size(ops′, 3)
    ops = reinterpret(SMatrix{N, N, ComplexF64, N*N}, reshape(ops′, N*N, num_ops))

    for (i, op) in enumerate(ops)
        for site in all_sites(sys)
            buf[i,site] = dot(Zs[site], op, Zs[site]) 
        end
    end

    return nothing
end

function expectation_trajectory(sys, integrator, nsnaps, ops; kwargs...)
    num_ops =  size(ops, 3)

    traj_buf = zeros(ComplexF64, num_ops, size(sys.coherents)..., nsnaps)
    expectation_trajectory!(traj_buf, sys, integrator, nsnaps, ops; kwargs...)

    return traj_buf
end

function expectation_trajectory!(buf, sys, integrator, nsnaps, ops; measperiod = 1)
    @assert size(ops, 3) == size(buf, 1)

    observable_expectations!(@view(buf[:,:,:,:,:,1]), sys, ops)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(sys, integrator)
        end
        observable_expectations!(@view(buf[:,:,:,:,:,n]), sys, ops)
    end

    return nothing
end


function compute_mag!(M, sys::System, apply_g = true)
    for site in all_sites(sys)
        if apply_g
            M[:, site] .= sys.gs[site] * sys.dipoles[site]
        else
            M[:, site] .= sys.dipoles[site]
        end
    end
end

function dipole_trajectory(sys, integrator, nsnaps; kwargs...)
    traj_buf = zeros(ComplexF64, 3, size(sys.dipoles)..., nsnaps)
    dipole_trajectory!(traj_buf, sys, integrator, nsnaps; kwargs...)
    return traj_buf
end

function dipole_trajectory!(buf, sys, integrator, nsnaps; measperiod = 1, apply_g = true)
    @assert size(buf, 1) == 3

    compute_mag!(@view(buf[:,:,:,:,:,1]), sys, apply_g)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(sys, integrator)
        end
        compute_mag!(@view(buf[:,:,:,:,:,n]), sys, apply_g)
    end

    return nothing
end

function new_sample!(sf::StructureFactor, sys::System; processtraj! = no_processing)
    (; integrator, samplebuf, measperiod, apply_g) = sf
    nsnaps = size(samplebuf, 6)

    @assert size(sys.dipoles) == size(samplebuf)[2:5] "`System` size not compatible with given `StructureFactor`"

    if size(sf.observables)[1] == 0
        dipole_trajectory!(samplebuf, sys, integrator, nsnaps; measperiod, apply_g)
    else
        expectation_trajectory!(samplebuf, sys, integrator, nsnaps, sf.observables; measperiod)
    end

    processtraj!(sf)

    return nothing
end

function symmetrize!(sf::StructureFactor)
    (; samplebuf) = sf
    nsteps = size(samplebuf, 6)
    for t in 1:nsteps
        selectdim(samplebuf, 6, t) .= 0.5*(selectdim(samplebuf, 6, t) + selectdim(samplebuf, 6, nsteps-t+1))
    end
end

function subtract_mean!(sf::StructureFactor)
    (; samplebuf) = sf
    nsteps = size(samplebuf, 6)
    meanvals = sum(samplebuf, dims=6) ./ nsteps
    samplebuf .-= meanvals
end

function no_processing(::StructureFactor)
    nothing
end

# ddtodo: Plan FFT
function accum_sample!(sf::StructureFactor)
    (; data, idxinfo, samplebuf, nsamples) = sf
    latsize = size(samplebuf)[2:4]
    na, nω = size(samplebuf)[5:6]

    FFTW.fft!(samplebuf, (2,3,4,6))
    samplebuf /= nω * √(prod(latsize))  # Normalize FFT according to physical convention
    nsamples[1] += 1

    # Transfer to final memory layout while accumulating new samplebuf
    for ω in 1:nω, cell in CartesianIndices(latsize), j in 1:na, i in 1:na
        for (ci, c) in idxinfo 
            α, β = ci.I
            old = data[c,i,j,cell,ω]
            data[c,i,j,cell,ω] = old + (samplebuf[α,cell,i,ω] * conj(samplebuf[β,cell,j,ω]) - old) / nsamples[1]
        end
    end

    return nothing
end


"""
    add_sample!(sf::StructureFactor, sys::System)

`add_trajectory` uses the spin configuration contained in the `System` to
generate a correlation data and accumulate it into `sf`. For static structure
factors, this involves analyzing the spin-spin correlations of the spin
configuration provided. For a dynamic structure factor, a trajectory is
calculated using the given spin configuration as an initial condition. The
spin-spin correlations are then calculating in time and accumulated into `sf`. 

This function will change the state of `sys` when calculating dynamical
structure factor data. To preserve the initial state of `sys`, it must be saved
separately prior to calling `add_sample!`. Alternatively, the initial spin
configuration may be copied into a new `System` and this new `System` can be
passed to `add_sample!`.
"""
function add_sample!(sf::StructureFactor, sys::System; processtraj! = no_processing) 
    new_sample!(sf, sys; processtraj!)
    accum_sample!(sf)
end