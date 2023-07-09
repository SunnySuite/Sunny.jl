function observable_values!(buf, sys::System{N}, ops; apply_g = true) where N
    if N == 0
        for site in all_sites(sys), (i, op) in enumerate(ops)
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
          matrix_operator = convert(Matrix,op)
            for site in all_sites(sys)
                buf[i,site] = dot(Zs[site], matrix_operator, Zs[site])
            end
        end
    end

    return nothing
end

function trajectory(sys::System{N}, integrator, nsnaps, ops; kwargs...) where N
    num_ops = length(ops)

    traj_buf = zeros(N == 0 ? Float64 : ComplexF64, num_ops, sys.latsize..., natoms(sys.crystal), nsnaps)
    trajectory!(traj_buf, sys, integrator, nsnaps, ops; kwargs...)

    return traj_buf
end

function trajectory!(buf, sys, integrator, nsnaps, ops; measperiod = 1, apply_g = true)
    @assert length(ops) == size(buf, 1)

    observable_values!(@view(buf[:,:,:,:,:,1]), sys, ops, apply_g = apply_g)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(sys, integrator)
        end
        observable_values!(@view(buf[:,:,:,:,:,n]), sys, ops, apply_g = apply_g)
    end

    return nothing
end

function new_sample!(sf::StructureFactor, sys::System; processtraj! = no_processing)
    (; integrator, samplebuf, measperiod, apply_g) = sf
    nsnaps = size(samplebuf, 6)

    @assert size(sys.dipoles) == size(samplebuf)[2:5] "`System` size not compatible with given `StructureFactor`"

    trajectory!(samplebuf, sys, integrator, nsnaps, sf.observables; measperiod = measperiod, apply_g = apply_g)

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

function accum_sample!(sf::StructureFactor)
    (; data, correlations, samplebuf, copybuf, nsamples, fft!) = sf
    natoms = size(samplebuf)[5]

    fft! * samplebuf # Apply pre-planned and pre-normalized FFT
    nsamples[1] += 1

    # Transfer to final memory layout while accumulating new samplebuf
    for j in 1:natoms, i in 1:natoms, (ci, c) in correlations 
        α, β = ci.I
        @. copybuf = @views samplebuf[α,:,:,:,i,:] * conj(samplebuf[β,:,:,:,j,:]) - data[c,i,j,:,:,:,:] 
        @views data[c,i,j,:,:,:,:] .+= copybuf .* (1/nsamples[1])
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
