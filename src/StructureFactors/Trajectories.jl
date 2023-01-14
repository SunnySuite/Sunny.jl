expectation(Z, op) = Z' * (op * Z)  # No need for real -- keeping buffer complex and FFTing in place

function observable_expectations!(buf, sys::SpinSystem{N}, ops′::Array{ComplexF64, 3}) where N
    Zs = sys.coherents
    num_ops =  size(ops′, 3)
    ops = reinterpret(SMatrix{N, N, ComplexF64, N*N}, reshape(ops′, N*N, num_ops))

    for (i, op) in enumerate(ops)
        for idx in CartesianIndices(Zs)
            κ = sys.κs[idx[4]]
            buf[i,idx] = κ * expectation(Zs[idx], op) 
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


function compute_mag!(M, sys::SpinSystem, gfactor = true)
    if gfactor
        for idx in CartesianIndices(sys.dipoles)
            g = sys.gs[idx[4]]
            M[:, idx] .= g * sys.dipoles[idx]
        end
    else
        for idx in CartesianIndices(sys.dipoles)
            M[:, idx] .= sys.dipoles[idx]
        end
    end
end

function dipole_trajectory(sys, integrator, nsnaps; kwargs...)
    traj_buf = zeros(ComplexF64, 3, size(sys.dipoles)..., nsnaps)
    dipole_trajectory!(traj_buf, sys, integrator, nsnaps; kwargs...)
    return traj_buf
end

function dipole_trajectory!(buf, sys, integrator, nsnaps; measperiod = 1, gfactor = true)
    @assert size(buf, 1) == 3

    compute_mag!(@view(buf[:,:,:,:,:,1]), sys, gfactor)
    for n in 2:nsnaps
        for _ in 1:measperiod
            step!(sys, integrator)
        end
        compute_mag!(@view(buf[:,:,:,:,:,n]), sys, gfactor)
    end

    return nothing
end

function new_trajectory!(sftraj::SFTrajectory, sys_original::SpinSystem)
    (; dipolemode, integrator, traj, measperiod, sys, gfactor) = sftraj
    nsnaps = size(traj, 6)
    sys.dipoles .= sys_original.dipoles
    sys.coherents .= sys_original.coherents

    if dipolemode
        dipole_trajectory!(traj, sys, integrator, nsnaps; measperiod, gfactor)
    else
        expectation_trajectory!(traj, sys, integrator, nsnaps, sftraj.ops; measperiod)
    end

    return nothing
end
new_trajectory!(sf::StructureFactor, sys::SpinSystem) = new_trajectory!(sf.sftraj, sys)


# ddtodo: Plan FFT
function accum_trajectory!(sfdata::SFData, sftraj::SFTrajectory, nsamples::Int64)
    (; data, idxinfo) = sfdata
    (; traj, sys) = sftraj
    nb, nω = size(traj)[5:6]

    FFTW.fft!(traj, (2,3,4,6))
    traj /= nω * √(prod(sys.latsize))  # Normalize FFT according to physical convention

    ## Transfer to final memory layout while accumulating new trajectory
    for cell in CartesianIndices(sys.latsize), b2 in 1:nb, b1 in 1:nb, ω in 1:nω
        for (ci, c) in idxinfo 
            α, β = ci.I
            old = data[c,b1,b2,cell,ω]
            data[c,b1,b2,cell,ω] =  old + (traj[α,cell,b1,ω] * conj(traj[β,cell,b2,ω]) - old) / nsamples
        end
    end

    return nothing
end

function accum_trajectory!(sf::StructureFactor)
    sf.nsamples += 1
    accum_trajectory!(sf.sfdata, sf.sftraj, sf.nsamples)
    return nothing
end

function add_trajectory!(sf::StructureFactor, sys)
    new_trajectory!(sf, sys)
    accum_trajectory!(sf)
end