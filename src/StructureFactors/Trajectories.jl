expectation(Z, op) = Z' * (op * Z)  # No need for real -- keeping buffer complex and FFTing in place

function observable_expectations!(buf, sys::SpinSystem{N}, ops′::Array{ComplexF64, 3}) where N
    Zs = sys.coherents
    num_ops =  size(ops′, 3)
    ops = reinterpret(SMatrix{N, N, ComplexF64, N*N}, reshape(ops′, N*N, num_ops))

    for (i, op) in enumerate(ops)
        for site in 1:nbasis(sys.crystal)
            κ = sys.site_infos[site].spin_rescaling
            for idx in CartesianIndices(sys.size)
                buf[i,idx,site] = κ * expectation(Zs[idx,site], op) 
            end
        end
    end

    return nothing
end

function expectation_trajectory(sys, integrator, num_snaps, ops; kwargs...)
    num_ops =  size(ops, 3)

    traj_buf = zeros(ComplexF64, num_ops, size(sys)..., num_snaps)
    expectation_trajectory!(traj_buf, sys, integrator, num_snaps, ops; kwargs...)

    return traj_buf
end

function expectation_trajectory!(buf, sys, integrator, num_snaps, ops; meas_period = 1)
    @assert size(ops, 3) == size(buf, 1)

    observable_expectations!(@view(buf[:,:,:,:,:,1]), sys, ops)
    for n in 2:num_snaps
        for _ in 1:meas_period
            step!(sys, integrator)
        end
        observable_expectations!(@view(buf[:,:,:,:,:,n]), sys, ops)
    end

    return nothing
end


function compute_mag!(M, sys::SpinSystem, g_factor = true)
    if g_factor
        for b in 1:nbasis(sys.crystal)
            gS = sys.site_infos[b].g 
            for idx in CartesianIndices(sys.size)
                M[:, idx, b] .= gS * sys.dipoles[idx, b]
            end
        end
    else
        for b in 1:nbasis(sys.crystal), idx in CartesianIndices(sys.size)
            M[:, idx, b] .= sys.dipoles[idx, b]
        end
    end
end

function dipole_trajectory(sys, integrator, num_snaps; kwargs...)
    traj_buf = zeros(ComplexF64, 3, size(sys)..., num_snaps)
    dipole_trajectory!(traj_buf, sys, integrator, num_snaps; kwargs...)
    return traj_buf
end

function dipole_trajectory!(buf, sys, integrator, num_snaps; meas_period = 1, g_factor = true)
    @assert size(buf, 1) == 3

    compute_mag!(@view(buf[:,:,:,:,:,1]), sys, g_factor)
    for n in 2:num_snaps
        for _ in 1:meas_period
            step!(sys, integrator)
        end
        compute_mag!(@view(buf[:,:,:,:,:,n]), sys, g_factor)
    end

    return nothing
end

function new_trajectory!(sftraj::SFTrajectory, sys_original::SpinSystem)
    (; dipolemode, integrator, traj, meas_period, sys) = sftraj
    num_snaps = size(traj, 6)
    sys.dipoles .= sys_original.dipoles
    sys.coherents .= sys_original.coherents

    if dipolemode
        dipole_trajectory!(traj, sys, integrator, num_snaps; meas_period, g_factor = sftraj.g_factor)
    else
        expectation_trajectory!(traj, sys, integrator, num_snaps, sftraj.ops; meas_period)
    end

    return nothing
end
new_trajectory!(sf::StructureFactor, sys::SpinSystem) = new_trajectory!(sf.sftraj, sys)


# TODO: Plan FFTs (put in SFTrajectory)
function accum_trajectory!(sfdata::SFData, sftraj::SFTrajectory, num_samples::Int64)
    (; data, idx_info) = sfdata
    (; traj) = sftraj
    nops, na, nb, nc, ns, nω = size(traj)

    FFTW.fft!(traj, (2,3,4,6))
    traj /= nω * √(na*nb*nc)  # Normalize FFT according to physical convention
    Sα = reshape(traj, (nops, na, nb, nc, 1, ns, nω)) 
    Sβ = reshape(traj, (nops, na, nb, nc, ns, 1, nω)) 

    for ((α, β), idx) in idx_info
        # Test if broadcasting eliminates need for views
        res = @view(data[idx,:,:,:,:,:,:])
        α = @view(Sα[α,:,:,:,:,:,:])
        β = @view(Sβ[β,:,:,:,:,:,:])
        @. res = res + (α * conj(β) - res) / num_samples
    end

    return nothing
end

function accum_trajectory!(sf::StructureFactor)
    sf.num_samples += 1
    accum_trajectory!(sf.sfdata, sf.sftraj, sf.num_samples)
    return nothing
end

function add_trajectory!(sf::StructureFactor, sys)
    new_trajectory!(sf, sys)
    accum_trajectory!(sf)
end